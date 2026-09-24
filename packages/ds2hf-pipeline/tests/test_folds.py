from __future__ import annotations

import random
import re
from collections import Counter
from pathlib import Path

import pytest
from ds2hf_pipeline.tcga import folds, verify

# Project sizes chosen to straddle both fold counts: one below 10 patients.
_SIZES = {"TCGA-AAA": 51, "TCGA-BBB": 23, "TCGA-CCC": 7}


def _cohort(seed: int = 0) -> list[dict]:
    rng = random.Random(seed)
    rows = []
    for project, n in _SIZES.items():
        for i in range(n):
            rows.append(
                {
                    "case_id": f"{project}-{i:04d}-{rng.random():.12f}",
                    "case_submitter_id": f"TCGA-{project[-2:]}-{i:04d}",
                    "project_id": project,
                    "pfi_status": rng.choice(folds.STATUSES),
                    "os_status": rng.choice(folds.STATUSES),
                    "status_source": rng.choice([folds.LIU, folds.REDERIVED]),
                }
            )
    return rows


def _spread(rows: list[dict], fold: str, k: int, group) -> int:
    counts = Counter((group(r), r[fold]) for r in rows)
    groups = {group(r) for r in rows}
    return max(
        max(counts[(g, f)] for f in range(k)) - min(counts[(g, f)] for f in range(k))
        for g in groups
    )


def test_deal_ignores_input_order() -> None:
    cohort = _cohort()
    shuffled = cohort[:]
    random.Random(1).shuffle(shuffled)
    assert folds.deal(cohort) == folds.deal(shuffled)


@pytest.mark.parametrize(("fold", "k"), [("fold_10", folds.K), ("fold_5", folds.K_NESTED)])
def test_assign_balances_each_project(fold: str, k: int) -> None:
    rows = folds.assign(_cohort())
    assert _spread(rows, fold, k, lambda r: r["project_id"]) <= 1
    assert _spread(rows, fold, k, lambda r: (r["project_id"], r["pfi_status"])) <= 1
    events = [r for r in rows if r["os_status"] == "event"]
    assert _spread(events, fold, k, lambda r: r["project_id"]) <= 3


def test_fold_5_nests_in_fold_10() -> None:
    rows = folds.assign(_cohort())
    assert all(r["fold_5"] == r["fold_10"] % folds.K_NESTED for r in rows)
    assert len({r["case_id"] for r in rows}) == len(rows) == sum(_SIZES.values())


def test_load_cohort_prefers_liu_unless_absent_or_redacted(tmp_path: Path, monkeypatch) -> None:
    (tmp_path / "TCGA-AAA").mkdir()
    (tmp_path / "TCGA-AAA" / "cases.json").write_text("[]")
    patients = [
        {"case_id": f"c{i}", "case_submitter_id": f"TCGA-AA-000{i}", "project_id": "TCGA-AAA"}
        for i in range(3)
    ]
    liu = {
        "TCGA-AA-0000": {"cdr_PFI": 1, "cdr_OS": 0, "cdr_redaction": None},
        "TCGA-AA-0001": {"cdr_PFI": 1, "cdr_OS": 1, "cdr_redaction": "Redacted"},
    }

    def rederive(rows):
        for row in rows:
            row["survival_derived"] = {"pfi_event": 0, "os_event": None}

    monkeypatch.setattr(folds.cdr, "load_cdr_index", lambda raw: liu)
    monkeypatch.setattr(
        folds.clinical, "to_patient_rows", lambda cases: [dict(p) for p in patients]
    )
    monkeypatch.setattr(folds.clinical_supplement, "load_supplements_for_project", lambda d: {})
    monkeypatch.setattr(folds.survival, "attach_survival", rederive)

    got = {r["case_submitter_id"]: r for r in folds.load_cohort(tmp_path)}

    assert (got["TCGA-AA-0000"]["pfi_status"], got["TCGA-AA-0000"]["os_status"]) == (
        "event",
        "censored",
    )
    assert got["TCGA-AA-0000"]["status_source"] == folds.LIU
    for redacted_or_new in ("TCGA-AA-0001", "TCGA-AA-0002"):
        row = got[redacted_or_new]
        assert (row["pfi_status"], row["os_status"], row["status_source"]) == (
            "censored",
            "missing",
            folds.REDERIVED,
        )


def test_load_cohort_requires_the_cdr_workbook(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(folds.cdr, "load_cdr_index", lambda raw: {})
    with pytest.raises(FileNotFoundError, match="fetch-cdr"):
        folds.load_cohort(tmp_path)


def _built(tmp_path: Path) -> tuple[Path, str]:
    from ds2hf_pipeline.tcga.dataset_card import write_folds_card

    rows = folds.assign(_cohort())
    # Real-looking barcodes and projects so the key checks apply unchanged.
    for i, row in enumerate(rows):
        row["case_submitter_id"] = f"TCGA-{i // 1000:02d}-{i % 1000:04d}"
    folds.write(rows, tmp_path)
    card = write_folds_card(tmp_path, rows, ["Data Release 0"]).read_text()
    return tmp_path, card


def test_card_reproduce_snippet_runs(tmp_path: Path) -> None:
    """The card's own code must re-derive the published folds, and catch a moved one."""
    import pandas as pd

    processed, card = _built(tmp_path)
    snippet = card.split("## Reproduce the assignment")[1].split("```python")[1].split("```")[0]
    frame = pd.read_parquet(processed / "folds" / "data.parquet")

    exec(snippet, {"folds": frame})

    frame.loc[0, "fold_10"] = (frame.loc[0, "fold_10"] + 1) % folds.K
    with pytest.raises(AssertionError):
        exec(snippet, {"folds": frame})


def test_folds_card_names_no_other_hub_dataset(tmp_path: Path) -> None:
    from ds2hf_pipeline.tcga.dataset_card import FOLDS_REPO_ID

    _, card = _built(tmp_path)
    assert set(re.findall(r"gabrielaltay/[\w.-]+", card)) == {FOLDS_REPO_ID}
    assert "huggingface.co/datasets" not in card


def test_local_checks_pass_on_a_clean_build(tmp_path: Path) -> None:
    processed, _ = _built(tmp_path)
    for check in (
        verify.check_folds_tree,
        verify.check_folds_keys,
        verify.check_folds_nested,
        verify.check_folds_reproduce,
        verify.check_folds_balance,
    ):
        result = check(processed)
        assert result.passed, (result.name, result.summary, result.details)


def test_local_checks_catch_a_moved_patient_and_a_stray_file(tmp_path: Path) -> None:
    import pandas as pd

    processed, _ = _built(tmp_path)
    path = processed / "folds" / "data.parquet"
    frame = pd.read_parquet(path)
    frame.loc[0, "fold_10"] = (frame.loc[0, "fold_10"] + 1) % folds.K
    frame.to_parquet(path, index=False)
    (processed / ".DS_Store").write_text("")

    assert not verify.check_folds_tree(processed).passed
    assert not verify.check_folds_nested(processed).passed
    assert not verify.check_folds_reproduce(processed).passed
