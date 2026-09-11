"""Tests for the GDC verification checks.

Only the checks that need no network are exercised here; the API-backed ones
are covered by running `verify-project` against a real project. What matters
is that a *failing* dataset actually fails — a check that always passes is
worse than no check.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
from tcga2hf_pipeline import verify

_ALIQUOT = "b3c68473-78a7-44d2-95bf-7997a7e6c0e8"


def _write_cases(processed: Path, aliquot_ids: list[str]) -> None:
    schema = pa.schema(
        [
            pa.field("case_id", pa.string()),
            pa.field(
                "samples",
                pa.list_(
                    pa.struct(
                        [
                            pa.field(
                                "portions",
                                pa.list_(
                                    pa.struct(
                                        [
                                            pa.field(
                                                "analytes",
                                                pa.list_(
                                                    pa.struct(
                                                        [
                                                            pa.field(
                                                                "aliquots",
                                                                pa.list_(
                                                                    pa.struct(
                                                                        [
                                                                            pa.field(
                                                                                "aliquot_id",
                                                                                pa.string(),
                                                                            )
                                                                        ]
                                                                    )
                                                                ),
                                                            )
                                                        ]
                                                    )
                                                ),
                                            )
                                        ]
                                    )
                                ),
                            )
                        ]
                    )
                ),
            ),
        ]
    )
    row = {
        "case_id": "c1",
        "samples": [
            {"portions": [{"analytes": [{"aliquots": [{"aliquot_id": a} for a in aliquot_ids]}]}]}
        ],
    }
    out = processed / "cases" / "data.parquet"
    out.parent.mkdir(parents=True)
    pq.write_table(pa.Table.from_pylist([row], schema=schema), out)


def _write_molecular(processed: Path, table: str, aliquot_ids: list[str]) -> None:
    out = processed / table / "data.parquet"
    out.parent.mkdir(parents=True)
    pq.write_table(pa.table({"aliquot_id": aliquot_ids}), out)


def test_fk_integrity_passes_when_every_aliquot_resolves(tmp_path: Path) -> None:
    _write_cases(tmp_path, [_ALIQUOT])
    _write_molecular(tmp_path, "gene_expression_quantification", [_ALIQUOT, _ALIQUOT])
    check = verify.check_fk_integrity(tmp_path)
    assert check.passed, check.details
    assert "0 orphan" in check.summary


def test_fk_integrity_catches_an_orphan_aliquot(tmp_path: Path) -> None:
    """The join every consumer relies on — a silent break here is the worst case."""
    _write_cases(tmp_path, [_ALIQUOT])
    _write_molecular(tmp_path, "gene_expression_quantification", [_ALIQUOT, "not-a-real-aliquot"])
    check = verify.check_fk_integrity(tmp_path)
    assert not check.passed
    assert any("not-a-real-aliquot" in d for d in check.details)


def test_fk_integrity_skips_tables_without_an_aliquot_column(tmp_path: Path) -> None:
    """`masked_somatic_mutation` keys on tumor_sample_id, and must not error."""
    _write_cases(tmp_path, [_ALIQUOT])
    out = tmp_path / "gene_expression_quantification" / "data.parquet"
    out.parent.mkdir(parents=True)
    pq.write_table(pa.table({"tumor_sample_id": ["s1"]}), out)
    check = verify.check_fk_integrity(tmp_path)
    assert check.passed
    assert any("no aliquot_id column" in d for d in check.details)


def _write_raw(raw: Path, modality: str, name: str, content: bytes, md5: str | None) -> None:
    d = raw / modality
    d.mkdir(parents=True, exist_ok=True)
    (d / name).write_bytes(content)
    digest = md5 if md5 is not None else hashlib.md5(content).hexdigest()  # noqa: S324
    (d / "manifest.json").write_text(
        json.dumps([{"file_name": name, "md5sum": digest, "file_id": "f1"}])
    )


def test_local_md5_passes_on_intact_bytes(tmp_path: Path) -> None:
    _write_raw(tmp_path, "expression", "a.tsv", b"gene\tvalue\n", None)
    check = verify.check_local_md5(tmp_path, sample=3)
    assert check.passed
    assert "1 raw files re-hashed, 0 mismatched" in check.summary


def test_local_md5_catches_a_corrupted_download(tmp_path: Path) -> None:
    """A truncated file keeps its name and manifest entry; only the hash moves."""
    _write_raw(tmp_path, "expression", "a.tsv", b"gene\tvalue\n", md5="0" * 32)
    check = verify.check_local_md5(tmp_path, sample=3)
    assert not check.passed
    assert any("md5 differs" in d for d in check.details)


def test_local_md5_samples_each_modality(tmp_path: Path) -> None:
    for modality in ("expression", "mutations", "methylation"):
        _write_raw(tmp_path, modality, "a.tsv", f"{modality}".encode(), None)
    check = verify.check_local_md5(tmp_path, sample=1)
    assert check.passed
    # One file per modality, so every modality is touched even at sample=1.
    assert "3 raw files re-hashed" in check.summary


def test_exclusions_are_declared_not_inferred() -> None:
    """A modality dropped from the pipeline must fail, not join the list."""
    assert set(verify.EXPECTED_EXCLUSIONS) == {"Slide Image", "Masked Intensities"}
    assert all(why for why in verify.EXPECTED_EXCLUSIONS.values())


# ---------------------------------------------------------------------------
# The `files` table as a complete index
# ---------------------------------------------------------------------------


def test_files_table_indexes_files_we_never_downloaded(tmp_path: Path) -> None:
    """A file in the GDC index but absent locally still gets a row and a URL."""
    from tcga2hf_pipeline import tabular

    raw = tmp_path / "TCGA-XYZ"
    (raw / "expression").mkdir(parents=True)
    (raw / "expression" / "a.tsv").write_text("x")
    (raw / "expression" / "manifest.json").write_text(
        json.dumps(
            [
                {
                    "file_id": "have-it",
                    "file_name": "a.tsv",
                    "_status": "downloaded",
                    "cases": [{"case_id": "c1", "submitter_id": "TCGA-XX-0001"}],
                }
            ]
        )
    )
    (raw / "files_index.json").write_text(
        json.dumps(
            [
                {
                    "file_id": "have-it",
                    "file_name": "a.tsv",
                    "data_type": "Gene Expression Quantification",
                    "cases": [{"case_id": "c1", "submitter_id": "TCGA-XX-0001"}],
                },
                {
                    "file_id": "never-fetched",
                    "file_name": "big.svs",
                    "data_type": "Slide Image",
                    "file_size": 900_000_000,
                    "cases": [{"case_id": "c1", "submitter_id": "TCGA-XX-0001"}],
                },
            ]
        )
    )
    rows = {r["file_id"]: r for r in tabular._files_rows(raw)}
    assert set(rows) == {"have-it", "never-fetched"}

    have = rows["have-it"]
    assert have["in_dataset"] is True
    assert have["dataset_table"] == "gene_expression_quantification"
    assert have["modality"] == "expression"

    absent = rows["never-fetched"]
    assert absent["in_dataset"] is False
    assert absent["dataset_table"] is None
    assert absent["modality"] is None
    # The whole point: you can still go get it.
    assert absent["gdc_download_url"] == "https://api.gdc.cancer.gov/data/never-fetched"
    assert absent["case_id"] == "c1"


def test_fetched_but_unpublished_files_report_in_dataset_false(tmp_path: Path) -> None:
    """BCR XML is downloaded yet redundant with the biotabs, so it isn't published."""
    from tcga2hf_pipeline import tabular

    raw = tmp_path / "TCGA-XYZ"
    (raw / "clinical_supplement_xml").mkdir(parents=True)
    (raw / "clinical_supplement_xml" / "c.xml").write_text("<x/>")
    (raw / "clinical_supplement_xml" / "manifest.json").write_text(
        json.dumps([{"file_id": "xml1", "file_name": "c.xml", "_status": "downloaded"}])
    )
    rows = tabular._files_rows(raw)
    assert len(rows) == 1
    assert rows[0]["modality"] == "clinical_supplement_xml"  # we do have the bytes
    assert rows[0]["in_dataset"] is False  # but they are not published
    assert rows[0]["dataset_table"] is None


def test_multi_case_files_report_a_null_case_id(tmp_path: Path) -> None:
    """A project-level biotab names every case; one row per case would lie."""
    from tcga2hf_pipeline import tabular

    raw = tmp_path / "TCGA-XYZ"
    raw.mkdir(parents=True)
    (raw / "files_index.json").write_text(
        json.dumps(
            [
                {
                    "file_id": "biotab",
                    "file_name": "clinical_patient_xyz.txt",
                    "data_type": "Clinical Supplement",
                    "cases": [{"case_id": "c1"}, {"case_id": "c2"}, {"case_id": "c3"}],
                }
            ]
        )
    )
    rows = tabular._files_rows(raw)
    assert len(rows) == 1
    assert rows[0]["case_id"] is None


def test_files_rows_still_work_without_an_index(tmp_path: Path) -> None:
    """Projects built before fetch-file-index existed must still build."""
    from tcga2hf_pipeline import tabular

    raw = tmp_path / "TCGA-XYZ"
    (raw / "mutations").mkdir(parents=True)
    (raw / "mutations" / "m.maf").write_text("x")
    (raw / "mutations" / "manifest.json").write_text(
        json.dumps(
            [{"file_id": "m1", "file_name": "m.maf", "_status": "downloaded", "cases": []}]
        )
    )
    rows = tabular._files_rows(raw)
    assert len(rows) == 1
    assert rows[0]["in_dataset"] is True
    assert rows[0]["dataset_table"] == "masked_somatic_mutation"
