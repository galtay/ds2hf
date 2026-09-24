"""Standard cross-validation folds for TCGA patients.

One row per TCGA patient with a fold for 10-fold and for 5-fold cross
validation. Folds are assigned to patients, never samples, so every sample
of a patient lands in the same fold and no patient straddles train and test.

Assignment is deterministic. Within each project, patients are sorted by
PFI status, then OS status, then a salted hash of `case_id`, and dealt to
folds round-robin from a per-project offset. Dealing a sorted sequence
keeps every project, and every PFI status within it, within one patient
of even across folds. `fold_5` is `fold_10 % 5`, so each 5-fold fold is
exactly two 10-fold folds and inherits the same balance.

Survival status comes from Liu et al. 2018's TCGA-CDR where Liu covers the
patient, because most published TCGA survival work uses those labels.
Patients added to GDC after Liu's 2018 freeze, and those Liu marks
`Redacted`, take the OS/PFI that `survival.attach_survival` re-derives from
current GDC records instead. `status_source` records which applies. The
build reads the CDR workbook and, for the re-derivation, each project's
`cases.json` plus its Clinical Supplements.

The statuses ship in their own `stratification` config, apart from the
folds, so nobody loading folds picks them up as training labels.

Each build deals the current cohort afresh, so folds can change between
releases as GDC adds patients or updates follow-up. A Hub revision pins
one assignment; results should cite it.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

from ds2hf_pipeline.tcga import cdr, clinical, clinical_supplement, survival

K = 10
K_NESTED = 5
# Changing the salt reshuffles every patient: it names the v1 assignment.
SALT = "tcga-patient-folds:v1"
# Sort order within a project. Events first so the rarest status is dealt
# as one contiguous run.
STATUSES = ("event", "censored", "missing")
LIU = "liu_2018"
REDERIVED = "rederived"

# Two configs, one row per patient each, in the same order.
FOLDS_COLUMNS = ["case_id", "case_submitter_id", "project_id", "fold_10", "fold_5"]
STRATIFICATION_COLUMNS = ["case_id", "pfi_status", "os_status", "status_source"]
CONFIGS = ("folds", "stratification")
COLUMNS = FOLDS_COLUMNS + STRATIFICATION_COLUMNS[1:]


def _hash(text: str) -> str:
    return hashlib.sha256(f"{SALT}:{text}".encode()).hexdigest()


def project_offset(project_id: str) -> int:
    """Fold the project's first patient is dealt to, so fold 0 isn't always largest."""
    return int(_hash(project_id), 16) % K


def status(event: int | None) -> str:
    if event is None:
        return "missing"
    return "event" if event == 1 else "censored"


def load_cohort(raw_dir: Path) -> list[dict[str, Any]]:
    """Every TCGA patient under `raw_dir`, with the survival status to stratify on."""
    liu = cdr.load_cdr_index(raw_dir)
    if not liu:
        raise FileNotFoundError(
            f"no CDR workbook under {raw_dir / 'cdr'}; run `ds2hf-pipeline tcga fetch-cdr`"
        )
    cohort: list[dict[str, Any]] = []
    for cases_path in sorted(raw_dir.glob("TCGA-*/cases.json")):
        rows = clinical.to_patient_rows(json.loads(cases_path.read_text()))
        supps = clinical_supplement.load_supplements_for_project(
            cases_path.parent / "clinical_supplement"
        )
        if supps:
            clinical_supplement.attach_supplements(rows, supps)
        survival.attach_survival(rows)
        for row in rows:
            record = liu.get(row["case_submitter_id"])
            if record is not None and record["cdr_redaction"] != "Redacted":
                pfi, os_, source = record["cdr_PFI"], record["cdr_OS"], LIU
            else:
                derived = row["survival_derived"]
                pfi, os_, source = derived["pfi_event"], derived["os_event"], REDERIVED
            cohort.append(
                {
                    "case_id": row["case_id"],
                    "case_submitter_id": row["case_submitter_id"],
                    "project_id": row["project_id"],
                    "pfi_status": status(pfi),
                    "os_status": status(os_),
                    "status_source": source,
                }
            )
    return cohort


def _sort_key(row: dict[str, Any]) -> tuple[int, int, str]:
    return (
        STATUSES.index(row["pfi_status"]),
        STATUSES.index(row["os_status"]),
        _hash(row["case_id"]),
    )


def deal(cohort: list[dict[str, Any]]) -> dict[str, int]:
    """The stratified assignment: case_id -> fold_10. Input order is irrelevant."""
    by_project: dict[str, list[dict[str, Any]]] = {}
    for row in cohort:
        by_project.setdefault(row["project_id"], []).append(row)
    folds: dict[str, int] = {}
    for project_id, rows in by_project.items():
        offset = project_offset(project_id)
        for position, row in enumerate(sorted(rows, key=_sort_key)):
            folds[row["case_id"]] = (position + offset) % K
    return folds


def assign(cohort: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Deal every patient; return rows in COLUMNS order, sorted by project and barcode."""
    folds = deal(cohort)
    rows = [
        {**row, "fold_10": folds[row["case_id"]], "fold_5": folds[row["case_id"]] % K_NESTED}
        for row in cohort
    ]
    rows.sort(key=lambda r: (r["project_id"], r["case_submitter_id"]))
    return [{c: row[c] for c in COLUMNS} for row in rows]


def write(rows: list[dict[str, Any]], out_dir: Path) -> list[Path]:
    """Write the `folds` and `stratification` configs; return their paths."""
    import pyarrow as pa
    import pyarrow.parquet as pq

    types = {"fold_10": pa.int8(), "fold_5": pa.int8()}
    paths = []
    for config, columns in (("folds", FOLDS_COLUMNS), ("stratification", STRATIFICATION_COLUMNS)):
        schema = pa.schema([(c, types.get(c, pa.string())) for c in columns])
        table = pa.Table.from_pylist([{c: r[c] for c in columns} for r in rows], schema=schema)
        path = out_dir / config / "data.parquet"
        path.parent.mkdir(parents=True, exist_ok=True)
        pq.write_table(table, path)
        paths.append(path)
    return paths
