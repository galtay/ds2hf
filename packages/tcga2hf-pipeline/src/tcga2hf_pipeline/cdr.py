"""Liu et al., Cell 2018 — TCGA Pan-Cancer Clinical Data Resource (CDR).

Curated survival endpoints (OS, DSS, DFI, PFI) for the TCGA cohort, frozen
at the 2018 data freeze. The file is hosted on the GDC's PanCanAtlas
auxiliary distribution under a stable UUID — fetchable via
`POST/GET /data/<uuid>` but **not indexed by `/files`**, so we don't go
through `genomic.fetch_files`. We download once into
`<data-dir>/raw/cdr/` and validate the md5 against the snapshot we tested
against.

The workbook has 8 sheets; only `TCGA-CDR` (the headline curated table)
is used. The sheet is keyed by `bcr_patient_barcode`, which equals our
`case_submitter_id`. Cases added to the GDC after Liu's 2018 freeze are
absent — `attach_cdr` flags those as `cdr_matched=False` rather than
silently dropping.

References:
  - Liu et al., Cell. 2018 Apr 5;173(2):400-416.e11
    https://doi.org/10.1016/j.cell.2018.02.052
  - Workbook: GDC PanCanAtlas, UUID 1b5f413e-a8d1-4d10-92eb-7c4ae739ed81
"""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any

import pandas as pd

from tcga2hf_pipeline.gdc import GDCClient

CDR_FILE_UUID = "1b5f413e-a8d1-4d10-92eb-7c4ae739ed81"
CDR_FILE_NAME = "TCGA-CDR-SupplementalTableS1.xlsx"
# md5 of the 2026-05 snapshot we tested against. Liu's CDR is frozen and
# this UUID points at a stable file, so a future hash mismatch would be
# noteworthy (likely the GDC re-hosted the same content with a different
# wrapper, or — much less likely — Liu issued a v2). Either way, surface
# the change rather than silently consume new bytes.
CDR_FILE_MD5 = "a4591b2dcee39591f59e5e25a6ce75fa"

# Sheet inside the workbook that holds the curated patient-level table.
_CDR_SHEET = "TCGA-CDR"

# Liu's join key on the curated sheet equals our `case_submitter_id`.
_CDR_JOIN_KEY = "bcr_patient_barcode"

# Survival columns we lift verbatim from the workbook into our schema.
# Column names in the workbook are e.g. `OS`, `OS.time`; we rename the
# `*.time` ones to `*_time` (parquet column names with `.` aren't a hard
# error but trip up several SQL engines).
_SURVIVAL_COLS: list[tuple[str, str]] = [
    ("OS", "cdr_OS"),
    ("OS.time", "cdr_OS_time"),
    ("DSS", "cdr_DSS"),
    ("DSS.time", "cdr_DSS_time"),
    ("DFI", "cdr_DFI"),
    ("DFI.time", "cdr_DFI_time"),
    ("PFI", "cdr_PFI"),
    ("PFI.time", "cdr_PFI_time"),
]
_REDACTION_COL = "Redaction"


def fetch_cdr_workbook(raw_dir: Path) -> Path:
    """Download the CDR workbook into `<raw_dir>/cdr/` if not already cached.

    The destination filename matches GDC's `Content-Disposition`
    (`TCGA-CDR-SupplementalTableS1.xlsx`). Returns the path. If the file is
    already on disk, only re-downloads when its md5 doesn't match
    `CDR_FILE_MD5`.
    """
    cdr_dir = raw_dir / "cdr"
    cdr_dir.mkdir(parents=True, exist_ok=True)
    out_path = cdr_dir / CDR_FILE_NAME

    if out_path.exists():
        existing_md5 = hashlib.md5(out_path.read_bytes()).hexdigest()
        if existing_md5 == CDR_FILE_MD5:
            return out_path
        # md5 drift — re-download.

    with GDCClient() as client:
        client.download(CDR_FILE_UUID, out_path)

    actual_md5 = hashlib.md5(out_path.read_bytes()).hexdigest()
    if actual_md5 != CDR_FILE_MD5:
        # Surface the drift but don't fail — downstream code will work with
        # whatever the workbook contains; we just want the operator to know.
        print(
            f"warning: CDR workbook md5 changed: expected {CDR_FILE_MD5}, "
            f"got {actual_md5}. The pinned snapshot may have been updated."
        )
    return out_path


def load_cdr_index(raw_dir: Path) -> dict[str, dict[str, Any]]:
    """Read the curated `TCGA-CDR` sheet, return one dict per patient.

    Keyed by `bcr_patient_barcode` (== `case_submitter_id`). Each value
    contains the renamed survival columns plus the Redaction flag, with
    NaN normalized to `None` so pyarrow nullable typing round-trips
    cleanly.
    """
    cdr_path = raw_dir / "cdr" / CDR_FILE_NAME
    if not cdr_path.exists():
        return {}

    df = pd.read_excel(cdr_path, sheet_name=_CDR_SHEET)
    if _CDR_JOIN_KEY not in df.columns:
        raise RuntimeError(f"CDR workbook sheet {_CDR_SHEET!r} missing join key {_CDR_JOIN_KEY!r}")

    out: dict[str, dict[str, Any]] = {}
    for _, row in df.iterrows():
        patient = row[_CDR_JOIN_KEY]
        if pd.isna(patient):
            continue
        record: dict[str, Any] = {}
        # Survival event columns are 0/1 stored as float because of NaNs;
        # cast to int when populated, leave None otherwise.
        for src, dst in _SURVIVAL_COLS:
            v = row.get(src)
            if pd.isna(v):
                record[dst] = None
            elif dst.endswith("_time"):
                record[dst] = float(v)
            else:
                record[dst] = int(v)
        red = row.get(_REDACTION_COL)
        record["cdr_redaction"] = None if pd.isna(red) else str(red)
        out[str(patient)] = record
    return out


def attach_cdr(
    rows: list[dict[str, Any]], index: dict[str, dict[str, Any]]
) -> list[dict[str, Any]]:
    """Mutate patient rows to populate the `cdr_*` columns + audit flags.

    For matched cases (those whose `case_submitter_id` is in `index`):
    - `cdr_matched = True`
    - the 8 survival columns and `cdr_redaction` get the workbook values
    - `cdr_survival_complete` is True iff every survival column is non-null

    For unmatched cases: `cdr_matched = False`, every other `cdr_*`
    column stays None.
    """
    cdr_cols = [dst for _, dst in _SURVIVAL_COLS]
    for row in rows:
        record = index.get(row.get("case_submitter_id"))
        if record is None:
            row["cdr_matched"] = False
            row["cdr_redaction"] = None
            for c in cdr_cols:
                row[c] = None
            row["cdr_survival_complete"] = False
            continue
        row["cdr_matched"] = True
        row["cdr_redaction"] = record["cdr_redaction"]
        for c in cdr_cols:
            row[c] = record[c]
        row["cdr_survival_complete"] = all(record[c] is not None for c in cdr_cols)
    return rows
