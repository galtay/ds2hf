"""Load GDC copy number segment files into consolidated patient records.

Two data types, one module, because they share everything but their value
columns:

  - **Allele-specific** (ASCAT2 / ASCAT3 / AscatNGS) — integer total copy
    number with its major/minor allelic split. ASCAT is a paired caller, so
    each file names two aliquots; the tumour is the one the file's own
    `GDC_Aliquot` column reports, and the other is the matched normal.
  - **Masked** (DNAcopy) — log2 ratio against a diploid reference with
    germline CNVs masked out. Single aliquot.

All three allele-specific callers ship for overlapping aliquots and fit
purity and ploidy independently, so one case can carry several records for
the same tumour aliquot that disagree with each other. `workflow_type` is
therefore scalar on the record — one record per (aliquot, workflow) — which
makes selecting a caller a filter on the struct rather than inside the
arrays.

Records are struct-of-arrays (see `ALLELE_SPECIFIC_CNV_FIELDS` /
`MASKED_CNV_FIELDS`), matching how `expression` shapes its per-aliquot
payload.

Layout on disk:

    <data-dir>/raw/<project_id>/copy_number_allele_specific/
    <data-dir>/raw/<project_id>/copy_number_masked/
        <file>.seg.txt
        manifest.json
"""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path
from typing import Any

import pandas as pd

ALLELE_SPECIFIC_DIR = "copy_number_allele_specific"
MASKED_DIR = "copy_number_masked"


def _case_id(entry: dict[str, Any]) -> str | None:
    """The single case this file belongs to; None if absent or ambiguous."""
    case_ids = {c["case_id"] for c in (entry.get("cases") or []) if c.get("case_id")}
    return next(iter(case_ids)) if len(case_ids) == 1 else None


def _aliquot_entities(entry: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        e
        for e in (entry.get("associated_entities") or [])
        if e.get("entity_type") == "aliquot" and e.get("entity_id")
    ]


def _seg_aliquot(df: pd.DataFrame) -> str | None:
    """The single aliquot a seg file's own `GDC_Aliquot` column reports."""
    values = set(df["GDC_Aliquot"].dropna().unique())
    return next(iter(values)) if len(values) == 1 else None


def _int_list(series: pd.Series) -> list[int | None]:
    return [None if pd.isna(v) else int(v) for v in series]


def _float_list(series: pd.Series) -> list[float | None]:
    return [None if pd.isna(v) else float(v) for v in series]


def _str_list(series: pd.Series) -> list[str | None]:
    return [None if pd.isna(v) else str(v) for v in series]


def _iter_files(project_raw_dir: Path, modality_dir: str) -> Any:
    """Yield `(entry, file_path, case_id)` for manifest entries present on disk."""
    mod_dir = project_raw_dir / modality_dir
    manifest_path = mod_dir / "manifest.json"
    if not manifest_path.exists():
        return
    for entry in json.loads(manifest_path.read_text()):
        file_path = mod_dir / entry["file_name"]
        if not file_path.exists():
            continue
        case_id = _case_id(entry)
        if not case_id:
            continue
        yield entry, file_path, case_id


def load_allele_specific_for_project(
    project_raw_dir: Path,
) -> dict[str, list[dict[str, Any]]]:
    """Return {case_id: [record, ...]} for raw/<PROJECT>/copy_number_allele_specific/.

    Empty dict if the directory or its manifest is missing. Files whose
    `GDC_Aliquot` column isn't a single consistent value are skipped rather
    than guessed at — without it there is no defensible tumour FK.
    """
    by_case: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for entry, file_path, case_id in _iter_files(project_raw_dir, ALLELE_SPECIFIC_DIR):
        df = pd.read_csv(file_path, sep="\t", low_memory=False)
        if df.empty:
            continue
        tumor_aliquot = _seg_aliquot(df)
        if tumor_aliquot is None:
            continue
        others = [
            e["entity_id"] for e in _aliquot_entities(entry) if e["entity_id"] != tumor_aliquot
        ]
        by_case[case_id].append(
            {
                "sample_id": None,  # resolved at attach time
                "aliquot_id": tumor_aliquot,
                "matched_normal_aliquot_id": others[0] if len(others) == 1 else None,
                "workflow_type": entry.get("workflow_type"),
                "experimental_strategy": entry.get("experimental_strategy"),
                "source_file_id": entry["file_id"],
                "chromosome": _str_list(df["Chromosome"]),
                "start": _int_list(df["Start"]),
                "end": _int_list(df["End"]),
                "copy_number": _int_list(df["Copy_Number"]),
                "major_copy_number": _int_list(df["Major_Copy_Number"]),
                "minor_copy_number": _int_list(df["Minor_Copy_Number"]),
            }
        )
    return dict(by_case)


def load_masked_for_project(project_raw_dir: Path) -> dict[str, list[dict[str, Any]]]:
    """Return {case_id: [record, ...]} for raw/<PROJECT>/copy_number_masked/."""
    by_case: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for entry, file_path, case_id in _iter_files(project_raw_dir, MASKED_DIR):
        df = pd.read_csv(file_path, sep="\t", low_memory=False)
        if df.empty:
            continue
        aliquot_id = _seg_aliquot(df)
        if aliquot_id is None:
            continue
        by_case[case_id].append(
            {
                "sample_id": None,
                "aliquot_id": aliquot_id,
                "workflow_type": entry.get("workflow_type"),
                "source_file_id": entry["file_id"],
                # Bare names in this data type; str() guards against pandas
                # reading a numeric-looking column as int64.
                "chromosome": _str_list(df["Chromosome"]),
                "start": _int_list(df["Start"]),
                "end": _int_list(df["End"]),
                "num_probes": _int_list(df["Num_Probes"]),
                "segment_mean": _float_list(df["Segment_Mean"]),
            }
        )
    return dict(by_case)


def aliquot_to_sample(row: dict[str, Any]) -> dict[str, str]:
    """Map aliquot_id -> sample_id from a built patient row's `samples` tree."""
    out: dict[str, str] = {}
    for s in row.get("samples") or []:
        sid = s.get("sample_id")
        if not sid:
            continue
        for portion in s.get("portions") or []:
            for analyte in portion.get("analytes") or []:
                for a in analyte.get("aliquots") or []:
                    aq = a.get("aliquot_id")
                    if aq:
                        out[aq] = sid
    return out


def attach(
    rows: list[dict[str, Any]],
    by_case: dict[str, list[dict[str, Any]]],
    column: str,
) -> list[dict[str, Any]]:
    """Mutate `rows` to populate `column` with the copy number records.

    Resolves each record's `sample_id` from the patient row's own `samples`
    tree (the source of truth for every aliquot the case has). Sorted by
    (aliquot_id, workflow_type) so a case carrying several callers for one
    aliquot has a deterministic order; rows with no records get [].
    """
    for row in rows:
        records = by_case.get(row["case_id"], [])
        a2s = aliquot_to_sample(row)
        for r in records:
            r["sample_id"] = a2s.get(r["aliquot_id"])
        records.sort(key=lambda r: (r.get("aliquot_id") or "", r.get("workflow_type") or ""))
        row[column] = records
    return rows
