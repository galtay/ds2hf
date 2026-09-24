"""Load GDC RPPA protein expression files into consolidated patient records.

Reverse Phase Protein Array, ~487 antibodies per file. This is the only
modality in the dataset that attaches to a **portion** rather than an
aliquot, so the sample FK is resolved by walking samples -> portions rather
than samples -> portions -> analytes -> aliquots.

Coverage is the narrowest we ship (7,827 of 11,428 TCGA cases) because RPPA
was only run on a subset, and the antibody panel grew over the project's
life — `set_id` records which version a measurement came from, so a target
absent for a portion may mean "not on that panel" rather than "measured as
zero".

`protein_expression` is null where the source says the literal string `NA`
— a failed measurement, not a zero. Around 5.5% of cells pan-cancer.

Records are struct-of-arrays (see `PROTEIN_EXPRESSION_FIELDS`).

Layout on disk:

    <data-dir>/raw/<project_id>/protein_expression/
        TCGA-W5-AA2Q-01A-21-A45N-20_RPPA_data.tsv
        manifest.json
"""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path
from typing import Any

import pandas as pd

PROTEIN_EXPRESSION_DIR = "protein_expression"


def _case_id(entry: dict[str, Any]) -> str | None:
    case_ids = {c["case_id"] for c in (entry.get("cases") or []) if c.get("case_id")}
    return next(iter(case_ids)) if len(case_ids) == 1 else None


def _single_portion(entry: dict[str, Any]) -> str | None:
    """The file's one associated portion; None if absent or ambiguous."""
    portions = [
        e
        for e in (entry.get("associated_entities") or [])
        if e.get("entity_type") == "portion" and e.get("entity_id")
    ]
    return portions[0]["entity_id"] if len(portions) == 1 else None


def _str_list(series: pd.Series) -> list[str | None]:
    """Identifier columns kept as strings so they never gain a decimal point."""
    return [None if pd.isna(v) else str(v) for v in series]


def load_for_project(project_raw_dir: Path) -> dict[str, list[dict[str, Any]]]:
    """Return {case_id: [record, ...]} for raw/<PROJECT>/protein_expression/.

    Empty dict if the directory or its manifest is missing.
    """
    rppa_dir = project_raw_dir / PROTEIN_EXPRESSION_DIR
    manifest_path = rppa_dir / "manifest.json"
    if not manifest_path.exists():
        return {}

    by_case: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for entry in json.loads(manifest_path.read_text()):
        file_path = rppa_dir / entry["file_name"]
        if not file_path.exists():
            continue
        case_id = _case_id(entry)
        portion_id = _single_portion(entry)
        if not case_id or not portion_id:
            continue
        # `NA` in protein_expression is pandas' default NaN token, so the
        # failed measurements arrive as NaN and become null below.
        df = pd.read_csv(file_path, sep="\t", low_memory=False)
        by_case[case_id].append(
            {
                "sample_id": None,  # resolved at attach time
                "portion_id": portion_id,
                "source_file_id": entry["file_id"],
                "agid": _str_list(df["AGID"]),
                "lab_id": _str_list(df["lab_id"]),
                "catalog_number": _str_list(df["catalog_number"]),
                "set_id": _str_list(df["set_id"]),
                "peptide_target": _str_list(df["peptide_target"]),
                "protein_expression": [
                    None if pd.isna(v) else float(v) for v in df["protein_expression"]
                ],
            }
        )
    return dict(by_case)


def portion_to_sample(row: dict[str, Any]) -> dict[str, str]:
    """Map portion_id -> sample_id from a built patient row's `samples` tree."""
    out: dict[str, str] = {}
    for s in row.get("samples") or []:
        sid = s.get("sample_id")
        if not sid:
            continue
        for p in s.get("portions") or []:
            pid = p.get("portion_id")
            if pid:
                out[pid] = sid
    return out


def attach(
    rows: list[dict[str, Any]],
    by_case: dict[str, list[dict[str, Any]]],
) -> list[dict[str, Any]]:
    """Mutate `rows` to populate `samples_protein_expression_quantification`.

    Resolves each record's `sample_id` through its portion. A portion GDC
    names but the case tree doesn't report leaves `sample_id` null rather
    than dropping the measurement. Sorted by portion_id for deterministic
    output; rows with no records get [].
    """
    for row in rows:
        records = by_case.get(row["case_id"], [])
        p2s = portion_to_sample(row)
        for r in records:
            r["sample_id"] = p2s.get(r["portion_id"])
        records.sort(key=lambda r: r.get("portion_id") or "")
        row["samples_protein_expression_quantification"] = records
    return rows
