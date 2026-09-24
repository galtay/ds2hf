"""Load GDC miRNA-Seq quantification files into consolidated patient records.

One file per aliquot, ~1,881 miRBase v21 mature miRNAs each. Unlike the seg
files these carry no aliquot column of their own, so the FK is the file's
single associated aliquot entity.

Two file spellings exist and both parse identically: `*.mirbase21.*` (TXT,
11,082 files) and `*.mirnaseq.*` (TSV, 359 files in GBM / OV / LUSC, from
the older Genome Analyzer platform). The fetch filter deliberately does not
pin `data_format` so the latter aren't silently dropped.

Records are struct-of-arrays (see `MIRNA_FIELDS`), matching `expression`.

Layout on disk:

    <data-dir>/raw/<project_id>/mirna/
        <uuid>.mirbase21.mirnas.quantification.txt
        manifest.json
"""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path
from typing import Any

import pandas as pd

MIRNA_DIR = "mirna"

# Source header spells this with a hyphen, which is neither a legal SQL
# identifier nor reachable as an attribute on a pandas row. Values untouched.
_CROSS_MAPPED_SOURCE = "cross-mapped"


def _case_id(entry: dict[str, Any]) -> str | None:
    case_ids = {c["case_id"] for c in (entry.get("cases") or []) if c.get("case_id")}
    return next(iter(case_ids)) if len(case_ids) == 1 else None


def _single_aliquot(entry: dict[str, Any]) -> str | None:
    """The file's one associated aliquot; None if absent or ambiguous."""
    aliquots = [
        e
        for e in (entry.get("associated_entities") or [])
        if e.get("entity_type") == "aliquot" and e.get("entity_id")
    ]
    return aliquots[0]["entity_id"] if len(aliquots) == 1 else None


def load_for_project(project_raw_dir: Path) -> dict[str, list[dict[str, Any]]]:
    """Return {case_id: [record, ...]} for raw/<PROJECT>/mirna/.

    Empty dict if the directory or its manifest is missing — letting `build`
    proceed without miRNA data when it hasn't been fetched. Manifest entries
    whose file isn't on disk are skipped, matching the other modalities.
    """
    mirna_dir = project_raw_dir / MIRNA_DIR
    manifest_path = mirna_dir / "manifest.json"
    if not manifest_path.exists():
        return {}

    by_case: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for entry in json.loads(manifest_path.read_text()):
        file_path = mirna_dir / entry["file_name"]
        if not file_path.exists():
            continue
        case_id = _case_id(entry)
        aliquot_id = _single_aliquot(entry)
        if not case_id or not aliquot_id:
            continue
        df = pd.read_csv(file_path, sep="\t", low_memory=False)
        df = df.rename(columns={_CROSS_MAPPED_SOURCE: "cross_mapped"})
        by_case[case_id].append(
            {
                "sample_id": None,  # resolved at attach time
                "aliquot_id": aliquot_id,
                "source_file_id": entry["file_id"],
                "mirna_id": df["miRNA_ID"].tolist(),
                "read_count": [None if pd.isna(v) else int(v) for v in df["read_count"]],
                "reads_per_million_mirna_mapped": [
                    None if pd.isna(v) else float(v) for v in df["reads_per_million_miRNA_mapped"]
                ],
                "cross_mapped": [None if pd.isna(v) else str(v) for v in df["cross_mapped"]],
            }
        )
    return dict(by_case)


def attach(
    rows: list[dict[str, Any]],
    by_case: dict[str, list[dict[str, Any]]],
) -> list[dict[str, Any]]:
    """Mutate `rows` to populate `samples_mirna_expression_quantification`.

    Resolves each record's `sample_id` from the patient row's `samples` tree.
    Sorted by aliquot_id for deterministic output; rows with no records get [].
    """
    from ds2hf_pipeline.tcga.copy_number import aliquot_to_sample

    for row in rows:
        records = by_case.get(row["case_id"], [])
        a2s = aliquot_to_sample(row)
        for r in records:
            r["sample_id"] = a2s.get(r["aliquot_id"])
        records.sort(key=lambda r: r.get("aliquot_id") or "")
        row["samples_mirna_expression_quantification"] = records
    return rows
