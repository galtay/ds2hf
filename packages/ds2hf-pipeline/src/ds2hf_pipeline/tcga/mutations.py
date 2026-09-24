"""Parse open-access GDC Masked Somatic Mutation MAFs into per-case variant lists.

Each MAF file is per-(tumor_aliquot, matched_normal_aliquot) pair. We keep all
140 MAF columns (preserving GDC source fidelity) and add three FK fields that
join back to the patient row's `samples[].aliquots[]` entries:
  - tumor_sample_id (parent sample of MAF.Tumor_Sample_UUID)
  - matched_normal_sample_id (parent sample of MAF.Matched_Norm_Sample_UUID)
  - source_file_id (the MAF's own GDC file UUID)

The MAF already carries `case_id` (col 132) and the aliquot UUIDs themselves.

GDC MAF format spec:
  https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/
"""

from __future__ import annotations

import gzip
import json
from collections import defaultdict
from pathlib import Path
from typing import Any

import pandas as pd
from ds2hf.tcga.schema import _MAF_COLUMNS, _MAF_FLOAT_COLS, _MAF_INT_COLS, MUTATION_FIELDS

_FK_NAMES = {f.name for f in MUTATION_FIELDS} - set(_MAF_COLUMNS)


def _aliquot_to_sample_from_patient(samples: list[dict[str, Any]]) -> dict[str, str]:
    """Build {aliquot_id: sample_id} by walking the patient's full GDC tree
    samples → portions → analytes → aliquots.

    The patient row carries every aliquot the GDC has for the case (from the
    /cases endpoint), unlike the per-MAF manifest which only carries samples
    directly associated with the MAF file. Use the patient row as source of
    truth for resolving Tumor_Sample_UUID -> tumor_sample_id.
    """
    out: dict[str, str] = {}
    for sample in samples:
        sample_id = sample.get("sample_id")
        if not sample_id:
            continue
        for portion in sample.get("portions") or []:
            for analyte in portion.get("analytes") or []:
                for aliquot in analyte.get("aliquots") or []:
                    aq_id = aliquot.get("aliquot_id")
                    if aq_id:
                        out[aq_id] = sample_id
    return out


def _read_maf(path: Path) -> pd.DataFrame:
    """Read a gzipped MAF, skipping its variable-length leading `#` header.

    Forces every MAF column to string at read time, then casts numeric columns
    to their schema-declared type. This avoids pandas' inferred-dtype
    inconsistency across files (esp. all-null columns becoming float64).
    """
    skip = 0
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                skip += 1
            else:
                break
    df = pd.read_csv(
        path,
        sep="\t",
        skiprows=skip,
        compression="gzip",
        dtype=str,
        keep_default_na=False,
        na_values=[""],
        low_memory=False,
    )
    # Cast numeric columns. Errors become NaN so all-null / blank stays null.
    for col in _MAF_INT_COLS & set(df.columns):
        df[col] = pd.to_numeric(df[col], errors="coerce").astype("Int64")
    for col in _MAF_FLOAT_COLS & set(df.columns):
        df[col] = pd.to_numeric(df[col], errors="coerce").astype("Float64")
    return df


def _row_to_variant(row: pd.Series, *, source_file_id: str) -> dict[str, Any]:
    """Build a variant record from one MAF row. FK sample_ids are resolved later
    in `attach` against the patient's full samples list."""
    record: dict[str, Any] = {}
    for col in _MAF_COLUMNS:
        v = row.get(col)
        # pandas NA / NaN -> None for pyarrow
        if v is pd.NA or (isinstance(v, float) and pd.isna(v)):
            record[col] = None
        else:
            record[col] = v
    record["source_file_id"] = source_file_id
    record["tumor_sample_id"] = None
    record["matched_normal_sample_id"] = None
    return record


def load_for_project(project_raw_dir: Path) -> dict[str, list[dict[str, Any]]]:
    """Return {case_id: [variant_dict, ...]} for all MAFs in raw/<PROJECT>/mutations/.

    Empty dict if the mutations directory or its manifest is missing — letting
    `build` proceed clinical-only when mutations haven't been fetched.
    """
    mutations_dir = project_raw_dir / "mutations"
    manifest_path = mutations_dir / "manifest.json"
    if not manifest_path.exists():
        return {}

    manifest = json.loads(manifest_path.read_text())
    by_case: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for entry in manifest:
        file_path = mutations_dir / entry["file_name"]
        if not file_path.exists():
            continue
        df = _read_maf(file_path)
        for _, row in df.iterrows():
            variant = _row_to_variant(row, source_file_id=entry["file_id"])
            case_id = variant.get("case_id")
            if case_id:
                by_case[case_id].append(variant)
    return dict(by_case)


def attach(
    rows: list[dict[str, Any]], by_case: dict[str, list[dict[str, Any]]]
) -> list[dict[str, Any]]:
    """Mutate `rows` to populate `samples_masked_somatic_mutation` from `by_case`.

    Resolves tumor_sample_id / matched_normal_sample_id by looking up the MAF's
    aliquot UUIDs in the patient row's `samples` field (the source of truth for
    every aliquot the case has). Sorts variants per row by
    (Chromosome, Start_Position, Tumor_Sample_UUID) for deterministic ordering.
    Rows with no entries get [].
    """
    for row in rows:
        variants = by_case.get(row["case_id"], [])
        aliquot_to_sample = _aliquot_to_sample_from_patient(row.get("samples") or [])
        for v in variants:
            v["tumor_sample_id"] = aliquot_to_sample.get(v.get("Tumor_Sample_UUID"))
            v["matched_normal_sample_id"] = aliquot_to_sample.get(v.get("Matched_Norm_Sample_UUID"))
        variants.sort(
            key=lambda v: (
                str(v.get("Chromosome") or ""),
                v.get("Start_Position") or 0,
                v.get("Tumor_Sample_UUID") or "",
            )
        )
        row["samples_masked_somatic_mutation"] = variants
    return rows
