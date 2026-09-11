"""Row emitters for the tabular HF dataset (`tcga2hf-tabular-open`).

Companion to `tcga2hf.clinical` (which produces the consolidated patient
row with both clinical and molecular vectors nested). The tabular layout
splits each project into four tables, named after stable GDC concepts:

  - `cases`                            GDC `cases.json` row, nested
                                       structure preserved (one row per
                                       patient)
  - `masked_somatic_mutation`          one row per MAF variant — snake-
                                       cased GDC data_type for somatic
                                       calls (today: ensemble-masked
                                       WXS-derived MAFs)
  - `gene_expression_quantification`   one row per (aliquot, gene) —
                                       snake-cased GDC data_type for
                                       RNA-Seq expression (today: STAR
                                       counts TSVs)
  - `pathology_report`                 one row per scanned report PDF —
                                       the document bytes verbatim, no
                                       text extraction (see
                                       `tcga2hf_pipeline.pathology`)
  - `ssgsea_scores_<collection>`       one row per (aliquot, pathway) of
                                       raw ssGSEA pathway activity, one
                                       table per MSigDB collection
  - `ssgsea_stats_<collection>`        reference distributions for those
                                       scores (this project + pan-cancer),
                                       written after every project scores
  - `files`                            one row per (file, case) from the
                                       per-modality manifests

Two further tables are built elsewhere and land in the same tree:
`survival_derived` (projected off the cases rows by `derived_survival_rows`
after the caller attaches survival) and the flex-schema
`clinical_supplement_*` set (one per BCR biotab form).

The `cases` table reuses `clinical.to_patient_rows` and drops the
per-file modality columns. The molecular, document, and provenance tables
read the same raw inputs the consolidated build does.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from tcga2hf import schema as _schema_mod
from tcga2hf.schema import SSGSEA_COLLECTIONS, TABULAR_CASES_FIELDS, TABULAR_TABLES

from tcga2hf_pipeline import biospecimen_supplement as _biospecimen_supplement_mod
from tcga2hf_pipeline import clinical as _clinical_mod
from tcga2hf_pipeline import clinical_supplement as _clinical_supplement_mod
from tcga2hf_pipeline import expression as _expression_mod
from tcga2hf_pipeline import msigdb as _msigdb_mod
from tcga2hf_pipeline import mutations as _mutations_mod
from tcga2hf_pipeline import pathology as _pathology_mod
from tcga2hf_pipeline import ssgsea as _ssgsea_mod

# Names of `cases` table columns, used to re-shape patient rows. Excludes
# the two molecular-vector columns that the consolidated row carries (they
# have their own flat tables in the tabular layout).
_CASES_COLS: list[str] = [f.name for f in TABULAR_CASES_FIELDS]


# ---------------------------------------------------------------------------
# Cases (one nested row per patient — mirrors `cases.json` shape)
# ---------------------------------------------------------------------------


def _cases_rows(cases: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Build one nested row per case for the `cases` tabular table.

    Reuses `clinical.to_patient_rows` (which already mirrors the GDC
    `cases.json` structure with demographic / diagnoses / follow_ups /
    exposures / family_histories / samples nested), then drops the two
    molecular-vector columns those rows carry by default.
    """
    rows = _clinical_mod.to_patient_rows(cases)
    return [{k: row[k] for k in _CASES_COLS} for row in rows]


# ---------------------------------------------------------------------------
# Molecular emitters (mutations + expression)
# ---------------------------------------------------------------------------


def _aliquot_to_sample(case: dict[str, Any]) -> dict[str, str]:
    """Map aliquot_id → sample_id for one case (walks samples → … → aliquots)."""
    out: dict[str, str] = {}
    for s in case.get("samples") or []:
        sid = s.get("sample_id")
        if not sid:
            continue
        for p in s.get("portions") or []:
            for an in p.get("analytes") or []:
                for aq in an.get("aliquots") or []:
                    aid = aq.get("aliquot_id")
                    if aid:
                        out[aid] = sid
    return out


def _aliquot_to_submitter(case: dict[str, Any]) -> dict[str, str]:
    """Map aliquot_id → aliquot.submitter_id (for browseability on the
    expression table, where we don't otherwise carry the submitter id)."""
    out: dict[str, str] = {}
    for s in case.get("samples") or []:
        for p in s.get("portions") or []:
            for an in p.get("analytes") or []:
                for aq in an.get("aliquots") or []:
                    aid = aq.get("aliquot_id")
                    if aid:
                        out[aid] = aq.get("submitter_id")
    return out


def _mutations_rows(
    cases: list[dict[str, Any]],
    by_case: dict[str, list[dict[str, Any]]],
) -> list[dict[str, Any]]:
    """One row per MAF variant, with FKs resolved against the case's samples."""
    rows: list[dict[str, Any]] = []
    case_index = {c.get("case_id"): c for c in cases}
    for case_id, variants in by_case.items():
        case = case_index.get(case_id)
        if case is None:
            continue
        case_submitter_id = case.get("submitter_id")
        a2s = _aliquot_to_sample(case)
        for v in variants:
            row = dict(v)
            row["case_submitter_id"] = case_submitter_id
            row["tumor_sample_id"] = a2s.get(v.get("Tumor_Sample_UUID"))
            row["matched_normal_sample_id"] = a2s.get(v.get("Matched_Norm_Sample_UUID"))
            rows.append(row)
    return rows


def _expression_batches(
    cases: list[dict[str, Any]],
    project_raw_dir: Path,
) -> Any:
    """Yield one `pa.RecordBatch` per aliquot's STAR expression TSV.

    Emitted as a stream for the same reason `_gene_level_copy_number_batches`
    is: TCGA-BRCA is 1,231 aliquots x 60,660 genes = 74.7M rows, which as
    Python dicts is tens of GB of interpreter objects held live while every
    other table is still being built. One batch per file caps peak memory at
    ~60k rows.

    The 4 N_* alignment-summary rows at the top of each TSV (`N_unmapped`,
    `N_multimapping`, `N_noFeature`, `N_ambiguous` — empty gene_name /
    gene_type, and not genes) are dropped: they share the column shape but
    not the semantics, and they complicate every downstream filter.

    `gene_name` / `gene_type` are not emitted; they are invariant across
    files and live in `gene_model`.
    """
    import pyarrow as pa

    schema = TABULAR_TABLES["gene_expression_quantification"]
    expr_dir = project_raw_dir / "expression"
    manifest_path = expr_dir / "manifest.json"
    if not manifest_path.exists():
        return

    case_index = {c.get("case_id"): c for c in cases}
    value_columns = [
        f.name
        for f in schema
        if f.name
        not in {
            "case_id",
            "case_submitter_id",
            "aliquot_id",
            "aliquot_submitter_id",
            "source_file_id",
            "gene_id",
        }
    ]
    for entry in json.loads(manifest_path.read_text()):
        file_path = expr_dir / entry["file_name"]
        if not file_path.exists():
            continue
        case_id, aliquot_id = _expression_mod._file_aliquot_and_case(entry)
        if not case_id or not aliquot_id:
            continue
        case = case_index.get(case_id)
        if case is None:
            continue
        df = pd.read_csv(
            file_path,
            sep="\t",
            comment="#",
            usecols=["gene_id", *value_columns],
            low_memory=False,
        )
        df = df[~df["gene_id"].str.startswith("N_", na=False)].reset_index(drop=True)
        if df.empty:
            continue
        n = len(df)
        common = {
            "case_id": case_id,
            "case_submitter_id": case.get("submitter_id"),
            "aliquot_id": aliquot_id,
            "aliquot_submitter_id": _aliquot_to_submitter(case).get(aliquot_id),
            "source_file_id": entry["file_id"],
        }
        arrays = []
        for field in schema:
            if field.name in common:
                arrays.append(pa.array([common[field.name]] * n, type=field.type))
            elif field.name == "gene_id":
                arrays.append(pa.array(df["gene_id"].astype("string"), type=field.type))
            else:
                arrays.append(pa.array(df[field.name], type=field.type, from_pandas=True))
        yield pa.RecordBatch.from_arrays(arrays, schema=schema)


def _na_to_none(v: Any) -> Any:
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return None
    return v


def _to_int_or_none(v: Any) -> int | None:
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return None
    return int(v)


def _to_float_or_none(v: Any) -> float | None:
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return None
    return float(v)


def _pathology_report_rows(
    cases: list[dict[str, Any]],
    project_raw_dir: Path,
) -> list[dict[str, Any]]:
    """One row per pathology report PDF, with case FKs prepended.

    Reuses the consolidated loader (which reads the bytes and resolves the
    sample FK from `associated_entities`), then applies the same filename-UUID
    fallback `pathology.attach` does — here against the raw GDC case dict
    rather than a built patient row.
    """
    by_case = _pathology_mod.load_for_project(project_raw_dir)
    if not by_case:
        return []

    rows: list[dict[str, Any]] = []
    case_index = {c.get("case_id"): c for c in cases}
    for case_id, reports in by_case.items():
        case = case_index.get(case_id)
        if case is None:
            continue
        uuid_to_sample = {
            s["pathology_report_uuid"]: s
            for s in (case.get("samples") or [])
            if s.get("pathology_report_uuid")
        }
        for r in reports:
            row = dict(r)
            if row["sample_id"] is None:
                fallback = uuid_to_sample.get(row["pathology_report_uuid"] or "")
                if fallback:
                    row["sample_id"] = fallback.get("sample_id")
                    row["sample_submitter_id"] = fallback.get("submitter_id")
            row["case_id"] = case_id
            row["case_submitter_id"] = case.get("submitter_id")
            rows.append(row)
    rows.sort(key=lambda r: (r.get("case_submitter_id") or "", r.get("file_name") or ""))
    return rows


# ---------------------------------------------------------------------------
# Copy number, miRNA, protein expression
#
# These four emitters share a shape: walk a modality's manifest, resolve the
# biospecimen FKs off the raw GDC case dict, parse one flat TSV per file.
# They differ only in which entity the file attaches to and what its columns
# mean, so the walking and FK resolution are factored out here.
# ---------------------------------------------------------------------------


def _iter_modality_files(
    project_raw_dir: Path,
    modality_dir: str,
    cases: list[dict[str, Any]],
) -> Any:
    """Yield `(entry, file_path, case)` for manifest entries present on disk.

    Entries whose bytes weren't downloaded (`_status="manifest_only"`, from a
    capped fetch) and entries that don't resolve to exactly one known case are
    skipped — the same contract every other modality loader honours.
    """
    mod_dir = project_raw_dir / modality_dir
    manifest_path = mod_dir / "manifest.json"
    if not manifest_path.exists():
        return
    case_index = {c.get("case_id"): c for c in cases}
    for entry in json.loads(manifest_path.read_text()):
        file_path = mod_dir / entry["file_name"]
        if not file_path.exists():
            continue
        case_ids = {c["case_id"] for c in (entry.get("cases") or []) if c.get("case_id")}
        if len(case_ids) != 1:
            continue
        case = case_index.get(next(iter(case_ids)))
        if case is None:
            continue
        yield entry, file_path, case


def _sample_for_aliquot(case: dict[str, Any], aliquot_id: str | None) -> dict[str, Any]:
    """The sample dict owning `aliquot_id`, or {} when it can't be resolved."""
    if not aliquot_id:
        return {}
    sample_id = _aliquot_to_sample(case).get(aliquot_id)
    if not sample_id:
        return {}
    return next(
        (s for s in (case.get("samples") or []) if s.get("sample_id") == sample_id),
        {},
    )


def _aliquot_fks(case: dict[str, Any], aliquot_id: str | None) -> dict[str, Any]:
    """Case + sample + aliquot FK columns for one aliquot-scoped measurement."""
    sample = _sample_for_aliquot(case, aliquot_id)
    return {
        "case_id": case.get("case_id"),
        "case_submitter_id": case.get("submitter_id"),
        "sample_id": sample.get("sample_id"),
        "sample_submitter_id": sample.get("submitter_id"),
        "sample_type": sample.get("sample_type"),
        "aliquot_id": aliquot_id,
        "aliquot_submitter_id": _aliquot_to_submitter(case).get(aliquot_id or ""),
    }


def _aliquot_entities(entry: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        e
        for e in (entry.get("associated_entities") or [])
        if e.get("entity_type") == "aliquot" and e.get("entity_id")
    ]


def _submitter_for_entity(entry: dict[str, Any], entity_id: str | None) -> str | None:
    for e in _aliquot_entities(entry):
        if e["entity_id"] == entity_id:
            return e.get("entity_submitter_id")
    return None


def _seg_aliquot(df: pd.DataFrame) -> str | None:
    """The single aliquot a seg file's own `GDC_Aliquot` column reports."""
    values = {v for v in df["GDC_Aliquot"].dropna().unique()}
    return next(iter(values)) if len(values) == 1 else None


def _allele_specific_copy_number_segment_rows(
    cases: list[dict[str, Any]],
    project_raw_dir: Path,
) -> list[dict[str, Any]]:
    """One row per allele-specific copy number segment.

    ASCAT is a paired caller, so each file names two aliquots. The tumour is
    identified positively as the aliquot the file's own `GDC_Aliquot` column
    reports; whichever other aliquot the file is associated with is the
    matched normal. Files whose `GDC_Aliquot` isn't a single consistent value
    are skipped rather than guessed at.
    """
    rows: list[dict[str, Any]] = []
    for entry, file_path, case in _iter_modality_files(
        project_raw_dir, "copy_number_allele_specific", cases
    ):
        df = pd.read_csv(file_path, sep="\t", low_memory=False)
        if df.empty:
            continue
        tumor_aliquot = _seg_aliquot(df)
        if tumor_aliquot is None:
            continue
        others = [
            e["entity_id"] for e in _aliquot_entities(entry) if e["entity_id"] != tumor_aliquot
        ]
        normal_aliquot = others[0] if len(others) == 1 else None
        common = {
            **_aliquot_fks(case, tumor_aliquot),
            "matched_normal_aliquot_id": normal_aliquot,
            "matched_normal_aliquot_submitter_id": _submitter_for_entity(entry, normal_aliquot),
            "workflow_type": entry.get("workflow_type"),
            "experimental_strategy": entry.get("experimental_strategy"),
            "source_file_id": entry["file_id"],
        }
        for r in df.itertuples(index=False):
            rows.append(
                {
                    **common,
                    "chromosome": _na_to_none(r.Chromosome),
                    "start": _to_int_or_none(r.Start),
                    "end": _to_int_or_none(r.End),
                    "copy_number": _to_int_or_none(r.Copy_Number),
                    "major_copy_number": _to_int_or_none(r.Major_Copy_Number),
                    "minor_copy_number": _to_int_or_none(r.Minor_Copy_Number),
                }
            )
    return rows


def _masked_copy_number_segment_rows(
    cases: list[dict[str, Any]],
    project_raw_dir: Path,
) -> list[dict[str, Any]]:
    """One row per masked (germline-CNV-removed) DNAcopy segment.

    Single-aliquot files, so the FK comes straight from `GDC_Aliquot`.
    Note this data type writes bare chromosome names ("1") where the
    allele-specific type writes "chr1"; both are carried as written.
    """
    rows: list[dict[str, Any]] = []
    for entry, file_path, case in _iter_modality_files(
        project_raw_dir, "copy_number_masked", cases
    ):
        df = pd.read_csv(file_path, sep="\t", low_memory=False)
        if df.empty:
            continue
        aliquot_id = _seg_aliquot(df)
        if aliquot_id is None:
            continue
        common = {
            **_aliquot_fks(case, aliquot_id),
            "workflow_type": entry.get("workflow_type"),
            "source_file_id": entry["file_id"],
        }
        for r in df.itertuples(index=False):
            rows.append(
                {
                    **common,
                    "chromosome": _na_to_none(str(r.Chromosome)),
                    "start": _to_int_or_none(r.Start),
                    "end": _to_int_or_none(r.End),
                    "num_probes": _to_int_or_none(r.Num_Probes),
                    "segment_mean": _to_float_or_none(r.Segment_Mean),
                }
            )
    return rows


def _seg_aliquot_barcode(df: pd.DataFrame) -> str | None:
    """The single aliquot *barcode* a GATK4 CNV seg file reports.

    GATK4 CNV is the one segment data type that names its aliquot by
    submitter id rather than UUID, under `GDC_Aliquot_ID` rather than
    `GDC_Aliquot`.
    """
    values = {v for v in df["GDC_Aliquot_ID"].dropna().unique()}
    return next(iter(values)) if len(values) == 1 else None


def _copy_number_segment_rows(
    cases: list[dict[str, Any]],
    project_raw_dir: Path,
) -> list[dict[str, Any]]:
    """One row per unmasked copy number segment, from both workflows.

    The two workflows are read from separate raw directories because they
    write incompatible headers, and are unioned into one table because they
    carry the same six measurements:

      - **DNAcopy** is single-aliquot with a UUID in `GDC_Aliquot`, exactly
        like the masked type it is the unmasked counterpart of.
      - **GATK4 CNV** is a paired WGS caller with a *barcode* in
        `GDC_Aliquot_ID`. The tumour is resolved positively by matching that
        barcode against the file's `associated_entities.entity_submitter_id`
        — never by reading sample-type digits out of it — and whichever
        other aliquot the file names is the matched normal.

    Files whose aliquot column isn't a single consistent value, or whose
    barcode doesn't match an associated entity, are skipped rather than
    guessed at.
    """
    rows: list[dict[str, Any]] = []
    for entry, file_path, case in _iter_modality_files(
        project_raw_dir, "copy_number_segment_dnacopy", cases
    ):
        df = pd.read_csv(file_path, sep="\t", low_memory=False)
        if df.empty:
            continue
        aliquot_id = _seg_aliquot(df)
        if aliquot_id is None:
            continue
        common = {
            **_aliquot_fks(case, aliquot_id),
            "matched_normal_aliquot_id": None,
            "matched_normal_aliquot_submitter_id": None,
            "workflow_type": entry.get("workflow_type"),
            "experimental_strategy": entry.get("experimental_strategy"),
            "source_file_id": entry["file_id"],
        }
        rows.extend(_segment_measurement_rows(df, common))

    for entry, file_path, case in _iter_modality_files(
        project_raw_dir, "copy_number_segment_gatk4", cases
    ):
        df = pd.read_csv(file_path, sep="\t", low_memory=False)
        if df.empty:
            continue
        barcode = _seg_aliquot_barcode(df)
        if barcode is None:
            continue
        entities = _aliquot_entities(entry)
        tumor = next(
            (e["entity_id"] for e in entities if e.get("entity_submitter_id") == barcode), None
        )
        if tumor is None:
            continue
        others = [e["entity_id"] for e in entities if e["entity_id"] != tumor]
        normal = others[0] if len(others) == 1 else None
        common = {
            **_aliquot_fks(case, tumor),
            "matched_normal_aliquot_id": normal,
            "matched_normal_aliquot_submitter_id": _submitter_for_entity(entry, normal),
            "workflow_type": entry.get("workflow_type"),
            "experimental_strategy": entry.get("experimental_strategy"),
            "source_file_id": entry["file_id"],
        }
        rows.extend(_segment_measurement_rows(df, common))
    return rows


def _segment_measurement_rows(
    df: pd.DataFrame, common: dict[str, Any]
) -> list[dict[str, Any]]:
    """The per-segment columns shared by both unmasked workflows."""
    return [
        {
            **common,
            "chromosome": _na_to_none(str(r.Chromosome)),
            "start": _to_int_or_none(r.Start),
            "end": _to_int_or_none(r.End),
            "num_probes": _to_int_or_none(r.Num_Probes),
            "segment_mean": _to_float_or_none(r.Segment_Mean),
        }
        for r in df.itertuples(index=False)
    ]


def _ascat_tumor_by_pair(project_raw_dir: Path) -> dict[tuple[str, frozenset[str]], str]:
    """Map (workflow_type, {aliquot ids}) -> tumour aliquot, from the seg files.

    The three ASCAT gene-level workflows are paired and name two aliquots,
    but — unlike every segment data type — their TSVs carry no aliquot
    column at all, so the tumour cannot be read out of the file. Filenames
    don't settle it either: ASCAT2 and ASCAT3 embed the tumour aliquot UUID,
    but AscatNGS embeds a UUID that is not an aliquot at all.

    What does settle it is that each paired gene-level file has exactly one
    counterpart in `allele_specific_copy_number_segment` with the same
    workflow and the same pair of associated aliquots, and *that* file names
    its tumour positively in `GDC_Aliquot`. Verified on TCGA-CHOL: 114 of
    114 paired gene-level files resolve, none ambiguous.

    Returns an empty map when the allele-specific directory is absent, which
    makes the paired gene-level files unresolvable and therefore skipped —
    the same "skip rather than guess" contract the segment loaders honour.
    """
    mod_dir = project_raw_dir / "copy_number_allele_specific"
    manifest_path = mod_dir / "manifest.json"
    if not manifest_path.exists():
        return {}
    out: dict[tuple[str, frozenset[str]], str] = {}
    for entry in json.loads(manifest_path.read_text()):
        file_path = mod_dir / entry["file_name"]
        if not file_path.exists():
            continue
        entities = frozenset(e["entity_id"] for e in _aliquot_entities(entry))
        if len(entities) != 2:
            continue
        df = pd.read_csv(file_path, sep="\t", usecols=["GDC_Aliquot"], low_memory=False)
        tumor = _seg_aliquot(df)
        if tumor is None or tumor not in entities:
            continue
        out[(entry.get("workflow_type") or "", entities)] = tumor
    return out


# Columns we read out of a gene-level TSV. The other four (`gene_name`,
# `chromosome`, `start`, `end`) are identical in every file and are emitted
# once into `gene_model` instead.
_GENE_LEVEL_VALUE_COLUMNS = ("copy_number", "min_copy_number", "max_copy_number")


def _gene_level_copy_number_batches(
    cases: list[dict[str, Any]],
    project_raw_dir: Path,
) -> Any:
    """Yield one `pa.RecordBatch` per gene-level copy number file.

    This table is emitted as a stream of Arrow batches rather than a list of
    dicts because of its size: TCGA-BRCA alone is 3,314 files x 60,623 genes
    = 201M rows, which as Python dicts would be roughly 100 GB of interpreter
    objects. One batch per file caps peak memory at ~60k rows regardless of
    project size, and `write_tables` streams the batches straight to Parquet.

    Tumour resolution differs by workflow: ABSOLUTE LiftOver files name one
    aliquot and need none; the three ASCAT workflows are paired and resolve
    through `_ascat_tumor_by_pair`.
    """
    import pyarrow as pa

    schema = TABULAR_TABLES["gene_level_copy_number"]
    tumor_by_pair = _ascat_tumor_by_pair(project_raw_dir)
    for entry, file_path, case in _iter_modality_files(
        project_raw_dir, "gene_level_copy_number", cases
    ):
        entities = _aliquot_entities(entry)
        if len(entities) == 1:
            tumor: str | None = entities[0]["entity_id"]
            normal: str | None = None
        else:
            key = (entry.get("workflow_type") or "", frozenset(e["entity_id"] for e in entities))
            tumor = tumor_by_pair.get(key)
            if tumor is None:
                continue
            others = [e["entity_id"] for e in entities if e["entity_id"] != tumor]
            normal = others[0] if len(others) == 1 else None
        df = pd.read_csv(
            file_path,
            sep="\t",
            usecols=["gene_id", *_GENE_LEVEL_VALUE_COLUMNS],
            low_memory=False,
        )
        if df.empty:
            continue
        n = len(df)
        common = {
            **_aliquot_fks(case, tumor),
            "matched_normal_aliquot_id": normal,
            "matched_normal_aliquot_submitter_id": _submitter_for_entity(entry, normal),
            "workflow_type": entry.get("workflow_type"),
            "experimental_strategy": entry.get("experimental_strategy"),
            "source_file_id": entry["file_id"],
        }
        arrays = []
        for field in schema:
            if field.name in common:
                arrays.append(pa.array([common[field.name]] * n, type=field.type))
            elif field.name == "gene_id":
                arrays.append(pa.array(df["gene_id"].astype("string"), type=field.type))
            else:
                # Nullable ints: the callers leave genes uncalled (all of chrM,
                # plus anything outside a segment), and pandas reads those as NaN.
                arrays.append(pa.array(df[field.name], type=field.type, from_pandas=True))
        yield pa.RecordBatch.from_arrays(arrays, schema=schema)


def _methylation_beta_value_batches(
    cases: list[dict[str, Any]],
    project_raw_dir: Path,
) -> Any:
    """Yield one `pa.RecordBatch` per aliquot's SeSAMe beta TXT.

    Streamed like the other per-probe/per-gene tables: a 450k file is
    486,427 rows, so even TCGA-CHOL's 45 files are 21.9M rows.

    The source file is **headerless** — two unnamed columns, probe id then
    beta — so the names are supplied here rather than read. Masked probes
    are written by SeSAMe as `NA` and land as null, which is a real value
    meaning "not trustworthy", not zero.
    """
    import pyarrow as pa

    schema = TABULAR_TABLES["methylation_beta_value"]
    for entry, file_path, case in _iter_modality_files(project_raw_dir, "methylation", cases):
        entities = _aliquot_entities(entry)
        if len(entities) != 1:
            continue
        aliquot_id = entities[0]["entity_id"]
        df = pd.read_csv(
            file_path,
            sep="\t",
            header=None,
            names=["probe_id", "beta_value"],
            low_memory=False,
        )
        if df.empty:
            continue
        n = len(df)
        common = {
            **_aliquot_fks(case, aliquot_id),
            "platform": entry.get("platform"),
            "source_file_id": entry["file_id"],
        }
        arrays = []
        for field in schema:
            if field.name in common:
                arrays.append(pa.array([common[field.name]] * n, type=field.type))
            elif field.name == "probe_id":
                arrays.append(pa.array(df["probe_id"].astype("string"), type=field.type))
            else:
                arrays.append(pa.array(df["beta_value"], type=field.type, from_pandas=True))
        yield pa.RecordBatch.from_arrays(arrays, schema=schema)


def _isoform_expression_quantification_rows(
    cases: list[dict[str, Any]],
    project_raw_dir: Path,
) -> list[dict[str, Any]]:
    """One row per (aliquot, miRNA isoform).

    Single-aliquot files with a header. `cross-mapped` is underscored to
    `cross_mapped` for the same reason it is in the miRNA table: the hyphen
    is not a legal bare SQL identifier.
    """
    rows: list[dict[str, Any]] = []
    for entry, file_path, case in _iter_modality_files(project_raw_dir, "mirna_isoform", cases):
        entities = _aliquot_entities(entry)
        if len(entities) != 1:
            continue
        aliquot_id = entities[0]["entity_id"]
        df = pd.read_csv(file_path, sep="\t", low_memory=False)
        if df.empty:
            continue
        common = {
            **_aliquot_fks(case, aliquot_id),
            "source_file_id": entry["file_id"],
        }
        for r in df.itertuples(index=False):
            rows.append(
                {
                    **common,
                    "mirna_id": _na_to_none(r.miRNA_ID),
                    "isoform_coords": _na_to_none(r.isoform_coords),
                    "read_count": _to_int_or_none(r.read_count),
                    "reads_per_million_mirna_mapped": _to_float_or_none(
                        r.reads_per_million_miRNA_mapped
                    ),
                    # `cross-mapped` isn't a valid attribute name on the
                    # namedtuple pandas yields, so it comes back positionally.
                    "cross_mapped": _na_to_none(r[df.columns.get_loc("cross-mapped")]),
                    "mirna_region": _na_to_none(r.miRNA_region),
                }
            )
    return rows


def _gene_model_rows(project_raw_dir: Path) -> list[dict[str, Any]]:
    """The GENCODE v36 gene model, assembled from the two GDC sources.

    `gene_name` / `gene_type` come from a STAR expression TSV and
    `chromosome` / `start` / `end` from a gene-level copy number TSV. Both
    halves are byte-identical across every file of their kind (the gene
    model is fixed for v36), so one file of each is read rather than all of
    them; `test_gene_model_is_invariant_across_files` guards that.

    Genes present in only one source keep the half that exists: the 37 chrM
    genes are expression-only and carry null coordinates rather than
    coordinates borrowed from outside the GDC.
    """
    by_gene: dict[str, dict[str, Any]] = {}

    expr_file = _first_modality_file(project_raw_dir, "expression")
    if expr_file is not None:
        df = pd.read_csv(
            expr_file, sep="\t", comment="#", usecols=["gene_id", "gene_name", "gene_type"]
        )
        df = df[~df["gene_id"].str.startswith("N_", na=False)]
        for r in df.itertuples(index=False):
            by_gene[r.gene_id] = {
                "gene_id": _na_to_none(r.gene_id),
                "gene_name": _na_to_none(r.gene_name),
                "gene_type": _na_to_none(r.gene_type),
                "chromosome": None,
                "start": None,
                "end": None,
            }

    gl_file = _first_modality_file(project_raw_dir, "gene_level_copy_number")
    if gl_file is not None:
        df = pd.read_csv(
            gl_file, sep="\t", usecols=["gene_id", "gene_name", "chromosome", "start", "end"]
        )
        for r in df.itertuples(index=False):
            row = by_gene.setdefault(
                r.gene_id,
                {
                    "gene_id": _na_to_none(r.gene_id),
                    "gene_name": _na_to_none(r.gene_name),
                    "gene_type": None,
                    "chromosome": None,
                    "start": None,
                    "end": None,
                },
            )
            row["chromosome"] = _na_to_none(r.chromosome)
            row["start"] = _to_int_or_none(r.start)
            row["end"] = _to_int_or_none(r.end)

    return sorted(by_gene.values(), key=lambda r: r["gene_id"] or "")


def _first_modality_file(project_raw_dir: Path, modality_dir: str) -> Path | None:
    """The first downloaded file in a modality dir, in manifest order."""
    mod_dir = project_raw_dir / modality_dir
    manifest_path = mod_dir / "manifest.json"
    if not manifest_path.exists():
        return None
    for entry in json.loads(manifest_path.read_text()):
        file_path = mod_dir / entry["file_name"]
        if file_path.exists():
            return file_path
    return None


def _mirna_expression_quantification_rows(
    cases: list[dict[str, Any]],
    project_raw_dir: Path,
) -> list[dict[str, Any]]:
    """One row per (aliquot, miRBase v21 mature miRNA).

    miRNA files carry no aliquot column of their own, so the FK is the file's
    single associated aliquot entity.
    """
    rows: list[dict[str, Any]] = []
    for entry, file_path, case in _iter_modality_files(project_raw_dir, "mirna", cases):
        entities = _aliquot_entities(entry)
        if len(entities) != 1:
            continue
        aliquot_id = entities[0]["entity_id"]
        df = pd.read_csv(file_path, sep="\t", low_memory=False)
        # Source header is `cross-mapped`. The hyphen makes it unusable as a
        # bare SQL identifier — and unreachable as an attribute on an
        # `itertuples` row — so rename it here. Values are untouched.
        df = df.rename(columns={"cross-mapped": "cross_mapped"})
        common = {**_aliquot_fks(case, aliquot_id), "source_file_id": entry["file_id"]}
        for r in df.itertuples(index=False):
            rows.append(
                {
                    **common,
                    "mirna_id": _na_to_none(r.miRNA_ID),
                    "read_count": _to_int_or_none(r.read_count),
                    "reads_per_million_mirna_mapped": _to_float_or_none(
                        r.reads_per_million_miRNA_mapped
                    ),
                    "cross_mapped": _na_to_none(r.cross_mapped),
                }
            )
    return rows


def _portion_to_sample(case: dict[str, Any]) -> dict[str, dict[str, Any]]:
    """Map portion_id -> {portion_submitter_id, sample dict} for one case."""
    out: dict[str, dict[str, Any]] = {}
    for s in case.get("samples") or []:
        for p in s.get("portions") or []:
            pid = p.get("portion_id")
            if pid:
                out[pid] = {"portion_submitter_id": p.get("submitter_id"), "sample": s}
    return out


def _protein_expression_quantification_rows(
    cases: list[dict[str, Any]],
    project_raw_dir: Path,
) -> list[dict[str, Any]]:
    """One row per (portion, antibody) RPPA measurement.

    RPPA is the one modality that attaches to a portion rather than an
    aliquot, so the sample FK is resolved through the portion.
    """
    rows: list[dict[str, Any]] = []
    for entry, file_path, case in _iter_modality_files(
        project_raw_dir, "protein_expression", cases
    ):
        portions = [
            e
            for e in (entry.get("associated_entities") or [])
            if e.get("entity_type") == "portion" and e.get("entity_id")
        ]
        if len(portions) != 1:
            continue
        portion_id = portions[0]["entity_id"]
        resolved = _portion_to_sample(case).get(portion_id, {})
        sample = resolved.get("sample") or {}
        df = pd.read_csv(file_path, sep="\t", low_memory=False)
        common = {
            "case_id": case.get("case_id"),
            "case_submitter_id": case.get("submitter_id"),
            "sample_id": sample.get("sample_id"),
            "sample_submitter_id": sample.get("submitter_id"),
            "sample_type": sample.get("sample_type"),
            "portion_id": portion_id,
            "portion_submitter_id": (
                resolved.get("portion_submitter_id") or portions[0].get("entity_submitter_id")
            ),
            "source_file_id": entry["file_id"],
        }
        for r in df.itertuples(index=False):
            rows.append(
                {
                    **common,
                    "agid": _na_to_none(r.AGID),
                    # Numeric in the source but an identifier, not a
                    # measurement; kept as a string so it never picks up a
                    # decimal point on a null-containing column.
                    "lab_id": None if pd.isna(r.lab_id) else str(r.lab_id),
                    "catalog_number": _na_to_none(r.catalog_number),
                    "set_id": _na_to_none(r.set_id),
                    "peptide_target": _na_to_none(r.peptide_target),
                    "protein_expression": _to_float_or_none(r.protein_expression),
                }
            )
    return rows


# ---------------------------------------------------------------------------
# ssGSEA pathway activity
# ---------------------------------------------------------------------------


def _ssgsea_scores_rows(
    cases: list[dict[str, Any]],
    project_raw_dir: Path,
    collection: str,
    msigdb_dir: Path,
    load_matrix: Any = None,
) -> list[dict[str, Any]]:
    """One row per (aliquot, pathway) of raw ssGSEA scores for one project.

    Depends on nothing outside this project: raw scores are a function of a
    sample's own expression, the gene universe and the gene set. That is
    what lets pathway activity be built inside the ordinary per-project
    loop; the cohort-level reference distributions are a separate pass (see
    `ssgsea_stats_rows`).
    """
    gmt_path = msigdb_dir / _msigdb_mod.COLLECTIONS[collection].file_name
    if not gmt_path.exists():
        return []
    # Reading every STAR TSV dominates the cost, and the matrix is the same
    # for every collection, so the caller passes a memoized loader rather
    # than paying for it once per collection.
    loader = load_matrix or (lambda: _expression_mod.tpm_matrix_for_project(project_raw_dir))
    genes, matrix, records = loader()
    if not records:
        return []

    kept, stats = _ssgsea_mod.map_gene_sets(_ssgsea_mod.load_gmt(gmt_path), genes)
    if not kept:
        return []
    names = list(kept)
    by_name = {s["pathway"]: s for s in stats}
    scores = _ssgsea_mod.ssgsea_scores(matrix.astype(np.float64), [kept[n] for n in names])

    case_index = {c.get("case_id"): c for c in cases}
    rows: list[dict[str, Any]] = []
    for j, rec in enumerate(records):
        case = case_index.get(rec["case_id"])
        if case is None:
            continue
        a2s = _aliquot_to_sample(case)
        a2sub = _aliquot_to_submitter(case)
        sample_id = a2s.get(rec["aliquot_id"])
        sample = next(
            (s for s in (case.get("samples") or []) if s.get("sample_id") == sample_id), {}
        )
        common = {
            "case_id": rec["case_id"],
            "case_submitter_id": case.get("submitter_id"),
            "sample_id": sample_id,
            "sample_submitter_id": sample.get("submitter_id"),
            "sample_type": sample.get("sample_type"),
            "aliquot_id": rec["aliquot_id"],
            "aliquot_submitter_id": a2sub.get(rec["aliquot_id"]),
            "source_file_id": rec["source_file_id"],
        }
        for i, name in enumerate(names):
            rows.append(
                {
                    **common,
                    "pathway": name,
                    "pathway_url": _msigdb_mod.geneset_url(name),
                    "matched_gene_count": int(by_name[name]["matched_gene_count"]),
                    "original_gene_count": int(by_name[name]["original_gene_count"]),
                    "score_raw": float(scores[i, j]),
                }
            )
    return rows


def _describe(series: pd.Series) -> dict[str, Any]:
    """Reference distribution for one (population, pathway) cell."""
    n = int(series.size)
    return {
        "n_aliquots": n,
        "mean": float(series.mean()),
        # ddof=1 is undefined for a single observation; emit null rather
        # than a misleading 0.0.
        "sd": float(series.std(ddof=1)) if n > 1 else None,
        "min": float(series.min()),
        "q25": float(series.quantile(0.25)),
        "median": float(series.median()),
        "q75": float(series.quantile(0.75)),
        "max": float(series.max()),
    }


def ssgsea_stats_rows(
    processed_dir: Path,
    collection: str,
    sample_types: tuple[str, ...] = ("Primary Tumor",),
) -> dict[str, list[dict[str, Any]]]:
    """Reference distributions per project, computed from the written scores.

    Returns `{project_id: rows}`. Every project's rows carry that project's
    own distributions **and** the pan-cancer ones, so a consumer who loads a
    single project can still normalize against all of TCGA without scanning
    every config.

    This reads the `ssgsea_scores_*` parquets rather than rescoring: the
    stats are pure aggregates of published values, and computing them from
    the published bytes is the strongest guarantee that they agree. It also
    means `--table ssgsea_stats_<collection>` can refresh them alone.

    Everything here is composition-dependent — adding a project changes the
    pan-cancer rows — which is precisely why it lives in this small table
    instead of as extra columns on the immutable scores.
    """
    import pyarrow.parquet as pq

    table_name = _schema_mod.ssgsea_scores_table(collection)
    frames: dict[str, pd.DataFrame] = {}
    for path in sorted(processed_dir.glob(f"*/{table_name}/data.parquet")):
        df = pq.read_table(
            path, columns=["pathway", "sample_type", "score_raw"]
        ).to_pandas()
        if len(df):
            frames[path.parent.parent.name] = df
    if not frames:
        return {}

    def _rows_for(df: pd.DataFrame, population: str, project_id: str | None) -> list[dict]:
        out = []
        for sample_type in (None, *sample_types):
            sub = df if sample_type is None else df[df["sample_type"] == sample_type]
            if not len(sub):
                continue
            for pathway, series in sub.groupby("pathway")["score_raw"]:
                out.append(
                    {
                        "population": population,
                        "project_id": project_id,
                        "sample_type": sample_type,
                        "pathway": pathway,
                        **_describe(series),
                    }
                )
        return out

    pan = _rows_for(pd.concat(frames.values(), ignore_index=True), "pan_cancer", None)
    return {proj: _rows_for(df, "project", proj) + pan for proj, df in frames.items()}


# ---------------------------------------------------------------------------
# Files (provenance from per-modality manifests)
# ---------------------------------------------------------------------------


# Which published table carries a given raw modality directory's content.
#
# Not every fetched modality is published: the four BCR XML directories are
# downloaded so the `files` table can describe them and so the bytes are
# local, but they are a second serialization of data the biotab tables
# already carry (verified value-for-value — 918/918 fields agree between
# `bcr ssf xml` and `ssf_tumor_samples`), so they map to None.
#
# A `*` suffix means a family of per-form tables rather than one table.
MODALITY_TABLE: dict[str, str | None] = {
    "mutations": "masked_somatic_mutation",
    "expression": "gene_expression_quantification",
    "mirna": "mirna_expression_quantification",
    "mirna_isoform": "isoform_expression_quantification",
    "protein_expression": "protein_expression_quantification",
    "methylation": "methylation_beta_value",
    "pathology_reports": "pathology_report",
    "copy_number_allele_specific": "allele_specific_copy_number_segment",
    "copy_number_masked": "masked_copy_number_segment",
    "copy_number_segment_dnacopy": "copy_number_segment",
    "copy_number_segment_gatk4": "copy_number_segment",
    "gene_level_copy_number": "gene_level_copy_number",
    "clinical_supplement": "clinical_supplement_*",
    "biospecimen_supplement": "biospecimen_supplement_*",
    "clinical_supplement_xml": None,
    "clinical_supplement_omf_xml": None,
    "biospecimen_supplement_xml": None,
    "biospecimen_supplement_ssf_xml": None,
}

_GDC_DATA_URL = "https://api.gdc.cancer.gov/data/{file_id}"


def _local_file_entries(project_raw_dir: Path) -> dict[str, dict[str, Any]]:
    """{file_id: manifest entry + modality} for everything downloaded locally."""
    out: dict[str, dict[str, Any]] = {}
    for modality_dir in sorted(p for p in project_raw_dir.iterdir() if p.is_dir()):
        manifest_path = modality_dir / "manifest.json"
        if not manifest_path.exists():
            continue
        for entry in json.loads(manifest_path.read_text()):
            if entry.get("_status") == "manifest_only":
                # Discovered but its bytes were never downloaded (a capped
                # fetch), so nothing here can be in a published table.
                continue
            out[entry["file_id"]] = {**entry, "_modality": modality_dir.name}
    return out


def _files_rows(project_raw_dir: Path) -> list[dict[str, Any]]:
    """One row per open-access GDC file for the project.

    The spine is `files_index.json` — every open file the GDC holds, whether
    or not this pipeline downloaded it — so the table describes the
    project's entire open footprint. Local manifests then enrich the rows
    for files we actually have, supplying `modality` and the fields the
    index doesn't carry.

    Falls back to local manifests alone when no index has been fetched, so
    a project built before `fetch-file-index` existed still produces a
    table; those rows simply describe less.

    One row per file. A file that names several cases (the project-level
    BCR biotabs name all of them) reports a null `case_id` rather than being
    exploded into one row per case, which would imply a per-patient file
    that does not exist.
    """
    rows: list[dict[str, Any]] = []
    if not project_raw_dir.exists():
        return rows
    local = _local_file_entries(project_raw_dir)

    index_path = project_raw_dir / "files_index.json"
    if index_path.exists():
        spine = json.loads(index_path.read_text())
    else:
        spine = [{**e, "cases": e.get("cases") or []} for e in local.values()]

    for entry in spine:
        file_id = entry["file_id"]
        got = local.get(file_id)
        modality = got["_modality"] if got else None
        table = MODALITY_TABLE.get(modality) if modality else None
        cases = entry.get("cases") or []
        case = cases[0] if len(cases) == 1 else {}
        # Prefer the index's values (fetched in one consistent query) and
        # fall back to the manifest for anything it lacks.
        src = {**(got or {}), **{k: v for k, v in entry.items() if v is not None}}
        rows.append(
            {
                "case_id": case.get("case_id"),
                "case_submitter_id": case.get("submitter_id"),
                "file_id": file_id,
                "file_name": src.get("file_name"),
                "file_size": src.get("file_size"),
                "md5sum": src.get("md5sum"),
                "data_category": src.get("data_category"),
                "data_type": src.get("data_type"),
                "data_format": src.get("data_format"),
                "experimental_strategy": src.get("experimental_strategy"),
                "workflow_type": src.get("workflow_type"),
                "access": src.get("access"),
                "gdc_version": src.get("gdc_version"),
                "gdc_first_release": src.get("gdc_first_release"),
                "gdc_superseded": src.get("gdc_superseded"),
                "platform": src.get("platform"),
                "modality": modality,
                "in_dataset": table is not None,
                "dataset_table": table,
                "gdc_download_url": _GDC_DATA_URL.format(file_id=file_id),
            }
        )
    rows.sort(key=lambda r: (r.get("data_type") or "", r.get("file_name") or ""))
    return rows


# ---------------------------------------------------------------------------
# Top-level builder
# ---------------------------------------------------------------------------


def build_tables(
    cases: list[dict[str, Any]],
    project_raw_dir: Path,
    only: set[str] | None = None,
    msigdb_dir: Path | None = None,
) -> dict[str, list[dict[str, Any]]]:
    """Build tabular tables for one project (all of them unless `only` is given).

    Returns one entry per table keyed by table name. The fixed-schema
    keys (`cases`, `masked_somatic_mutation`, `gene_expression_quantification`,
    `pathology_report`, `files`, `survival_derived`) match `TABULAR_TABLES`. The
    `clinical_supplement_*` keys are flex-schema tables — their parquet
    column set is inferred per project, since BCR biotab forms vary by
    cancer type (e.g. BLCA's clinical_patient form has bladder-specific
    fields that don't exist in CHOL's). Cross-project queries union via
    HF `concatenate_datasets` with NULL padding — the same pattern
    cBioPortal and the GDC use for their per-study clinical exports.

    `only` restricts the build to the named tables. Each table's rows are
    produced by a thunk that is called only when its table is wanted, so
    the expensive emitters — `gene_expression_quantification` re-reads
    every STAR TSV in the project, `masked_somatic_mutation` every MAF —
    cost nothing when a caller is appending one new table to an existing
    tree. `cases` is still built whenever `survival_derived` is requested,
    since the derived endpoints are projected off those rows.

    Note: the `survival_derived` table starts empty here. The caller must
    `survival.attach_survival(tables["cases"])` and then call
    `derived_survival_rows(tables["cases"])` to populate it. Splitting
    that off keeps `build_tables` purely deterministic on raw inputs;
    survival re-derivation depends on Clinical Supplement data the caller
    may attach beforehand.
    """

    def _mutation_rows() -> list[dict[str, Any]]:
        return _mutations_rows(cases, _mutations_mod.load_for_project(project_raw_dir))

    thunks: dict[str, Any] = {
        "cases": lambda: _cases_rows(cases),
        "masked_somatic_mutation": _mutation_rows,
        # Yields Arrow batches rather than dicts; see the emitter's docstring.
        "gene_expression_quantification": lambda: _expression_batches(cases, project_raw_dir),
        "files": lambda: _files_rows(project_raw_dir),
        "survival_derived": list,
        "pathology_report": lambda: _pathology_report_rows(cases, project_raw_dir),
        "allele_specific_copy_number_segment": lambda: (
            _allele_specific_copy_number_segment_rows(cases, project_raw_dir)
        ),
        "masked_copy_number_segment": lambda: (
            _masked_copy_number_segment_rows(cases, project_raw_dir)
        ),
        "copy_number_segment": lambda: _copy_number_segment_rows(cases, project_raw_dir),
        # Yields Arrow batches rather than dicts; see the emitter's docstring.
        "gene_level_copy_number": lambda: _gene_level_copy_number_batches(
            cases, project_raw_dir
        ),
        "gene_model": lambda: _gene_model_rows(project_raw_dir),
        # Yields Arrow batches rather than dicts; see the emitter's docstring.
        "methylation_beta_value": lambda: _methylation_beta_value_batches(
            cases, project_raw_dir
        ),
        "isoform_expression_quantification": lambda: (
            _isoform_expression_quantification_rows(cases, project_raw_dir)
        ),
        "mirna_expression_quantification": lambda: (
            _mirna_expression_quantification_rows(cases, project_raw_dir)
        ),
        "protein_expression_quantification": lambda: (
            _protein_expression_quantification_rows(cases, project_raw_dir)
        ),
    }
    # ssGSEA scores are per-project pure; the matching stats tables are a
    # cohort-level aggregate written by the caller after every project's
    # scores exist (see `ssgsea_stats_rows`), so they start empty here for
    # the same reason `survival_derived` does.
    msigdb = msigdb_dir if msigdb_dir is not None else project_raw_dir.parent / "msigdb"
    _matrix_cache: list[Any] = []

    def _matrix():
        if not _matrix_cache:
            _matrix_cache.append(_expression_mod.tpm_matrix_for_project(project_raw_dir))
        return _matrix_cache[0]

    for coll in SSGSEA_COLLECTIONS:
        thunks[_schema_mod.ssgsea_scores_table(coll)] = (
            lambda c=coll: _ssgsea_scores_rows(cases, project_raw_dir, c, msigdb, _matrix)
        )
        thunks[_schema_mod.ssgsea_stats_table(coll)] = list

    wanted = set(thunks) if only is None else set(only)
    # survival_derived is projected off the cases rows, so it drags `cases`
    # in even when the caller didn't ask for it directly.
    if "survival_derived" in wanted:
        wanted.add("cases")

    tables: dict[str, list[dict[str, Any]]] = {
        name: thunk() for name, thunk in thunks.items() if name in wanted
    }

    if only is None or any(t.startswith("clinical_supplement_") for t in wanted):
        supp_rows = _clinical_supplement_mod.build_tabular_rows(
            project_raw_dir / "clinical_supplement"
        )
        for suffix, rows in supp_rows.items():
            name = f"clinical_supplement_{suffix}"
            if only is None or name in wanted:
                tables[name] = rows

    if only is None or any(t.startswith("biospecimen_supplement_") for t in wanted):
        bio_rows = _biospecimen_supplement_mod.build_tabular_rows(
            project_raw_dir / "biospecimen_supplement"
        )
        for suffix, rows in bio_rows.items():
            name = f"biospecimen_supplement_{suffix}"
            if only is None or name in wanted:
                tables[name] = rows
    return tables


def derived_survival_rows(cases: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Project the `survival_derived` struct off each case row into a flat
    table. Caller must run `survival.attach_survival(cases)` first.

    Output: one row per case with `case_submitter_id` plus the 8 derived
    survival fields (os_event, os_time, dss_event, dss_time, pfi_event,
    pfi_time, dfi_event, dfi_time). Rows whose `survival_derived` is None
    or all-null are skipped — they have nothing to contribute.
    """
    out: list[dict[str, Any]] = []
    for r in cases:
        sd = r.get("survival_derived") or {}
        if all(sd.get(k) is None for k in (
            "os_event", "os_time", "dss_event", "dss_time",
            "pfi_event", "pfi_time", "dfi_event", "dfi_time",
        )):
            continue
        out.append({"case_submitter_id": r["case_submitter_id"], **sd})
    return out


# ---------------------------------------------------------------------------
# Parquet writer
# ---------------------------------------------------------------------------

# Smaller row groups for the "expression" table since each row is narrow
# (~12 cols of small scalars) but the count is huge — 60660 rows per
# aliquot × 600 aliquots ≈ 36M rows for LUAD. 100k rows per group keeps
# each group well under the HF Data Studio 300 MB scan limit and gives
# the SQL viewer fast random access. The `cases` table's rows are wide
# (full nested case JSON); 50 rows per group there matches the
# consolidated build's sizing.
# A row group is the unit of column-chunk metadata, dictionary pages and
# page-index entries, so its size should follow how wide a row is. 50 rows
# is right when one row is a whole nested case record or an embedded PDF;
# it is pathological for a narrow measurement row, where it means thousands
# of tiny groups each paying that fixed overhead in full. TCGA-CHOL's
# `isoform_expression_quantification` was 3,897 groups over 194,815 rows.
#
# So the default is sized for the narrow long tables that make up nearly
# every entry in `TABULAR_TABLES`, and the genuinely wide ones are pinned
# small below.
_ROW_GROUP_SIZE_DEFAULT = 100_000
# Row groups for the per-gene tables must comfortably exceed the gene
# cardinality, not just "be large". Parquet writes a fresh dictionary page
# per column per row group, and with 60,660 distinct `gene_id`s a 100k-row
# group amortises a ~1 MB dictionary over only ~1.6 aliquots — so the
# dictionary gets rewritten 747 times across TCGA-BRCA's expression table
# and `gene_id` alone costs 497 MB, 6.6 bytes for every row. Raising the
# group to 1M rows (~16 aliquots) cut the whole table from 1.88 to 1.13 GiB
# with no schema change at all.
_ROW_GROUP_SIZE_GENE_SCALE = 1_000_000
# zstd over snappy is a further ~25% on these tables and costs only write
# time; every reader in the ecosystem handles it.
_COMPRESSION = "zstd"
_ROW_GROUP_SIZE_WIDE_ROW = 50
_ROW_GROUP_SIZE_BY_TABLE: dict[str, int] = {
    # One row is the patient's full nested GDC case tree.
    "cases": _ROW_GROUP_SIZE_WIDE_ROW,
    # One row embeds a scanned PDF — hundreds of KB of bytes.
    "pathology_report": _ROW_GROUP_SIZE_WIDE_ROW,
    "gene_expression_quantification": _ROW_GROUP_SIZE_GENE_SCALE,
    "gene_level_copy_number": _ROW_GROUP_SIZE_GENE_SCALE,
    # ~486k probes per 450k array file — same dictionary-page problem as the
    # per-gene tables, same fix.
    "methylation_beta_value": _ROW_GROUP_SIZE_GENE_SCALE,
    # The gene model is one small table read whole; a single row group is
    # both the smallest and the fastest thing to hand a reader.
    "gene_model": 100_000,
    # ssGSEA scores are the same shape as expression -- narrow rows, very
    # many of them (575k per project for Hallmark, ~17M for a Reactome-sized
    # collection) -- so they get the same row-group treatment.
    **{f"ssgsea_scores_{c}": _ROW_GROUP_SIZE_GENE_SCALE for c in SSGSEA_COLLECTIONS},
}


def write_tables(
    tables: dict[str, Any],
    processed_dir: Path,
    project_id: str,
    counts: dict[str, int] | None = None,
) -> dict[str, Path]:
    """Write each table to `<processed_dir>/<project_id>/<table>/data.parquet`.

    Tables in `TABULAR_TABLES` use their fixed pan-cancer schema. Tables
    not in that map (today: `clinical_supplement_*`) get schema inferred
    from rows — each project ships only the columns its biotab forms
    contain.

    A table's value is either a list of row dicts or an iterator of
    `pa.RecordBatch`. The per-gene tables use the latter — they are far too
    large to materialise as dicts — and are streamed straight to Parquet.
    Because a streamed table's row count isn't knowable until it has been
    written, `counts` is filled in with the rows written per table so
    callers can report totals without holding the rows.

    Uses the same row-group + page-index settings as the consolidated build
    so HF Data Studio can scan large tables (notably `expression`) without
    hitting its scan-size limit.
    """
    import pyarrow as pa
    import pyarrow.parquet as pq

    out_paths: dict[str, Path] = {}
    project_dir = processed_dir / project_id
    project_dir.mkdir(parents=True, exist_ok=True)
    for table_name, rows in tables.items():
        out_path = project_dir / table_name / "data.parquet"
        row_group_size = _ROW_GROUP_SIZE_BY_TABLE.get(table_name, _ROW_GROUP_SIZE_DEFAULT)
        if not isinstance(rows, list):
            # A stream of `pa.RecordBatch` — used by tables too large to
            # materialise as Python dicts (see
            # `_gene_level_copy_number_batches`). Written incrementally so
            # peak memory stays at one batch.
            n_written = _write_batches(
                rows, out_path, TABULAR_TABLES[table_name], row_group_size
            )
            if n_written is not None:
                out_paths[table_name] = out_path
                if counts is not None:
                    counts[table_name] = n_written
            continue
        out_path.parent.mkdir(parents=True, exist_ok=True)
        if table_name in TABULAR_TABLES:
            table = pa.Table.from_pylist(rows, schema=TABULAR_TABLES[table_name])
        elif rows:
            # Flex-schema table — let pyarrow infer column types per project.
            # Empty-table case skipped: zero-row biotab forms (e.g. LAML's
            # missing follow_up) write nothing rather than a 0×0 parquet.
            table = pa.Table.from_pylist(rows)
        else:
            continue
        pq.write_table(
            table,
            out_path,
            compression=_COMPRESSION,
            row_group_size=row_group_size,
            write_page_index=True,
        )
        out_paths[table_name] = out_path
        if counts is not None:
            counts[table_name] = len(rows)
    return out_paths


def _write_batches(
    batches: Any,
    out_path: Path,
    schema: Any,
    row_group_size: int,
) -> int | None:
    """Stream record batches to one Parquet file; None if there were none.

    Batches arrive one per source file (~60k rows), which is far too small
    to be a row group of its own — see `_ROW_GROUP_SIZE_BY_TABLE` for why
    small groups are expensive when the gene cardinality is 60,660. They are
    accumulated until they reach `row_group_size` and flushed together, so
    the on-disk layout matches what the non-streaming path produces.

    Nothing is created until the first batch arrives, so a project with no
    files for this modality leaves no empty parquet behind.
    """
    import pyarrow as pa
    import pyarrow.parquet as pq

    writer: pq.ParquetWriter | None = None
    pending: list[pa.RecordBatch] = []
    pending_rows = 0
    total_rows = 0

    def flush() -> None:
        nonlocal pending, pending_rows
        if not pending:
            return
        assert writer is not None
        # Explicit row_group_size rather than pyarrow's default, so the
        # on-disk grouping doesn't drift if that default ever changes.
        writer.write_table(
            pa.Table.from_batches(pending, schema=schema), row_group_size=row_group_size
        )
        pending = []
        pending_rows = 0

    try:
        for batch in batches:
            if writer is None:
                out_path.parent.mkdir(parents=True, exist_ok=True)
                writer = pq.ParquetWriter(
                    out_path,
                    schema,
                    compression=_COMPRESSION,
                    write_page_index=True,
                )
            pending.append(batch)
            pending_rows += batch.num_rows
            total_rows += batch.num_rows
            if pending_rows >= row_group_size:
                flush()
        flush()
    finally:
        if writer is not None:
            writer.close()
    return total_rows if writer is not None else None
