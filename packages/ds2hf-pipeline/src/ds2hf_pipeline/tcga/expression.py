"""Parse open-access GDC Gene Expression Quantification TSVs (RNA-Seq STAR counts).

Each TSV is one RNA-Seq aliquot's STAR-aligned gene quantifications: ~60,660
GENCODE v36 gene rows × 9 columns, plus 4 N_* summary rows at the top. We:

  - Lift the 4 N_* QC rows out as scalar fields (using their `unstranded` value).
  - Keep all 4 unstranded value columns (raw counts + TPM + FPKM + FPKM-UQ).
  - Drop `stranded_first` and `stranded_second`. The GDC pipeline (see
    reference below) harmonizes by treating all reads as unstranded, so
    `unstranded` is the canonical column for cross-sample analysis.
  - Preserve full gene metadata (gene_id with GENCODE version suffix, gene_name,
    gene_type) so the dataset is self-describing without a side vocab file.

GDC pipeline reference:
  https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/
"""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

_QC_NAMES = ["N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous"]


def _file_aliquot_and_case(entry: dict[str, Any]) -> tuple[str | None, str | None]:
    """Return (case_id, aliquot_id) for one manifest entry.

    Expression files are per-aliquot, so we expect exactly one case and one
    aliquot. Returns (None, None) if the manifest's nested cases data is
    ambiguous — caller should skip the file.
    """
    case_ids: set[str] = set()
    aliquot_ids: set[str] = set()
    for case in entry.get("cases", []):
        case_ids.add(case["case_id"])
        for sample in case.get("samples", []):
            for portion in sample.get("portions") or []:
                for analyte in portion.get("analytes") or []:
                    for aliquot in analyte.get("aliquots") or []:
                        aliquot_ids.add(aliquot["aliquot_id"])
    if len(case_ids) != 1 or len(aliquot_ids) != 1:
        return None, None
    return next(iter(case_ids)), next(iter(aliquot_ids))


def _parse_expression_tsv(path: Path) -> tuple[dict[str, int | None], pd.DataFrame]:
    """Read the TSV; return ({QC_name: int_or_None}, gene_rows_df).

    The TSV starts with one `# gene-model: GENCODE v36` line (skipped via
    `comment="#"`), then a header line, then 4 N_* summary rows, then the
    ENSG gene rows. N_* rows have empty `gene_name`/`gene_type` cells.
    """
    df = pd.read_csv(path, sep="\t", comment="#", low_memory=False)
    qc_mask = df["gene_id"].str.startswith("N_", na=False)
    qc_rows = df[qc_mask].set_index("gene_id")
    gene_rows = df[~qc_mask].reset_index(drop=True)

    qc: dict[str, int | None] = {}
    for name in _QC_NAMES:
        if name in qc_rows.index and pd.notna(qc_rows.loc[name, "unstranded"]):
            qc[name] = int(qc_rows.loc[name, "unstranded"])
        else:
            qc[name] = None
    return qc, gene_rows


def _to_int_list(series: pd.Series) -> list[int | None]:
    """Pandas float-with-NaN -> Python list of int-or-None for pyarrow."""
    out: list[int | None] = []
    for v in series:
        if pd.isna(v):
            out.append(None)
        else:
            out.append(int(v))
    return out


def _to_float_list(series: pd.Series) -> list[float | None]:
    out: list[float | None] = []
    for v in series:
        if pd.isna(v):
            out.append(None)
        else:
            out.append(float(v))
    return out


def _build_record(
    qc: dict[str, int | None],
    gene_rows: pd.DataFrame,
    *,
    aliquot_id: str,
    source_file_id: str,
) -> dict[str, Any]:
    return {
        "sample_id": None,  # resolved at attach time from the patient's samples
        "aliquot_id": aliquot_id,
        "source_file_id": source_file_id,
        "N_unmapped": qc["N_unmapped"],
        "N_multimapping": qc["N_multimapping"],
        "N_noFeature": qc["N_noFeature"],
        "N_ambiguous": qc["N_ambiguous"],
        "gene_id": gene_rows["gene_id"].tolist(),
        "gene_name": gene_rows["gene_name"].tolist(),
        "gene_type": gene_rows["gene_type"].tolist(),
        "unstranded": _to_int_list(gene_rows["unstranded"]),
        "tpm_unstranded": _to_float_list(gene_rows["tpm_unstranded"]),
        "fpkm_unstranded": _to_float_list(gene_rows["fpkm_unstranded"]),
        "fpkm_uq_unstranded": _to_float_list(gene_rows["fpkm_uq_unstranded"]),
    }


def load_for_project(project_raw_dir: Path) -> dict[str, list[dict[str, Any]]]:
    """Return {case_id: [expression_record, ...]} for raw/<PROJECT>/expression/.

    Empty dict if the directory or its manifest is missing — letting `build`
    proceed without expression data when it hasn't been fetched.
    """
    expr_dir = project_raw_dir / "expression"
    manifest_path = expr_dir / "manifest.json"
    if not manifest_path.exists():
        return {}

    manifest = json.loads(manifest_path.read_text())
    by_case: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for entry in manifest:
        file_path = expr_dir / entry["file_name"]
        if not file_path.exists():
            continue
        case_id, aliquot_id = _file_aliquot_and_case(entry)
        if not case_id or not aliquot_id:
            continue
        qc, gene_rows = _parse_expression_tsv(file_path)
        record = _build_record(
            qc, gene_rows, aliquot_id=aliquot_id, source_file_id=entry["file_id"]
        )
        by_case[case_id].append(record)
    return dict(by_case)


# GENCODE gives immunoglobulin and T-cell-receptor segments their own
# biotypes rather than `protein_coding`, so a plain protein-coding filter
# silently drops all of them. That matters for tumour immunology: IGKC is
# among the ~400 highest-expressed genes in an infiltrated tumour and
# TRAC/TRBC2 are top-2500, and MSigDB's B-cell-receptor, Fc-receptor and
# complement pathways are largely built from these segments — without them
# REACTOME_CD22_MEDIATED_BCR_REGULATION matches 5 of its 61 genes.
#
# We include the functional segments and exclude the ~237 pseudogene
# biotypes, which are near-uniformly zero and would only pad the low ranks.
#
# Caveat worth carrying: V/D/J segments are somatically rearranged and
# clonally expanded, so their bulk expression reports lymphocyte
# infiltration and clonality rather than regulation of a fixed locus. That
# is what makes them informative about the tumour immune microenvironment,
# but it is a different kind of measurement from the rest of the matrix.
IMMUNE_RECEPTOR_GENE_TYPES: tuple[str, ...] = (
    "IG_V_gene",
    "IG_C_gene",
    "IG_D_gene",
    "IG_J_gene",
    "TR_V_gene",
    "TR_C_gene",
    "TR_D_gene",
    "TR_J_gene",
)


def tpm_matrix_for_project(
    project_raw_dir: Path,
) -> tuple[list[str], np.ndarray, list[dict[str, Any]]]:
    """Return `(gene_symbols, tpm_matrix, aliquot_records)` for ssGSEA scoring.

    The matrix is genes x aliquots of `tpm_unstranded`, restricted to
    **protein-coding genes plus the functional immunoglobulin / T-cell
    receptor segments** (see `IMMUNE_RECEPTOR_GENE_TYPES`), with the
    `_PAR_Y` duplicates dropped.

    These are load-bearing choices rather than tidying: ssGSEA ranks genes
    *within the supplied matrix*, so the gene universe changes every score
    and cannot be revised after publication without invalidating them. The
    ~40k lncRNA / pseudogene rows are mostly zero and would dominate the low
    ranks; the `_PAR_Y` entries are duplicate annotations of the X copy that
    GDC's STAR pipeline leaves at exactly 0.0, so dropping them is lossless
    and resolves all but a handful of duplicate symbols. Those few remaining
    duplicates are collapsed by max TPM.

    Adding the 409 immune-receptor segments moves existing scores by a
    median of 0.21 of a pathway's cross-sample SD (p95 0.54, per-pathway
    Spearman >= 0.986 across samples) and takes the affected MSigDB immune
    pathways from 8-28% gene coverage to 99-100%.

    All GDC STAR files share one gene model (GENCODE v36) regardless of the
    release they first appeared in, so a single gene index is valid across
    the whole cohort; the ordering here comes from the first file read and
    every other file is reindexed onto it.

    Files whose manifest entry doesn't resolve to exactly one case and one
    aliquot are skipped, matching `load_for_project`.
    """
    expr_dir = project_raw_dir / "expression"
    manifest_path = expr_dir / "manifest.json"
    if not manifest_path.exists():
        return [], np.empty((0, 0), dtype=np.float32), []

    entries = []
    for entry in json.loads(manifest_path.read_text()):
        file_path = expr_dir / entry["file_name"]
        if not file_path.exists():
            continue
        case_id, aliquot_id = _file_aliquot_and_case(entry)
        if not case_id or not aliquot_id:
            continue
        entries.append((file_path, case_id, aliquot_id, entry["file_id"]))
    if not entries:
        return [], np.empty((0, 0), dtype=np.float32), []

    cols = ["gene_id", "gene_name", "gene_type", "tpm_unstranded"]
    gene_index: pd.Index | None = None
    matrix: np.ndarray | None = None
    records: list[dict[str, Any]] = []
    for i, (file_path, case_id, aliquot_id, source_file_id) in enumerate(entries):
        df = pd.read_csv(file_path, sep="\t", comment="#", low_memory=False, usecols=cols)
        in_universe = df["gene_type"].eq("protein_coding") | df["gene_type"].isin(
            IMMUNE_RECEPTOR_GENE_TYPES
        )
        df = df[in_universe & (~df["gene_id"].str.contains("PAR_Y"))]
        series = df.groupby("gene_name")["tpm_unstranded"].max()
        if gene_index is None:
            gene_index = series.index
            matrix = np.empty((len(gene_index), len(entries)), dtype=np.float32)
        matrix[:, i] = series.reindex(gene_index).to_numpy(dtype=np.float32)
        records.append(
            {"case_id": case_id, "aliquot_id": aliquot_id, "source_file_id": source_file_id}
        )
    assert gene_index is not None and matrix is not None
    return list(gene_index), matrix, records


def attach(
    rows: list[dict[str, Any]], by_case: dict[str, list[dict[str, Any]]]
) -> list[dict[str, Any]]:
    """Mutate `rows` to populate samples_gene_expression_quantification.

    Resolves each record's sample_id by looking up its aliquot_id in the
    patient row's `samples` field (source of truth for every aliquot the case
    has). Sorts records by aliquot_id for deterministic output.
    """
    for row in rows:
        records = by_case.get(row["case_id"], [])
        # Walk the full GDC tree: samples → portions → analytes → aliquots
        aliquot_to_sample: dict[str, str] = {}
        for s in row.get("samples") or []:
            sid = s.get("sample_id")
            if not sid:
                continue
            for portion in s.get("portions") or []:
                for analyte in portion.get("analytes") or []:
                    for a in analyte.get("aliquots") or []:
                        aq = a.get("aliquot_id")
                        if aq:
                            aliquot_to_sample[aq] = sid
        for r in records:
            r["sample_id"] = aliquot_to_sample.get(r["aliquot_id"])
        records.sort(key=lambda r: r.get("aliquot_id") or "")
        row["samples_gene_expression_quantification"] = records
    return rows
