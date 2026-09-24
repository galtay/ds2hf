"""single-sample GSEA (ssGSEA), transcribed from Bioconductor GSVA.

Pathway activity per RNA-Seq aliquot: expression in, one score per gene set
out. Barbie et al. 2009's method as implemented by GSVA's `ssgsea.R`.

## Why a Python transcription rather than calling GSVA

The canonical implementation is R. Rather than make Bioconductor a pipeline
dependency, this module reproduces it and is validated against it: on a
100-sample × 19,938-gene TCGA matrix scored with MSigDB Hallmark, agreement
with GSVA 2.6.6 is Pearson/Spearman 1.0000000000, max *relative* difference
4.8e-13 (raw) and 6.0e-13 (normalized), with identical pathway ordering in
100/100 samples. Residuals are float64 summation-order noise. The R side is
kept as an offline oracle under `dev_research/ssgsea/`, in the same spirit
as the Liu 2018 CDR validation.

## The algorithm

For each sample independently:

    R    = integer rank of expression within the sample. GSVA computes
           `as.integer(rank(x))` — average ties, then truncate toward zero.
    ord  = genes ordered by decreasing rank (highest expressed first)
    in   = cumsum(R[ord]^alpha * inSet) / sum(R[ord]^alpha * inSet)
    out  = cumsum(!inSet) / sum(!inSet)
    ES   = sum(in - out)

`alpha` weights the walk by **rank**, not by expression. A consequence
worth knowing: any strictly monotonic transform of expression leaves ranks —
and therefore scores — bit-identical, so ssGSEA(TPM) == ssGSEA(log1p(TPM)).
There is no reason to log-transform before scoring.

## Normalization is deliberately separate

GSVA's `normalize=TRUE` divides every score by the range of the *entire*
score matrix. That is a single global scalar, so it is exactly decomposable
across chunks — but it also makes each score depend on the cohort and the
gene-set collection it was computed alongside, not on the sample alone.
Adding MSigDB Reactome to a Hallmark run widens the divisor by ~49% on TCGA
data, silently rescaling every previously published Hallmark value.

So `ssgsea_scores` returns raw scores, and `normalize_global` returns the
divisor alongside the normalized values. Publish the raw score and the
divisor; the normalized view is then reconstructable, while the reverse is
not.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

DEFAULT_ALPHA = 0.25

# A floor is defensible on stability: leave-one-gene-out sensitivity falls
# smoothly with set size (a dropped gene moves the score by ~98% of its
# cross-sample SD at n=5, 62% at n=10, 17% at n=100), so very small sets are
# noise-dominated. There is no natural breakpoint, and 10 is a judgement
# call rather than a standard.
DEFAULT_MIN_SIZE = 10

# No ceiling. `maxSize=500` is inherited convention (the Broad's GSEA
# desktop uses 15/500, gseapy 15/2000) and is *not* a GSVA default — GSVA
# ships minSize=1, maxSize=Inf. The rationale behind it belongs to
# permutation-based GSEA, where huge sets are non-specific hypotheses and
# inflate multiple testing; ssGSEA runs no permutation test and we publish
# descriptive scores. Empirically large sets are the *most* stable, and set
# size does not drive score magnitude (Spearman -0.001 against the mean
# score across 50 Hallmark sets). Dropping a set at build time is
# unrecoverable, whereas consumers can filter on the per-pathway gene counts
# we ship on every row.
DEFAULT_MAX_SIZE: int | None = None


def average_rank(x: np.ndarray) -> np.ndarray:
    """Ranks with ties averaged — equivalent to R's `rank()` default.

    Vectorised so we don't need scipy for `rankdata`; verified identical to
    `scipy.stats.rankdata(method="average")` on real TCGA expression
    columns (which are ~25% ties, since unexpressed genes are all 0.0).
    """
    n = x.size
    order = np.argsort(x, kind="stable")
    xs = x[order]
    is_new = np.empty(n, dtype=bool)
    is_new[0] = True
    np.not_equal(xs[1:], xs[:-1], out=is_new[1:])
    group = np.cumsum(is_new) - 1
    starts = np.flatnonzero(is_new)
    ends = np.append(starts[1:], n) - 1
    # Mean of the 1-based ranks spanned by each tie group.
    avg = (starts + ends + 2) / 2.0
    out = np.empty(n, dtype=np.float64)
    out[order] = avg[group]
    return out


def _sample_ranks(x: np.ndarray) -> np.ndarray:
    """GSVA's `as.integer(rank(x))`: average ties, then truncate toward zero."""
    return average_rank(x).astype(np.int64)


def ssgsea_scores(
    expression: np.ndarray,
    gene_sets: list[np.ndarray],
    alpha: float = DEFAULT_ALPHA,
) -> np.ndarray:
    """Raw (un-normalized) ssGSEA scores for a genes x samples matrix.

    `gene_sets` are arrays of row indices into `expression` — use
    `map_gene_sets` to build them from symbols. Returns a
    `len(gene_sets)` x `n_samples` array.

    Scores are computed one sample at a time and depend on no other sample,
    so callers may chunk the sample axis freely: concatenating chunked
    results is bit-identical to a single call. Only `normalize_global`
    couples samples together.
    """
    if expression.ndim != 2:
        raise ValueError(f"expression must be 2-D (genes x samples), got {expression.shape}")
    n_genes, n_samples = expression.shape
    for i, gs in enumerate(gene_sets):
        if gs.size and (gs.min() < 0 or gs.max() >= n_genes):
            raise ValueError(f"gene_sets[{i}] has indices outside the expression matrix")

    # Closed form of sum(walkStat), which is what GSVA itself evaluates.
    # Writing p for a gene's 1-based position in the descending-rank order,
    # w for its weight, and summing the two step CDFs over all k:
    #
    #   ES = sum_{p in set} w_p * (n - p + 1) / sum_{p in set} w_p
    #        - (n(n+1)/2 - sum_{p in set} (n - p + 1)) / (n - |set|)
    #
    # The out-of-set half telescopes into the total, so every term depends
    # only on the set's own members. That makes each gene set O(|set|)
    # rather than O(n_genes) -- the difference between minutes and hours
    # once a 1,450-set collection like Reactome is in play.
    total_positions = n_genes * (n_genes + 1) / 2.0
    out = np.empty((len(gene_sets), n_samples), dtype=np.float64)
    for j in range(n_samples):
        ranks = _sample_ranks(expression[:, j])
        # Descending rank; stable so tied genes keep matrix order, matching
        # R's `order()`, which is stable for the integer ranks GSVA passes it.
        order = np.argsort(-ranks, kind="stable")
        position = np.empty(n_genes, dtype=np.float64)
        position[order] = np.arange(1, n_genes + 1)
        weight = np.abs(ranks).astype(np.float64) ** alpha
        for i, gs in enumerate(gene_sets):
            tail = n_genes - position[gs] + 1.0  # (n - p + 1) per member
            w = weight[gs]
            out[i, j] = (w * tail).sum() / w.sum() - (total_positions - tail.sum()) / (
                n_genes - gs.size
            )
    return out


def normalize_global(raw: np.ndarray) -> tuple[np.ndarray, float]:
    """GSVA's `normalize=TRUE`: divide by the range of the whole matrix.

    Returns `(normalized, divisor)`. The divisor is returned rather than
    swallowed because it is a property of the run — the cohort and the
    gene-set collection that were scored together — not of the samples. A
    published dataset should carry the raw scores and this number, so the
    normalized view stays reconstructable and a later release that adds
    gene sets doesn't silently restate earlier values.
    """
    divisor = float(raw.max() - raw.min())
    if divisor == 0.0:
        raise ValueError("degenerate score matrix: max == min, nothing to normalize by")
    return raw / divisor, divisor


# ---------------------------------------------------------------------------
# Gene sets (GMT)
# ---------------------------------------------------------------------------


def load_gmt(path: Path | str) -> dict[str, list[str]]:
    """Parse a GMT file into `{set_name: [gene_symbol, ...]}`.

    GMT is tab-separated: name, description/URL, then the members.
    """
    sets: dict[str, list[str]] = {}
    with Path(path).open(encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) > 2 and parts[0]:
                sets[parts[0]] = [g for g in parts[2:] if g]
    return sets


def map_gene_sets(
    sets: dict[str, list[str]],
    genes: list[str] | np.ndarray,
    min_size: int = DEFAULT_MIN_SIZE,
    max_size: int | None = DEFAULT_MAX_SIZE,
) -> tuple[dict[str, np.ndarray], list[dict[str, object]]]:
    """Intersect each gene set with `genes`, then filter by post-mapping size.

    `max_size=None` means no upper bound. Returns `(kept, stats)`. GSVA
    applies the size filter *after* mapping onto the expression matrix, so
    the sizes that matter are the matched ones — that is what `stats`
    reports, per set, for provenance: `original_gene_count`,
    `matched_gene_count`, `match_fraction`, `kept`.

    `match_fraction` is the more informative filter of the two for
    consumers. A set can be small because the pathway genuinely has few
    genes (Reactome's leaf nodes are often 5-9, with matched == declared)
    or because most of its genes fall outside the expression matrix's gene
    universe — two different problems that a size floor alone conflates.
    """
    index = {g: i for i, g in enumerate(genes)}
    kept: dict[str, np.ndarray] = {}
    stats: list[dict[str, object]] = []
    for name, members in sets.items():
        unique = set(members)
        idx = np.array(sorted({index[m] for m in unique if m in index}), dtype=np.int64)
        keep = idx.size >= min_size and (max_size is None or idx.size <= max_size)
        stats.append(
            {
                "pathway": name,
                "original_gene_count": len(unique),
                "matched_gene_count": int(idx.size),
                "match_fraction": (idx.size / len(unique)) if unique else 0.0,
                "kept": keep,
            }
        )
        if keep:
            kept[name] = idx
    return kept, stats


# ---------------------------------------------------------------------------
# Consolidated (patient) layout
#
# The tabular layout emits one row per (aliquot, pathway); the patient layout
# wants one struct-of-arrays record per aliquot, matching how expression is
# stored there. Both derive from the same `ssgsea_scores` call, so the two
# published views carry identical numbers by construction.
# ---------------------------------------------------------------------------


def score_project(
    project_raw_dir,
    collection: str,
    msigdb_dir,
    load_matrix=None,
) -> tuple[list[str], np.ndarray, list[dict[str, object]], dict[str, dict[str, object]]]:
    """Score one project against one collection.

    Returns `(pathway_names, scores, aliquot_records, geneset_stats)` where
    `scores` is pathways x aliquots. Imports the expression reader lazily so
    this module stays importable with numpy alone.
    """
    from ds2hf_pipeline.tcga import expression as _expression
    from ds2hf_pipeline.tcga import msigdb as _msigdb

    gmt_path = msigdb_dir / _msigdb.COLLECTIONS[collection].file_name
    if not gmt_path.exists():
        return [], np.empty((0, 0)), [], {}
    loader = load_matrix or (lambda: _expression.tpm_matrix_for_project(project_raw_dir))
    genes, matrix, records = loader()
    if not records:
        return [], np.empty((0, 0)), [], {}
    kept, stats = map_gene_sets(load_gmt(gmt_path), genes)
    if not kept:
        return [], np.empty((0, 0)), [], {}
    names = list(kept)
    scores = ssgsea_scores(matrix.astype(np.float64), [kept[n] for n in names])
    return names, scores, records, {s["pathway"]: s for s in stats}


def load_for_project(
    project_raw_dir,
    collection: str,
    msigdb_dir,
    load_matrix=None,
) -> dict[str, list[dict[str, object]]]:
    """Return `{case_id: [per-aliquot struct-of-arrays record, ...]}`.

    Shaped for the consolidated patient row: one record per scored aliquot,
    with `pathway` / `pathway_url` / gene counts / `score_raw` as
    index-aligned arrays.
    """
    from ds2hf_pipeline.tcga import msigdb as _msigdb

    names, scores, records, stats = score_project(
        project_raw_dir, collection, msigdb_dir, load_matrix
    )
    if not names:
        return {}
    urls = [_msigdb.geneset_url(n) for n in names]
    matched = [int(stats[n]["matched_gene_count"]) for n in names]
    original = [int(stats[n]["original_gene_count"]) for n in names]

    by_case: dict[str, list[dict[str, object]]] = {}
    for j, rec in enumerate(records):
        by_case.setdefault(rec["case_id"], []).append(
            {
                "sample_id": None,  # resolved in `attach` from the patient's samples
                "aliquot_id": rec["aliquot_id"],
                "source_file_id": rec["source_file_id"],
                "pathway": names,
                "pathway_url": urls,
                "matched_gene_count": matched,
                "original_gene_count": original,
                "score_raw": [float(v) for v in scores[:, j]],
            }
        )
    return by_case


def attach(
    rows: list[dict[str, object]],
    by_case: dict[str, list[dict[str, object]]],
    column: str,
) -> list[dict[str, object]]:
    """Populate `column` on each patient row, resolving each record's sample_id.

    Mirrors `expression.attach`: the patient row's `samples` tree is the
    source of truth for which sample an aliquot belongs to.
    """
    for row in rows:
        records = by_case.get(row["case_id"], [])
        aliquot_to_sample: dict[str, str] = {}
        for s in row.get("samples") or []:
            sid = s.get("sample_id")
            if not sid:
                continue
            for portion in s.get("portions") or []:
                for analyte in portion.get("analytes") or []:
                    for a in analyte.get("aliquots") or []:
                        if a.get("aliquot_id"):
                            aliquot_to_sample[a["aliquot_id"]] = sid
        for r in records:
            r["sample_id"] = aliquot_to_sample.get(r["aliquot_id"])
        records.sort(key=lambda r: r.get("aliquot_id") or "")
        row[column] = records
    return rows
