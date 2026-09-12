"""Tests for the ssGSEA transcription.

The canonical implementation is Bioconductor GSVA. Rather than depend on R,
we validate against it: `dev_research/ssgsea/` holds a small deterministic
fixture and the scores GSVA 2.6.6 produced for it, and
`test_matches_gsva_reference` asserts we reproduce them. The remaining
tests pin properties of the algorithm that the fixture alone wouldn't
catch — the rank-invariance, the chunk-decomposability, and the size
filtering GSVA applies after mapping.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from tcga2hf_pipeline import ssgsea

FIXTURES = Path(__file__).resolve().parents[3] / "dev_research" / "ssgsea"
EXPR_TSV = FIXTURES / "fixture_expr.tsv"
SETS_GMT = FIXTURES / "fixture_sets.gmt"
EXPECTED_CSV = FIXTURES / "fixture_expected_raw.csv"


def _load_fixture() -> tuple[np.ndarray, list[str], dict[str, np.ndarray]]:
    import csv

    with EXPR_TSV.open() as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    genes = [r[0] for r in rows[1:]]
    expr = np.array([[float(v) for v in r[1:]] for r in rows[1:]], dtype=np.float64)
    sets = ssgsea.load_gmt(SETS_GMT)
    # Pinned to the filter the GSVA reference run used, so the comparison
    # keeps testing the same thing if the pipeline defaults change.
    kept, _ = ssgsea.map_gene_sets(sets, genes, min_size=10, max_size=500)
    return expr, genes, kept


# ---------------------------------------------------------------------------
# Agreement with the GSVA reference
# ---------------------------------------------------------------------------


def test_matches_gsva_reference() -> None:
    """Reproduce GSVA 2.6.6's raw ssGSEA scores on the committed fixture.

    Tolerance is relative and tight: the two implementations differ only in
    float64 summation order, so anything above ~1e-9 relative means the
    algorithm itself has drifted.
    """
    import csv

    expr, _, kept = _load_fixture()
    with EXPECTED_CSV.open() as fh:
        rows = list(csv.reader(fh))
    header = rows[0]
    assert header[0] in ("", '""'), "expected an R-style row-name column"
    expected_names = [r[0] for r in rows[1:]]
    expected = np.array([[float(v) for v in r[1:]] for r in rows[1:]], dtype=np.float64)

    names = list(kept)
    assert sorted(names) == sorted(expected_names), (
        "gene sets kept by map_gene_sets differ from those GSVA kept"
    )
    order = [names.index(n) for n in expected_names]
    got = ssgsea.ssgsea_scores(expr, [kept[n] for n in names])[order]

    rel = np.abs(got - expected) / np.maximum(np.abs(expected), 1e-12)
    assert rel.max() < 1e-9, f"max relative deviation from GSVA: {rel.max():.3e}"


def test_reference_fixture_excludes_undersized_set() -> None:
    """GSVA filters on post-mapping size; the 5-gene set must not appear."""
    _, genes, kept = _load_fixture()
    assert "SET_TOO_SMALL_5" not in kept
    assert "SET_A_12" in kept


# ---------------------------------------------------------------------------
# Algorithmic properties
# ---------------------------------------------------------------------------


def test_monotonic_transform_leaves_scores_identical() -> None:
    """alpha weights ranks, not expression, so log1p is a no-op.

    This is why the pipeline scores TPM directly instead of log-transforming:
    the transform provably cannot change the answer.
    """
    expr, _, kept = _load_fixture()
    sets = list(kept.values())
    a = ssgsea.ssgsea_scores(expr, sets)
    b = ssgsea.ssgsea_scores(np.log1p(expr), sets)
    assert np.array_equal(a, b)
    # Any strictly increasing map, not just log1p.
    c = ssgsea.ssgsea_scores(np.sqrt(expr) * 3.0 + 1.0, sets)
    assert np.array_equal(a, c)


def test_samples_are_independent_so_chunking_is_exact() -> None:
    """Concatenated chunks must be bit-identical to one call.

    The pipeline relies on this to score the cohort in pieces.
    """
    expr, _, kept = _load_fixture()
    sets = list(kept.values())
    whole = ssgsea.ssgsea_scores(expr, sets)
    chunked = np.hstack(
        [ssgsea.ssgsea_scores(expr[:, s], sets) for s in (slice(0, 3), slice(3, 5), slice(5, None))]
    )
    assert np.array_equal(whole, chunked)


def test_normalization_is_a_single_global_scalar() -> None:
    """Normalizing is pure scaling: it preserves sign, zero and all ratios."""
    expr, _, kept = _load_fixture()
    raw = ssgsea.ssgsea_scores(expr, list(kept.values()))
    norm, divisor = ssgsea.normalize_global(raw)

    assert divisor == pytest.approx(raw.max() - raw.min())
    assert np.allclose(norm * divisor, raw)
    assert np.array_equal(np.sign(norm), np.sign(raw))
    # A divisor accumulated over chunks reproduces the single-call result,
    # which is what lets the cohort be scored in pieces and normalized once.
    lo = min(raw[:, :4].min(), raw[:, 4:].min())
    hi = max(raw[:, :4].max(), raw[:, 4:].max())
    assert np.array_equal(raw / (hi - lo), norm)


def test_normalize_rejects_degenerate_matrix() -> None:
    with pytest.raises(ValueError, match="degenerate"):
        ssgsea.normalize_global(np.full((3, 4), 2.5))


# ---------------------------------------------------------------------------
# Ranking
# ---------------------------------------------------------------------------


def test_average_rank_matches_r_semantics() -> None:
    """Ties average, as R's rank() does — TPM data is ~30% ties (zeros)."""
    assert np.array_equal(ssgsea.average_rank(np.array([3.0, 1.0, 2.0])), [3.0, 1.0, 2.0])
    assert np.array_equal(ssgsea.average_rank(np.array([1.0, 1.0, 1.0])), [2.0, 2.0, 2.0])
    assert np.array_equal(
        ssgsea.average_rank(np.array([2.0, 2.0, 1.0, 1.0, 3.0])), [3.5, 3.5, 1.5, 1.5, 5.0]
    )
    # All-zero genes must share one averaged rank rather than an arbitrary order.
    r = ssgsea.average_rank(np.array([0.0, 0.0, 0.0, 0.0, 5.0]))
    assert np.array_equal(r, [2.5, 2.5, 2.5, 2.5, 5.0])


# ---------------------------------------------------------------------------
# Gene-set mapping
# ---------------------------------------------------------------------------


def test_map_gene_sets_filters_on_matched_size_not_declared_size(tmp_path: Path) -> None:
    """A large set that barely intersects the matrix must be dropped.

    GSVA applies minSize/maxSize after mapping onto the expression genes, so
    a 200-gene pathway with 3 genes present is a 3-gene set, not a 200-gene
    one. Getting this backwards would silently score near-empty pathways.
    """
    gmt = tmp_path / "s.gmt"
    gmt.write_text(
        "BIG_BUT_ABSENT\tdesc\t" + "\t".join(f"MISSING{i}" for i in range(200)) + "\tA\tB\tC\n"
        "PRESENT\tdesc\t" + "\t".join(f"GENE{i}" for i in range(12)) + "\n"
    )
    genes = [f"GENE{i}" for i in range(50)] + ["A", "B", "C"]
    kept, stats = ssgsea.map_gene_sets(ssgsea.load_gmt(gmt), genes, min_size=10, max_size=500)

    assert "BIG_BUT_ABSENT" not in kept
    assert "PRESENT" in kept
    by_name = {s["pathway"]: s for s in stats}
    assert by_name["BIG_BUT_ABSENT"]["original_gene_count"] == 203
    assert by_name["BIG_BUT_ABSENT"]["matched_gene_count"] == 3
    assert by_name["BIG_BUT_ABSENT"]["kept"] is False
    assert by_name["PRESENT"]["match_fraction"] == 1.0


def test_map_gene_sets_deduplicates_members(tmp_path: Path) -> None:
    """Repeated symbols in a GMT line count once, and indices are unique."""
    gmt = tmp_path / "s.gmt"
    gmt.write_text("S\tdesc\t" + "\t".join(["A"] * 5 + ["B", "C"]) + "\n")
    kept, stats = ssgsea.map_gene_sets(ssgsea.load_gmt(gmt), ["A", "B", "C", "D"], min_size=1)
    assert kept["S"].size == 3
    assert stats[0]["original_gene_count"] == 3


def test_ssgsea_rejects_out_of_range_indices() -> None:
    """A mis-built gene set should fail loudly, not index into the wrong gene."""
    expr = np.arange(20, dtype=np.float64).reshape(10, 2)
    with pytest.raises(ValueError, match="outside the expression matrix"):
        ssgsea.ssgsea_scores(expr, [np.array([0, 99])])


# ---------------------------------------------------------------------------
# Tabular emitters
# ---------------------------------------------------------------------------


def test_stats_rows_carry_project_and_pan_cancer(tmp_path: Path) -> None:
    """Each project's stats must include the pan-cancer reference too.

    That duplication is the point: a consumer who loads one project's config
    can still z-score against all of TCGA without scanning 33 of them.
    """
    import pyarrow as pa
    import pyarrow.parquet as pq
    from tcga2hf.schema import TABULAR_TABLES
    from tcga2hf_pipeline import tabular

    def _write(proj: str, scores: list[float], sample_types: list[str]) -> None:
        rows = [
            {
                "case_id": f"c{i}",
                "case_submitter_id": f"C{i}",
                "sample_id": f"s{i}",
                "sample_submitter_id": f"S{i}",
                "sample_type": st,
                "aliquot_id": f"a{i}",
                "aliquot_submitter_id": f"A{i}",
                "source_file_id": f"f{i}",
                "pathway": "P1",
                "matched_gene_count": 10,
                "original_gene_count": 12,
                "score_raw": v,
            }
            for i, (v, st) in enumerate(zip(scores, sample_types, strict=True))
        ]
        out = tmp_path / proj / "ssgsea_scores_hallmark" / "data.parquet"
        out.parent.mkdir(parents=True, exist_ok=True)
        pq.write_table(
            pa.Table.from_pylist(rows, schema=TABULAR_TABLES["ssgsea_scores_hallmark"]), out
        )

    _write("TCGA-AA", [1.0, 3.0], ["Primary Tumor", "Solid Tissue Normal"])
    _write("TCGA-BB", [5.0, 7.0], ["Primary Tumor", "Primary Tumor"])

    by_project = tabular.ssgsea_stats_rows(tmp_path, "hallmark")
    assert set(by_project) == {"TCGA-AA", "TCGA-BB"}

    aa = by_project["TCGA-AA"]
    pan = [r for r in aa if r["population"] == "pan_cancer" and r["sample_type"] is None]
    own = [r for r in aa if r["population"] == "project" and r["sample_type"] is None]
    assert len(pan) == 1 and len(own) == 1
    # pan-cancer spans both projects; the project row sees only its own.
    assert pan[0]["n_aliquots"] == 4 and pan[0]["mean"] == pytest.approx(4.0)
    assert own[0]["n_aliquots"] == 2 and own[0]["mean"] == pytest.approx(2.0)
    assert own[0]["project_id"] == "TCGA-AA" and pan[0]["project_id"] is None
    # Both projects must carry byte-equal pan-cancer rows.
    assert pan == [
        r
        for r in by_project["TCGA-BB"]
        if r["population"] == "pan_cancer" and r["sample_type"] is None
    ]


def test_stats_expose_the_gsva_divisor(tmp_path: Path) -> None:
    """MAX(max) - MIN(min) over pan_cancer rows must equal GSVA's divisor.

    This is what lets us publish raw scores only: the normalized view stays
    recoverable from the dataset itself rather than from prose in the card.
    """
    import numpy as np
    import pyarrow as pa
    import pyarrow.parquet as pq
    from tcga2hf.schema import TABULAR_TABLES
    from tcga2hf_pipeline import tabular

    rng = np.random.default_rng(0)
    rows = []
    for i in range(20):
        for p in ("P1", "P2"):
            rows.append(
                {
                    "case_id": f"c{i}",
                    "case_submitter_id": None,
                    "sample_id": None,
                    "sample_submitter_id": None,
                    "sample_type": "Primary Tumor",
                    "aliquot_id": f"a{i}",
                    "aliquot_submitter_id": None,
                    "source_file_id": None,
                    "pathway": p,
                    "matched_gene_count": 10,
                    "original_gene_count": 10,
                    "score_raw": float(rng.normal(0, 100)),
                }
            )
    out = tmp_path / "TCGA-XX" / "ssgsea_scores_hallmark" / "data.parquet"
    out.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.Table.from_pylist(rows, schema=TABULAR_TABLES["ssgsea_scores_hallmark"]), out)

    raw = np.array([r["score_raw"] for r in rows])
    _, expected = ssgsea.normalize_global(raw.reshape(2, -1))

    stats = tabular.ssgsea_stats_rows(tmp_path, "hallmark")["TCGA-XX"]
    pan = [r for r in stats if r["population"] == "pan_cancer" and r["sample_type"] is None]
    assert max(r["max"] for r in pan) - min(r["min"] for r in pan) == pytest.approx(expected)


def test_stats_sd_is_null_for_single_observation(tmp_path: Path) -> None:
    """ddof=1 is undefined at n=1; emit null rather than a misleading 0.0."""
    import pyarrow as pa
    import pyarrow.parquet as pq
    from tcga2hf.schema import TABULAR_TABLES
    from tcga2hf_pipeline import tabular

    rows = [
        {
            "case_id": "c",
            "case_submitter_id": None,
            "sample_id": None,
            "sample_submitter_id": None,
            "sample_type": "Primary Tumor",
            "aliquot_id": "a",
            "aliquot_submitter_id": None,
            "source_file_id": None,
            "pathway": "P1",
            "matched_gene_count": 10,
            "original_gene_count": 10,
            "score_raw": 1.0,
        }
    ]
    out = tmp_path / "TCGA-XX" / "ssgsea_scores_hallmark" / "data.parquet"
    out.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.Table.from_pylist(rows, schema=TABULAR_TABLES["ssgsea_scores_hallmark"]), out)

    stats = tabular.ssgsea_stats_rows(tmp_path, "hallmark")["TCGA-XX"]
    assert all(r["sd"] is None for r in stats)
    assert all(r["n_aliquots"] == 1 for r in stats)


def test_stats_rows_absent_when_no_scores(tmp_path: Path) -> None:
    from tcga2hf_pipeline import tabular

    assert tabular.ssgsea_stats_rows(tmp_path, "hallmark") == {}


def test_no_maximum_size_by_default(tmp_path: Path) -> None:
    """A very large gene set must survive the default filter.

    maxSize=500 is inherited GSEA-desktop convention, not a GSVA default,
    and dropping a set at build time is unrecoverable. Consumers filter on
    the gene counts we ship instead.
    """
    gmt = tmp_path / "s.gmt"
    gmt.write_text("HUGE\tdesc\t" + "\t".join(f"G{i}" for i in range(1200)) + "\n")
    genes = [f"G{i}" for i in range(2000)]
    kept, _ = ssgsea.map_gene_sets(ssgsea.load_gmt(gmt), genes)
    assert kept["HUGE"].size == 1200
    # ...but an explicit ceiling is still honoured when asked for.
    kept2, _ = ssgsea.map_gene_sets(ssgsea.load_gmt(gmt), genes, max_size=500)
    assert "HUGE" not in kept2
