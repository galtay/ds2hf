from __future__ import annotations

import json
import shutil
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from tcga2hf.schema import TABULAR_TABLES
from tcga2hf_pipeline import gene_expression_quantification as ed

GENES = ["ENSG00000000003.15", "ENSG00000000005.6", "ENSG00000000419.13"]


def _write_project(
    root: Path,
    project_id: str,
    aliquots: list[str],
    *,
    gene_ids: list[str] | None = None,
) -> None:
    """A minimal built-project tree: gene_model + expression + cases.json.

    Values are deterministic functions of (aliquot position, gene position)
    so a test can assert the matrix landed transposed correctly rather than
    merely landing.
    """
    gene_ids = gene_ids or GENES
    project_dir = root / "processed_project_tabular" / project_id

    model = pa.table(
        {
            "gene_id": gene_ids,
            "gene_name": [f"G{i}" for i in range(len(gene_ids))],
            # A real null, as the chrM genes have — this must survive.
            "gene_type": ["protein_coding"] * (len(gene_ids) - 1) + [None],
            "chromosome": ["chr1"] * (len(gene_ids) - 1) + [None],
            "start": [100 * i for i in range(len(gene_ids))],
            "end": [100 * i + 50 for i in range(len(gene_ids))],
        }
    )
    (project_dir / "gene_model").mkdir(parents=True)
    pq.write_table(model, project_dir / "gene_model" / "data.parquet")

    rows = []
    for a_i, aliquot in enumerate(aliquots):
        for g_i, gene_id in enumerate(gene_ids):
            base = a_i * 100 + g_i
            rows.append(
                {
                    "case_id": f"case-{a_i}",
                    "case_submitter_id": f"TCGA-XX-{a_i:04d}",
                    "aliquot_id": aliquot,
                    "aliquot_submitter_id": f"{aliquot}-sub",
                    "source_file_id": f"file-{a_i}",
                    "gene_id": gene_id,
                    "unstranded": base,
                    "stranded_first": base + 1,
                    "stranded_second": base + 2,
                    "tpm_unstranded": base + 0.5,
                    "fpkm_unstranded": base + 0.25,
                    "fpkm_uq_unstranded": base + 0.75,
                }
            )
    (project_dir / "gene_expression_quantification").mkdir(parents=True)
    pq.write_table(
        pa.Table.from_pylist(rows, schema=TABULAR_TABLES["gene_expression_quantification"]),
        project_dir / "gene_expression_quantification" / "data.parquet",
    )

    cases = [
        {
            "case_id": f"case-{a_i}",
            "submitter_id": f"TCGA-XX-{a_i:04d}",
            "samples": [
                {
                    "sample_id": f"sample-{a_i}",
                    "submitter_id": f"TCGA-XX-{a_i:04d}-01A",
                    "sample_type": "Solid Tissue Normal" if a_i % 2 else "Primary Tumor",
                    "portions": [
                        {
                            "analytes": [
                                {
                                    "aliquots": [
                                        {
                                            "aliquot_id": aliquot,
                                            "submitter_id": f"{aliquot}-sub",
                                        }
                                    ]
                                }
                            ]
                        }
                    ],
                }
            ],
        }
        for a_i, aliquot in enumerate(aliquots)
    ]
    raw = root / "raw" / project_id
    raw.mkdir(parents=True)
    (raw / "cases.json").write_text(json.dumps(cases))

    # Raw STAR TSVs, needed because the N_* tallies are read from `raw/`
    # rather than from the per-project table. Only the head is parsed, so a
    # single gene row after the QC block is enough to exercise the early
    # break.
    expr = raw / "expression"
    expr.mkdir()
    manifest = []
    for a_i, aliquot in enumerate(aliquots):
        name = f"{aliquot}.rna_seq.augmented_star_gene_counts.tsv"
        header = (
            "gene_id\tgene_name\tgene_type\tunstranded\tstranded_first"
            "\tstranded_second\ttpm_unstranded\tfpkm_unstranded\tfpkm_uq_unstranded"
        )
        base = (a_i + 1) * 1000
        (expr / name).write_text(
            "# gene-model: GENCODE v36\n"
            + header
            + "\n"
            + f"N_unmapped\t\t\t{base}\t{base}\t{base}\t\t\t\n"
            + f"N_multimapping\t\t\t{base + 1}\t{base + 1}\t{base + 1}\t\t\t\n"
            + f"N_noFeature\t\t\t{base + 2}\t{base + 2}\t{base + 2}\t\t\t\n"
            + f"N_ambiguous\t\t\t{base + 3}\t{base + 3}\t{base + 3}\t\t\t\n"
            + f"{gene_ids[0]}\tG0\tprotein_coding\t7\t3\t4\t1.0\t0.5\t0.6\n"
        )
        manifest.append(
            {
                "file_id": f"file-{a_i}",
                "file_name": name,
                "cases": [
                    {
                        "case_id": f"case-{a_i}",
                        "samples": [
                            {"portions": [{"analytes": [{"aliquots": [{"aliquot_id": aliquot}]}]}]}
                        ],
                    }
                ],
            }
        )
    (expr / "manifest.json").write_text(json.dumps(manifest))


def samples_have_balance(out: Path) -> bool:
    df = pq.read_table(out / "samples" / "data.parquet").to_pandas()
    return "strand_balance" in df.columns and df.strand_balance.notna().all()


def test_values_column_round_trips_a_matrix() -> None:
    matrix = np.arange(12, dtype=np.float32).reshape(3, 4)
    col = ed._values_column(matrix, pa.float32())
    assert len(col) == 3
    np.testing.assert_array_equal(np.stack(col.to_numpy(zero_copy_only=False)), matrix)


def test_gene_axis_keeps_order_and_nulls(tmp_path: Path) -> None:
    _write_project(tmp_path, "TCGA-AA", ["al-0", "al-1"])
    table = ed.gene_axis(tmp_path / "processed_project_tabular" / "TCGA-AA")

    assert table.column("gene_id").to_pylist() == GENES
    assert table.column("gene_index").to_pylist() == [0, 1, 2]
    # The null chromosome must stay null rather than becoming a NaN string.
    assert table.column("chromosome").to_pylist()[-1] is None
    assert table.column("gene_type").to_pylist()[-1] is None


def test_build_aligns_axes_positionally(tmp_path: Path) -> None:
    """genes[i] <-> values[i], and samples[j] <-> row j, across projects."""
    _write_project(tmp_path, "TCGA-AA", ["al-a0", "al-a1"])
    _write_project(tmp_path, "TCGA-BB", ["al-b0"])
    out = tmp_path / "processed_gene_expression_quantification"

    counts, strand = ed.build(tmp_path / "processed_project_tabular", tmp_path / "raw", out)
    assert counts["samples"] == 3
    assert counts["genes"] == len(GENES)
    assert counts["tpm_unstranded"] == 3
    assert strand["n_samples"] == 3
    # Every sample carries its own measured balance, not just a cohort stat.
    assert samples_have_balance(out)

    samples = pq.read_table(out / "samples" / "data.parquet").to_pandas()
    # Projects are visited in sorted order, so sample_index is 0..n-1 with
    # TCGA-AA's two aliquots first.
    assert samples.sample_index.tolist() == [0, 1, 2]
    assert samples.project_id.tolist() == ["TCGA-AA", "TCGA-AA", "TCGA-BB"]
    assert samples.aliquot_id.tolist() == ["al-a0", "al-a1", "al-b0"]
    assert samples.sample_type.tolist() == [
        "Primary Tumor",
        "Solid Tissue Normal",
        "Primary Tumor",
    ]

    for name in ed.QUANTIFICATIONS:
        table = pq.read_table(out / name / "data.parquet")
        matrix = np.stack(table.column("values").to_numpy(zero_copy_only=False))
        assert matrix.shape == (3, len(GENES))
        # Each project restarts its aliquot numbering, so row 2 (TCGA-BB's
        # only aliquot) carries the a_i=0 values, not a_i=2's.
        offset = {
            "unstranded": 0,
            "stranded_first": 1,
            "stranded_second": 2,
            "tpm_unstranded": 0.5,
            "fpkm_unstranded": 0.25,
            "fpkm_uq_unstranded": 0.75,
        }[name]
        expected = np.array(
            [
                [0 + offset, 1 + offset, 2 + offset],
                [100 + offset, 101 + offset, 102 + offset],
                [0 + offset, 1 + offset, 2 + offset],
            ]
        )
        np.testing.assert_allclose(matrix, expected, rtol=1e-6)
        # The inline label columns must agree with the `samples` config,
        # since a training loop trusts them instead of joining.
        assert table.column("sample_index").to_pylist() == [0, 1, 2]
        assert table.column("aliquot_id").to_pylist() == samples.aliquot_id.tolist()


def test_build_rejects_a_divergent_gene_model(tmp_path: Path) -> None:
    """A positional gene axis is only sound if every project shares the model."""
    _write_project(tmp_path, "TCGA-AA", ["al-a0"])
    _write_project(tmp_path, "TCGA-BB", ["al-b0"], gene_ids=list(reversed(GENES)))

    with pytest.raises(ValueError, match="gene_model differs"):
        ed.build(
            tmp_path / "processed_project_tabular",
            tmp_path / "raw",
            tmp_path / "processed_gene_expression_quantification",
        )


def test_build_needs_a_built_project(tmp_path: Path) -> None:
    (tmp_path / "processed_project_tabular").mkdir()
    with pytest.raises(ValueError, match="no built projects"):
        ed.build(
            tmp_path / "processed_project_tabular",
            tmp_path / "raw",
            tmp_path / "processed_gene_expression_quantification",
        )


def test_strand_summary_flags_an_unstranded_library() -> None:
    """A 50/50 split is the signature of a non-strand-specific protocol."""
    first = np.array([500.0, 1010.0])
    second = np.array([500.0, 990.0])

    stats = ed.strand_summary(first, second)
    assert stats["n_samples"] == 2
    assert stats["balance_median"] == pytest.approx(0.5, abs=0.01)
    assert stats["n_strand_specific"] == 0


def test_strand_summary_flags_a_stranded_library() -> None:
    """A genuinely stranded library puts nearly every read on one side."""
    stats = ed.strand_summary(np.array([980.0, 20.0]), np.array([20.0, 980.0]))
    assert stats["n_strand_specific"] == 2


def test_strand_summary_survives_an_empty_library() -> None:
    """A sample with no reads must not produce a divide-by-zero NaN."""
    stats = ed.strand_summary(np.array([0.0, 500.0]), np.array([0.0, 500.0]))
    assert stats["n_samples"] == 1
    assert stats["balance_median"] == pytest.approx(0.5)


def test_strand_summary_counts_a_single_outlier() -> None:
    """One lopsided library is one strand-specific sample."""
    stats = ed.strand_summary(np.array([980.0]), np.array([20.0]))
    assert stats["n_strand_specific"] == 1


def test_build_reports_projects_in_the_strand_breakdown(tmp_path: Path) -> None:
    """The per-project accounting must cover every project.

    Guards a real bug: the breakdown was silently absent and the card
    rendered "0 of the 0 projects", contradicting the outlier count printed
    two lines above it.
    """
    _write_project(tmp_path, "TCGA-AA", ["al-a0", "al-a1"])
    _write_project(tmp_path, "TCGA-BB", ["al-b0"])

    _, strand = ed.build(
        tmp_path / "processed_project_tabular",
        tmp_path / "raw",
        tmp_path / "processed_gene_expression_quantification",
    )
    assert "outlier_projects" in strand
    assert "n_projects_clean" in strand
    covered = strand["n_projects_clean"] + len(strand["outlier_projects"])
    assert covered == 2, "every project must appear in exactly one bucket"
    # And the two accountings of the same fact must agree.
    assert sum(n for _, n, _ in strand["outlier_projects"]) == strand["n_strand_specific"]


def test_qc_tallies_read_only_the_head(tmp_path: Path) -> None:
    """The N_* rows come off the top of each TSV, keyed by aliquot."""
    _write_project(tmp_path, "TCGA-AA", ["al-a0", "al-a1"])
    tallies = ed.qc_tallies(tmp_path / "raw" / "TCGA-AA")

    assert set(tallies) == {"al-a0", "al-a1"}
    assert tallies["al-a0"] == {
        "n_unmapped": 1000,
        "n_multimapping": 1001,
        "n_nofeature": 1002,
        "n_ambiguous": 1003,
    }
    assert tallies["al-a1"]["n_unmapped"] == 2000


def test_qc_tallies_are_null_without_raw_files(tmp_path: Path) -> None:
    """A missing library is not a library of zero unmapped reads."""
    _write_project(tmp_path, "TCGA-AA", ["al-a0"])
    shutil.rmtree(tmp_path / "raw" / "TCGA-AA" / "expression")

    _, _ = ed.build(
        tmp_path / "processed_project_tabular",
        tmp_path / "raw",
        tmp_path / "processed_gene_expression_quantification",
    )
    df = pq.read_table(
        tmp_path / "processed_gene_expression_quantification" / "samples" / "data.parquet"
    ).to_pandas()
    assert df.n_unmapped.isna().all()


def test_build_carries_the_qc_tallies(tmp_path: Path) -> None:
    _write_project(tmp_path, "TCGA-AA", ["al-a0", "al-a1"])
    _write_project(tmp_path, "TCGA-BB", ["al-b0"])
    out = tmp_path / "processed_gene_expression_quantification"

    ed.build(tmp_path / "processed_project_tabular", tmp_path / "raw", out)
    df = pq.read_table(out / "samples" / "data.parquet").to_pandas()
    assert df.n_unmapped.tolist() == [1000, 2000, 1000]
    assert df.n_ambiguous.tolist() == [1003, 2003, 1003]


# --- verification -----------------------------------------------------------
# Each of these breaks the built tree in one specific way. A check that only
# ever passes tells you nothing.


def _built(tmp_path: Path) -> tuple[Path, Path]:
    """A two-project tree, built. Returns (expression_dir, project_dir)."""
    _write_project(tmp_path, "TCGA-AA", ["al-a0", "al-a1"])
    _write_project(tmp_path, "TCGA-BB", ["al-b0"])
    out = tmp_path / "processed_gene_expression_quantification"
    ed.build(tmp_path / "processed_project_tabular", tmp_path / "raw", out)
    return out, tmp_path / "processed_project_tabular"


def test_row_group_lookup_handles_uneven_groups() -> None:
    """Row groups are per-project batches, so index // size is wrong."""
    from tcga2hf_pipeline.verify import _row_group_for

    class FakeMeta:
        num_row_groups = 3

        def row_group(self, i):
            return type("RG", (), {"num_rows": [64, 15, 64][i]})()

    fake = type("PF", (), {"metadata": FakeMeta()})()
    assert _row_group_for(fake, 0) == (0, 0)
    assert _row_group_for(fake, 63) == (0, 63)
    assert _row_group_for(fake, 64) == (1, 0)  # a fixed stride would say (1, 0) too
    assert _row_group_for(fake, 79) == (2, 0)  # but here it would wrongly say (1, 15)
    assert _row_group_for(fake, 100) == (2, 21)


def test_axis_alignment_passes_on_a_good_tree(tmp_path: Path) -> None:
    from tcga2hf_pipeline.verify import check_axis_alignment

    out, _ = _built(tmp_path)
    assert check_axis_alignment(out).passed


def test_axis_alignment_catches_a_reordered_config(tmp_path: Path) -> None:
    """The failure the positional contract exists to prevent."""
    from tcga2hf_pipeline.verify import check_axis_alignment

    out, _ = _built(tmp_path)
    path = out / "tpm_unstranded" / "data.parquet"
    table = pq.read_table(path)
    pq.write_table(table.take([2, 1, 0]), path)

    check = check_axis_alignment(out)
    assert not check.passed
    assert any("aliquot order differs" in d for d in check.details)


def test_gene_axis_catches_a_divergent_model(tmp_path: Path) -> None:

    from tcga2hf_pipeline.verify import check_gene_axis

    out, projects = _built(tmp_path)
    assert check_gene_axis(out, projects).passed

    genes = out / "genes" / "data.parquet"
    table = pq.read_table(genes)
    swapped = table.set_column(
        table.schema.get_field_index("gene_id"),
        "gene_id",
        pa.array(list(reversed(table.column("gene_id").to_pylist()))),
    )
    pq.write_table(swapped, genes)
    assert not check_gene_axis(out, projects).passed


def test_sample_metadata_catches_a_null_label(tmp_path: Path) -> None:

    from tcga2hf_pipeline.verify import check_sample_metadata

    out, _ = _built(tmp_path)
    assert check_sample_metadata(out).passed

    path = out / "samples" / "data.parquet"
    table = pq.read_table(path)
    nulled = table.set_column(
        table.schema.get_field_index("sample_type"),
        "sample_type",
        pa.array([None] * table.num_rows, type=pa.string()),
    )
    pq.write_table(nulled, path)

    check = check_sample_metadata(out)
    assert not check.passed
    assert any("sample_type" in d for d in check.details)


def test_expression_values_catches_a_corrupted_cell(tmp_path: Path) -> None:

    from tcga2hf_pipeline.verify import check_expression_values

    out, projects = _built(tmp_path)
    assert check_expression_values(out, projects, genes_per_sample=3).passed

    path = out / "tpm_unstranded" / "data.parquet"
    table = pq.read_table(path)
    matrix = np.stack(table.column("values").to_numpy(zero_copy_only=False))
    matrix[0, 0] += 999.0
    patched = table.set_column(
        table.schema.get_field_index("values"),
        "values",
        ed._values_column(matrix.astype(np.float32), pa.float32()),
    )
    pq.write_table(patched, path)

    check = check_expression_values(out, projects, genes_per_sample=3)
    assert not check.passed
    assert check.details
