from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from ds2hf.tcga.schema import TABULAR_TABLES
from ds2hf_pipeline.tcga import gene_expression_quantification as ed

GENES = ["ENSG00000000003.15", "ENSG00000000005.6", "ENSG00000000419.13"]

STATUS = {
    "data_release": "Data Release 46.0 - August 10, 2026",
    "data_release_version": {"major": 46, "minor": 0, "release_date": "2026-08-10"},
}


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
            "start": [100 * i for i in range(len(gene_ids) - 1)] + [None],
            "end": [100 * i + 50 for i in range(len(gene_ids) - 1)] + [None],
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
                    "case_id": f"case-{project_id}-{a_i}",
                    "case_submitter_id": f"TCGA-XX-{a_i:04d}",
                    "aliquot_id": aliquot,
                    "aliquot_submitter_id": f"{aliquot}-sub",
                    "source_file_id": f"file-{project_id}-{a_i}",
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
            "case_id": f"case-{project_id}-{a_i}",
            "submitter_id": f"TCGA-XX-{a_i:04d}",
            "samples": [
                {
                    "sample_id": f"sample-{project_id}-{a_i}",
                    "submitter_id": f"TCGA-XX-{a_i:04d}-01A",
                    "sample_type": "Solid Tissue Normal" if a_i % 2 else "Primary Tumor",
                    "portions": [
                        {
                            "portion_id": f"portion-{project_id}-{a_i}",
                            "submitter_id": f"TCGA-XX-{a_i:04d}-01A-01",
                            "analytes": [
                                {
                                    "analyte_id": f"analyte-{project_id}-{a_i}",
                                    # Aliquot 1's analyte barcode does not
                                    # nest, as GDC's -12R / -11H ones do not.
                                    "submitter_id": f"TCGA-XX-{a_i:04d}-01A-"
                                    + ("11H" if a_i == 1 else "01R"),
                                    "analyte_type": "RNA",
                                    "aliquots": [
                                        {
                                            "aliquot_id": aliquot,
                                            "submitter_id": f"TCGA-XX-{a_i:04d}-01A-01R-A000-07",
                                        }
                                    ],
                                }
                            ],
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

    # Raw STAR TSVs, the GDC files the per-project rows above were read
    # from: the same values, preceded by the N_* tally block. The tallies
    # are read from here, and verification re-reads whole rows from here.
    expr = raw / "expression"
    expr.mkdir()
    (expr / "gdc_status.json").write_text(json.dumps(STATUS))
    manifest = []
    for a_i, aliquot in enumerate(aliquots):
        name = f"{aliquot}.rna_seq.augmented_star_gene_counts.tsv"
        header = (
            "gene_id\tgene_name\tgene_type\tunstranded\tstranded_first"
            "\tstranded_second\ttpm_unstranded\tfpkm_unstranded\tfpkm_uq_unstranded"
        )
        tally = (a_i + 1) * 1000
        gene_rows = "".join(
            f"{gene_id}\tG{g_i}\t{'protein_coding' if g_i < len(gene_ids) - 1 else ''}"
            f"\t{a_i * 100 + g_i}\t{a_i * 100 + g_i + 1}\t{a_i * 100 + g_i + 2}"
            f"\t{a_i * 100 + g_i + 0.5}\t{a_i * 100 + g_i + 0.25}\t{a_i * 100 + g_i + 0.75}\n"
            for g_i, gene_id in enumerate(gene_ids)
        )
        content = (
            "# gene-model: GENCODE v36\n"
            + header
            + "\n"
            + f"N_unmapped\t\t\t{tally}\t{tally}\t{tally}\t\t\t\n"
            + f"N_multimapping\t\t\t{tally + 1}\t{tally + 1}\t{tally + 1}\t\t\t\n"
            + f"N_noFeature\t\t\t{tally + 2}\t{tally + 5}\t{tally + 7}\t\t\t\n"
            + f"N_ambiguous\t\t\t{tally + 3}\t{tally + 3}\t{tally + 3}\t\t\t\n"
            + gene_rows
        ).encode()
        (expr / name).write_bytes(content)
        manifest.append(
            {
                "file_id": f"file-{project_id}-{a_i}",
                "file_name": name,
                "md5sum": hashlib.md5(content).hexdigest(),
                "gdc_version": "1",
                "gdc_first_release": "32.0",
                "cases": [
                    {
                        "case_id": f"case-{project_id}-{a_i}",
                        "samples": [
                            {"portions": [{"analytes": [{"aliquots": [{"aliquot_id": aliquot}]}]}]}
                        ],
                    }
                ],
            }
        )
    (expr / "manifest.json").write_text(json.dumps(manifest))

    # A gene-level copy number TSV, the source of the `genes` coordinates.
    # The last gene is absent, as the chrM genes are in the real files.
    cnv = raw / "gene_level_copy_number"
    cnv.mkdir()
    (cnv / "one.gene_level.copy_number_variation.tsv").write_text(
        "gene_id\tgene_name\tchromosome\tstart\tend\tcopy_number\n"
        + "".join(
            f"{gene_id}\tG{g_i}\tchr1\t{100 * g_i}\t{100 * g_i + 50}\t2\n"
            for g_i, gene_id in enumerate(gene_ids[:-1])
        )
    )


def aliquots_have_balance(out: Path) -> bool:
    df = pq.read_table(out / "aliquots" / "data.parquet").to_pandas()
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
    """genes[i] <-> values[i], and aliquots[j] <-> row j, across projects."""
    _write_project(tmp_path, "TCGA-AA", ["al-a0", "al-a1"])
    _write_project(tmp_path, "TCGA-BB", ["al-b0"])
    out = tmp_path / "processed_gene_expression_quantification"

    counts, strand = ed.build(tmp_path / "processed_project_tabular", tmp_path / "raw", out)
    assert counts["aliquots"] == 3
    assert counts["genes"] == len(GENES)
    assert counts["tpm_unstranded"] == 3
    assert strand["n_samples"] == 3
    # Every aliquot carries its own measured balance, not just a cohort stat.
    assert aliquots_have_balance(out)

    aliquots = pq.read_table(out / "aliquots" / "data.parquet").to_pandas()
    # Projects are visited in sorted order, so aliquot_index is 0..n-1 with
    # TCGA-AA's two aliquots first.
    assert aliquots.aliquot_index.tolist() == [0, 1, 2]
    assert aliquots.project_id.tolist() == ["TCGA-AA", "TCGA-AA", "TCGA-BB"]
    assert aliquots.aliquot_id.tolist() == ["al-a0", "al-a1", "al-b0"]
    assert aliquots.sample_id.tolist() == [
        "sample-TCGA-AA-0",
        "sample-TCGA-AA-1",
        "sample-TCGA-BB-0",
    ]
    assert aliquots.sample_type.tolist() == [
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
        # The inline label columns must agree with the `aliquots` config,
        # since a training loop trusts them instead of joining.
        assert table.column("aliquot_index").to_pylist() == [0, 1, 2]
        assert table.column("aliquot_id").to_pylist() == aliquots.aliquot_id.tolist()


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
    for tsv in (tmp_path / "raw" / "TCGA-AA" / "expression").glob("*.tsv"):
        tsv.unlink()

    _, _ = ed.build(
        tmp_path / "processed_project_tabular",
        tmp_path / "raw",
        tmp_path / "processed_gene_expression_quantification",
    )
    df = pq.read_table(
        tmp_path / "processed_gene_expression_quantification" / "aliquots" / "data.parquet"
    ).to_pandas()
    assert df.n_unmapped.isna().all()


def test_build_carries_the_qc_tallies(tmp_path: Path) -> None:
    _write_project(tmp_path, "TCGA-AA", ["al-a0", "al-a1"])
    _write_project(tmp_path, "TCGA-BB", ["al-b0"])
    out = tmp_path / "processed_gene_expression_quantification"

    ed.build(tmp_path / "processed_project_tabular", tmp_path / "raw", out)
    df = pq.read_table(out / "aliquots" / "data.parquet").to_pandas()
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
    from ds2hf_pipeline.tcga.verify import _row_group_for

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
    from ds2hf_pipeline.tcga.verify import check_axis_alignment

    out, _ = _built(tmp_path)
    assert check_axis_alignment(out).passed


def test_axis_alignment_catches_a_reordered_config(tmp_path: Path) -> None:
    """The failure the positional contract exists to prevent."""
    from ds2hf_pipeline.tcga.verify import check_axis_alignment

    out, _ = _built(tmp_path)
    path = out / "tpm_unstranded" / "data.parquet"
    table = pq.read_table(path)
    pq.write_table(table.take([2, 1, 0]), path)

    check = check_axis_alignment(out)
    assert not check.passed
    assert any("aliquot order differs" in d for d in check.details)


def test_gene_axis_catches_a_divergent_model(tmp_path: Path) -> None:
    from ds2hf_pipeline.tcga.verify import check_gene_axis

    out, _ = _built(tmp_path)
    assert check_gene_axis(out, tmp_path / "raw").passed

    genes = out / "genes" / "data.parquet"
    table = pq.read_table(genes)
    swapped = table.set_column(
        table.schema.get_field_index("gene_id"),
        "gene_id",
        pa.array(list(reversed(table.column("gene_id").to_pylist()))),
    )
    pq.write_table(swapped, genes)
    assert not check_gene_axis(out, tmp_path / "raw").passed


def test_aliquot_metadata_catches_a_null_label(tmp_path: Path) -> None:
    from ds2hf_pipeline.tcga.verify import check_aliquot_metadata

    out, _ = _built(tmp_path)
    assert check_aliquot_metadata(out).passed

    path = out / "aliquots" / "data.parquet"
    table = pq.read_table(path)
    nulled = table.set_column(
        table.schema.get_field_index("sample_type"),
        "sample_type",
        pa.array([None] * table.num_rows, type=pa.string()),
    )
    pq.write_table(nulled, path)

    check = check_aliquot_metadata(out)
    assert not check.passed
    assert any("sample_type" in d for d in check.details)


def test_source_bytes_passes_on_a_good_tree(tmp_path: Path) -> None:
    from ds2hf_pipeline.tcga.verify import check_source_bytes

    out, _ = _built(tmp_path)
    check = check_source_bytes(out, tmp_path / "raw", per_project=5)
    assert check.passed, check.details


def test_source_bytes_catches_a_corrupted_cell(tmp_path: Path) -> None:
    from ds2hf_pipeline.tcga.verify import check_source_bytes

    out, _ = _built(tmp_path)
    path = out / "tpm_unstranded" / "data.parquet"
    table = pq.read_table(path)
    matrix = np.stack(table.column("values").to_numpy(zero_copy_only=False))
    matrix[0, 0] += 1e-3  # far below the old relative tolerance, still wrong
    patched = table.set_column(
        table.schema.get_field_index("values"),
        "values",
        ed._values_column(matrix.astype(np.float32), pa.float32()),
    )
    pq.write_table(patched, path)

    check = check_source_bytes(out, tmp_path / "raw", per_project=5)
    assert not check.passed
    assert any("tpm_unstranded differs in 1 of" in d for d in check.details)


def test_source_bytes_catches_a_file_that_is_not_the_one_named(tmp_path: Path) -> None:
    """Bytes that do not hash to the row's md5 are not its source."""
    from ds2hf_pipeline.tcga.verify import check_source_bytes

    out, _ = _built(tmp_path)
    raw = next((tmp_path / "raw" / "TCGA-BB" / "expression").glob("*.tsv"))
    raw.write_bytes(raw.read_bytes() + b"\n")

    check = check_source_bytes(out, tmp_path / "raw", per_project=5)
    assert not check.passed
    assert any("md5" in d for d in check.details)


def test_source_bytes_catches_a_wrong_tally(tmp_path: Path) -> None:
    from ds2hf_pipeline.tcga.verify import check_source_bytes

    out, _ = _built(tmp_path)
    path = out / "aliquots" / "data.parquet"
    table = pq.read_table(path)
    index = table.schema.get_field_index("n_nofeature")
    # The stranded_first tally, which the TSV also carries, is the wrong one.
    wrong = pa.array([1005, 2005, 1005], type=pa.int64())
    pq.write_table(table.set_column(index, "n_nofeature", wrong), path)

    check = check_source_bytes(out, tmp_path / "raw", per_project=5)
    assert not check.passed
    assert any("n_nofeature" in d for d in check.details)


def test_aliquot_metadata_catches_a_url_for_another_file(tmp_path: Path) -> None:
    from ds2hf_pipeline.tcga.verify import check_aliquot_metadata

    out, _ = _built(tmp_path)
    path = out / "aliquots" / "data.parquet"
    table = pq.read_table(path)
    urls = table.column("source_file_url").to_pylist()
    urls[0], urls[1] = urls[1], urls[0]
    index = table.schema.get_field_index("source_file_url")
    pq.write_table(table.set_column(index, "source_file_url", pa.array(urls)), path)

    check = check_aliquot_metadata(out)
    assert not check.passed
    assert any("source_file_url" in d for d in check.details)


class _FakeGDC:
    """Answers the one `/files` query `check_gdc_current` makes."""

    def __init__(self, hits: list[dict]) -> None:
        self.hits = hits

    def files(self, **_kwargs) -> list[dict]:
        return self.hits


def _gdc_hits(out: Path) -> list[dict]:
    frame = pq.read_table(out / "aliquots" / "data.parquet").to_pandas()
    return [
        {"file_id": r.source_file_id, "md5sum": r.source_file_md5sum, "version": "1"}
        for r in frame.itertuples()
    ]


def test_gdc_current_passes_when_gdc_serves_the_same_files(tmp_path: Path) -> None:
    from ds2hf_pipeline.tcga.verify import check_gdc_current

    out, _ = _built(tmp_path)
    assert check_gdc_current(_FakeGDC(_gdc_hits(out)), out).passed


def test_gdc_current_reports_each_kind_of_drift(tmp_path: Path) -> None:
    from ds2hf_pipeline.tcga.verify import check_gdc_current

    out, _ = _built(tmp_path)
    hits = _gdc_hits(out)
    hits[0]["md5sum"] = "0" * 32  # replaced under the same id
    del hits[1]  # withdrawn
    hits.append({"file_id": "new-file", "md5sum": "1" * 32, "version": "1"})  # added

    check = check_gdc_current(_FakeGDC(hits), out)
    assert not check.passed
    assert "1 added, 1 withdrawn, 1 changed" in check.summary


def test_build_records_each_rows_source_file(tmp_path: Path) -> None:
    out, _ = _built(tmp_path)
    df = pq.read_table(out / "aliquots" / "data.parquet").to_pandas()
    assert df.source_file_id.tolist() == ["file-TCGA-AA-0", "file-TCGA-AA-1", "file-TCGA-BB-0"]
    assert df.source_file_url.tolist() == [
        f"https://api.gdc.cancer.gov/data/{f}" for f in df.source_file_id
    ]
    assert df.source_file_md5sum.str.len().eq(32).all()
    assert df.source_file_version.eq("1").all()
    assert df.source_file_first_release.eq("32.0").all()


def test_build_rejects_a_row_whose_file_is_not_in_the_manifest(tmp_path: Path) -> None:
    """A value whose source cannot be named cannot be verified."""
    _write_project(tmp_path, "TCGA-AA", ["al-a0"])
    manifest = tmp_path / "raw" / "TCGA-AA" / "expression" / "manifest.json"
    manifest.write_text("[]")

    with pytest.raises(ValueError, match="not in the expression manifest"):
        ed.build(
            tmp_path / "processed_project_tabular",
            tmp_path / "raw",
            tmp_path / "processed_gene_expression_quantification",
        )


def test_download_release_reads_the_expression_fetch(tmp_path: Path) -> None:
    _write_project(tmp_path, "TCGA-AA", ["al-a0"])
    # The project-level status belongs to the clinical fetch; it must not
    # be the one reported.
    (tmp_path / "raw" / "TCGA-AA" / "gdc_status.json").write_text(
        json.dumps({"data_release": "Data Release 45.0"})
    )
    release = ed.download_release(tmp_path / "raw", ["TCGA-AA"])
    assert release["data_release"] == STATUS["data_release"]


def test_download_release_rejects_mixed_releases(tmp_path: Path) -> None:
    _write_project(tmp_path, "TCGA-AA", ["al-a0"])
    _write_project(tmp_path, "TCGA-BB", ["al-b0"])
    status = tmp_path / "raw" / "TCGA-BB" / "expression" / "gdc_status.json"
    status.write_text(json.dumps({**STATUS, "data_release": "Data Release 47.0"}))

    with pytest.raises(ValueError, match="more than one GDC release"):
        ed.download_release(tmp_path / "raw", ["TCGA-AA", "TCGA-BB"])


def test_card_stands_alone(tmp_path: Path) -> None:
    """The card links GDC and the pipeline source, never a sibling dataset."""
    from ds2hf_pipeline.tcga import dataset_card

    out, _ = _built(tmp_path)
    dataset_card.write_expression_card(
        out, {"aliquots": 3, "genes": len(GENES)}, ["TCGA-AA", "TCGA-BB"], STATUS
    )
    text = (out / "README.md").read_text()
    assert "huggingface.co/datasets" not in text
    others = {
        m
        for m in __import__("re").findall(r"gabrielaltay/[\w-]+", text)
        if m != "gabrielaltay/tcga-gene-expression-quantification-open"
    }
    assert not others, f"card names other datasets: {others}"
    assert "-tabular-open" not in text


def test_card_states_what_it_measured(tmp_path: Path) -> None:
    from ds2hf_pipeline.tcga import dataset_card

    out, _ = _built(tmp_path)
    dataset_card.write_expression_card(
        out, {"aliquots": 3, "genes": len(GENES)}, ["TCGA-AA", "TCGA-BB"], STATUS
    )
    text = (out / "README.md").read_text()
    # The header's structural counts, read off the built table.
    assert "3 aliquots × 3 genes, from 3 patients in 2 TCGA projects" in text
    assert "GDC Data Release 46.0 (2026-08-10)" in text
    # Nothing that would need rechecking when the data changes.
    assert "sha256:" not in text
    assert "%" not in text.split("## Data dictionary")[0].replace("% ", "")


def test_generated_card_frontmatter_parses_as_yaml(tmp_path: Path) -> None:
    """The end-to-end guard: HF parses this block to decide what to serve.

    The frontmatter is isolated between the opening and closing `---`
    delimiters, and the closing one is asserted to exist. An earlier version
    split on "---\n" and took element [1], which silently returns the whole
    document when the closing delimiter is missing -- so a card whose
    frontmatter never terminated still parsed, still yielded the keys it
    checked, and shipped to the Hub as "empty or missing yaml metadata".
    """
    import yaml
    from ds2hf_pipeline.tcga import dataset_card

    out, _ = _built(tmp_path)
    counts = {"aliquots": 3, "genes": len(GENES)}
    dataset_card.write_expression_card(out, counts, ["TCGA-AA", "TCGA-BB"], STATUS)

    text = (out / "README.md").read_text()
    assert text.startswith("---\n"), "card must open with a frontmatter delimiter"
    closing = text.find("\n---\n", 3)
    assert closing != -1, "frontmatter has no closing --- on its own line"

    front = yaml.safe_load(text[4 : closing + 1])
    assert isinstance(front, dict), f"frontmatter is {type(front).__name__}, not a mapping"

    assert front["license"] == "other"
    assert front["pretty_name"].startswith("TCGA Gene Expression Quantification")
    assert isinstance(front["tags"], list) and "rna-seq" in front["tags"]

    declared = {c["config_name"] for c in front["configs"]}
    on_disk = {d.name for d in out.iterdir() if (d / "data.parquet").exists()}
    assert declared == on_disk, f"declared {declared} != on disk {on_disk}"

    # The body must start after the frontmatter, not be swallowed by it.
    assert text[closing:].lstrip("-\n").startswith("# TCGA")


def test_generated_card_has_well_formed_block_boundaries(tmp_path: Path) -> None:
    """Headings and link definitions each need a blank line before them.

    Markdown is whitespace-significant in ways that are invisible until
    rendered: a `##` heading glued to the previous line is not a heading,
    and a `[label]: url` that lands inside a paragraph is absorbed into it
    as literal text instead of defining a link. Both shipped -- the second
    put the whole reference block on screen as prose.
    """
    from ds2hf_pipeline.tcga import dataset_card

    out, _ = _built(tmp_path)
    dataset_card.write_expression_card(
        out, {"aliquots": 3, "genes": len(GENES)}, ["TCGA-AA"], STATUS
    )
    lines = (out / "README.md").read_text().split("\n")

    problems = []
    for i, line in enumerate(lines):
        if not i or lines[i - 1].strip() == "":
            continue
        is_definition = line.startswith("[") and "]: " in line
        # Definitions may stack; only the first of a run needs the blank.
        after_definition = lines[i - 1].startswith("[") and "]: " in lines[i - 1]
        if line.startswith("## ") or (is_definition and not after_definition):
            problems.append(f"line {i + 1}: {line[:60]!r} follows {lines[i - 1][:40]!r}")
    assert not problems, "markdown blocks need a blank line before them:\n  " + "\n  ".join(
        problems
    )


def test_card_dictionary_documents_every_column() -> None:
    """Every published column has a dictionary row, in schema order.

    A column added to the schema without a definition and a source would
    ship undocumented; this fails until the dictionary says what it is and
    where it comes from.
    """
    from ds2hf_pipeline.tcga.dataset_card import _EXPRESSION_DICTIONARY

    documented = {k: [row[0] for row in rows] for k, rows in _EXPRESSION_DICTIONARY.items()}
    assert documented["aliquots"] == [f.name for f in ed.ALIQUOTS_FIELDS]
    assert documented["genes"] == [f.name for f in ed.GENES_FIELDS]
    assert documented["measure"] == [f.name for f in ed.INLINE_ALIQUOT_FIELDS] + ["values"]


def test_card_dictionary_marks_exactly_the_computed_columns() -> None:
    """The card says everything not marked computed is GDC's, verbatim."""
    from ds2hf_pipeline.tcga.dataset_card import (
        _EXPRESSION_DICTIONARY,
        CARD_MODULE_COMPUTED_COLUMNS,
    )

    computed = {
        column
        for rows in _EXPRESSION_DICTIONARY.values()
        for column, _, source in rows
        if source == "computed"
    }
    assert computed == CARD_MODULE_COMPUTED_COLUMNS
    assert computed == {"aliquot_index", "gene_index", "strand_balance", "source_file_url"}


def test_card_dictionary_links_resolve_to_real_gdc_fields() -> None:
    """Every GDC dictionary link names an entity and field that exist.

    The GDC dictionary viewer is a single-page app, so a link to a missing
    field still returns 200 and simply shows nothing. Checked instead
    against the dictionary snapshot the fetch saves, when one is present.
    """
    import os

    from ds2hf_pipeline.tcga.dataset_card import _EXPRESSION_DICTIONARY, _EXPRESSION_SOURCES

    root = Path(os.environ.get("TCGA2HF_DATA_DIR", Path.home() / "data" / "tcga2hf"))
    snapshots = sorted((root / "raw").glob("gdc_dictionary.*.json"))
    if not snapshots:
        pytest.skip("no GDC dictionary snapshot under the data dir")
    dictionary = json.loads(snapshots[-1].read_text())

    missing = []
    for rows in _EXPRESSION_DICTIONARY.values():
        for column, _, source in rows:
            if source in _EXPRESSION_SOURCES:
                continue
            entity, _, field = source.partition(".")
            if entity not in dictionary or (
                field and field not in dictionary[entity]["properties"]
            ):
                missing.append(f"{column} -> {source}")
    assert not missing, f"not in {snapshots[-1].name}: {missing}"
