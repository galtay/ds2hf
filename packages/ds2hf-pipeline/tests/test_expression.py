from __future__ import annotations

import json
from pathlib import Path

import pytest
from ds2hf.tcga.schema import EXPRESSION_FIELDS
from ds2hf_pipeline.tcga import expression


def _build_synthetic_expression_project(tmp_path: Path) -> Path:
    """Synthesize a tiny GDC-style expression dir: one TSV + manifest."""
    project_dir = tmp_path / "TCGA-XYZ"
    expr_dir = project_dir / "expression"
    expr_dir.mkdir(parents=True)

    # GDC TSV layout: 1 `# gene-model` line, header, 4 N_* QC rows, then gene rows.
    # N_* rows have empty gene_name / gene_type and only the first 3 count cols.
    tsv = (
        "# gene-model: GENCODE v36\n"
        "gene_id\tgene_name\tgene_type\tunstranded\tstranded_first\tstranded_second"
        "\ttpm_unstranded\tfpkm_unstranded\tfpkm_uq_unstranded\n"
        "N_unmapped\t\t\t1000\t1000\t1000\t\t\t\n"
        "N_multimapping\t\t\t2000\t2000\t2000\t\t\t\n"
        "N_noFeature\t\t\t3000\t3500\t3700\t\t\t\n"
        "N_ambiguous\t\t\t4000\t900\t910\t\t\t\n"
        "ENSG00000000003.15\tTSPAN6\tprotein_coding\t1375\t692\t683\t31.4\t12.1\t15.6\n"
        "ENSG00000000005.6\tTNMD\tprotein_coding\t1\t1\t0\t0.07\t0.027\t0.034\n"
        "ENSG00000000419.13\tDPM1\tprotein_coding\t500\t250\t250\t10.5\t4.2\t5.3\n"
    )
    file_name = "synthetic.rna_seq.augmented_star_gene_counts.tsv"
    (expr_dir / file_name).write_text(tsv)

    manifest = [
        {
            "file_id": "synthetic-expr-file-id",
            "file_name": file_name,
            "cases": [
                {
                    "case_id": "case-uuid",
                    "submitter_id": "TCGA-XX-1",
                    "samples": [
                        {
                            "sample_id": "tumor-sample-id",
                            "portions": [
                                {"analytes": [{"aliquots": [{"aliquot_id": "rna-aliquot-uuid"}]}]}
                            ],
                        }
                    ],
                }
            ],
        }
    ]
    (expr_dir / "manifest.json").write_text(json.dumps(manifest))
    return project_dir


def test_load_for_project_parses_qc_and_genes(tmp_path: Path) -> None:
    project_dir = _build_synthetic_expression_project(tmp_path)
    by_case = expression.load_for_project(project_dir)
    assert set(by_case) == {"case-uuid"}
    records = by_case["case-uuid"]
    assert len(records) == 1
    rec = records[0]

    # FK + provenance (sample_id resolved later in attach)
    assert rec["aliquot_id"] == "rna-aliquot-uuid"
    assert rec["source_file_id"] == "synthetic-expr-file-id"
    assert rec["sample_id"] is None

    # QC scalars lifted from N_* rows (using the unstranded count)
    assert rec["N_unmapped"] == 1000
    assert rec["N_multimapping"] == 2000
    assert rec["N_noFeature"] == 3000
    assert rec["N_ambiguous"] == 4000

    # Gene metadata — only the 3 ENSG rows, N_* rows excluded
    assert rec["gene_id"] == ["ENSG00000000003.15", "ENSG00000000005.6", "ENSG00000000419.13"]
    assert rec["gene_name"] == ["TSPAN6", "TNMD", "DPM1"]
    assert all(t == "protein_coding" for t in rec["gene_type"])

    # Numeric arrays — index-aligned to gene_id
    assert rec["unstranded"] == [1375, 1, 500]
    assert rec["tpm_unstranded"][0] == pytest.approx(31.4)
    assert rec["fpkm_uq_unstranded"][2] == pytest.approx(5.3)


def test_expression_schema_drops_stranded_columns() -> None:
    # Schema must NOT include stranded_first / stranded_second
    names = {f.name for f in EXPRESSION_FIELDS}
    assert "stranded_first" not in names
    assert "stranded_second" not in names
    # But must include all the unstranded value columns
    for required in [
        "unstranded",
        "tpm_unstranded",
        "fpkm_unstranded",
        "fpkm_uq_unstranded",
        "N_unmapped",
        "N_multimapping",
        "N_noFeature",
        "N_ambiguous",
        "gene_id",
        "gene_name",
        "gene_type",
        "sample_id",
        "aliquot_id",
        "source_file_id",
    ]:
        assert required in names


def test_attach_resolves_sample_id_from_patient_samples(tmp_path: Path) -> None:
    project_dir = _build_synthetic_expression_project(tmp_path)
    by_case = expression.load_for_project(project_dir)

    rows = [
        {
            "case_id": "case-uuid",
            "samples": [
                {
                    "sample_id": "tumor-sample-id",
                    "portions": [
                        {"analytes": [{"aliquots": [{"aliquot_id": "rna-aliquot-uuid"}]}]}
                    ],
                }
            ],
        },
        {"case_id": "no-rna-here", "samples": []},
    ]
    expression.attach(rows, by_case)

    assert len(rows[0]["samples_gene_expression_quantification"]) == 1
    assert rows[0]["samples_gene_expression_quantification"][0]["sample_id"] == "tumor-sample-id"
    assert rows[1]["samples_gene_expression_quantification"] == []


def test_attach_leaves_sample_id_null_when_aliquot_unknown(tmp_path: Path) -> None:
    project_dir = _build_synthetic_expression_project(tmp_path)
    by_case = expression.load_for_project(project_dir)

    rows = [{"case_id": "case-uuid", "samples": []}]
    expression.attach(rows, by_case)

    rec = rows[0]["samples_gene_expression_quantification"][0]
    assert rec["sample_id"] is None
    # But other fields are still populated
    assert rec["aliquot_id"] == "rna-aliquot-uuid"
    assert len(rec["gene_id"]) == 3


def test_load_for_project_missing_dir_returns_empty(tmp_path: Path) -> None:
    assert expression.load_for_project(tmp_path / "nonexistent") == {}
