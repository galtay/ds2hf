from __future__ import annotations

import gzip
import json
from pathlib import Path

import pytest
from ds2hf.tcga.schema import _MAF_COLUMNS, MUTATION_FIELDS
from ds2hf_pipeline.tcga import mutations


def _build_mini_maf_tarball(tmp_path: Path) -> Path:
    """Synthesize a tiny MAF + manifest in a project-shaped directory."""
    project_dir = tmp_path / "TCGA-XYZ"
    mut_dir = project_dir / "mutations"
    mut_dir.mkdir(parents=True)

    # Write a minimal MAF: 4 `#` lines + header + 2 variants
    header = "\t".join(_MAF_COLUMNS)
    rows = []
    for i, (chrom, pos) in enumerate([("chr2", 200), ("chr1", 100)]):
        cells = [""] * len(_MAF_COLUMNS)
        cells[_MAF_COLUMNS.index("Hugo_Symbol")] = f"GENE{i}"
        cells[_MAF_COLUMNS.index("Entrez_Gene_Id")] = "1234"
        cells[_MAF_COLUMNS.index("Chromosome")] = chrom
        cells[_MAF_COLUMNS.index("Start_Position")] = str(pos)
        cells[_MAF_COLUMNS.index("End_Position")] = str(pos)
        cells[_MAF_COLUMNS.index("Variant_Classification")] = "Missense_Mutation"
        cells[_MAF_COLUMNS.index("Reference_Allele")] = "A"
        cells[_MAF_COLUMNS.index("Tumor_Seq_Allele2")] = "G"
        cells[_MAF_COLUMNS.index("t_depth")] = "100"
        cells[_MAF_COLUMNS.index("t_alt_count")] = "50"
        cells[_MAF_COLUMNS.index("gnomAD_AF")] = "0.001"
        cells[_MAF_COLUMNS.index("Tumor_Sample_UUID")] = "tumor-aliquot-uuid"
        cells[_MAF_COLUMNS.index("Matched_Norm_Sample_UUID")] = "normal-aliquot-uuid"
        cells[_MAF_COLUMNS.index("case_id")] = "case-uuid"
        rows.append("\t".join(cells))
    body = (
        "\n".join(
            [
                "#version gdc-1.0.0",
                "#filedate 20260504",
                "#tumor.aliquot x",
                "#normal.aliquot y",
                header,
                *rows,
            ]
        )
        + "\n"
    )

    maf_path = mut_dir / "synthetic.maf.gz"
    with gzip.open(maf_path, "wt") as fh:
        fh.write(body)

    manifest = [
        {
            "file_id": "synthetic-file-id",
            "file_name": "synthetic.maf.gz",
            "cases": [
                {
                    "case_id": "case-uuid",
                    "submitter_id": "TCGA-XX-1",
                    "samples": [
                        {
                            "sample_id": "tumor-sample-id",
                            "portions": [
                                {"analytes": [{"aliquots": [{"aliquot_id": "tumor-aliquot-uuid"}]}]}
                            ],
                        },
                        {
                            "sample_id": "normal-sample-id",
                            "portions": [
                                {
                                    "analytes": [
                                        {"aliquots": [{"aliquot_id": "normal-aliquot-uuid"}]}
                                    ]
                                }
                            ],
                        },
                    ],
                }
            ],
        }
    ]
    (mut_dir / "manifest.json").write_text(json.dumps(manifest))
    return project_dir


def test_load_for_project_parses_mini_maf(tmp_path: Path) -> None:
    project_dir = _build_mini_maf_tarball(tmp_path)
    by_case = mutations.load_for_project(project_dir)
    assert set(by_case) == {"case-uuid"}
    variants = by_case["case-uuid"]
    assert len(variants) == 2

    # source_file_id is set at parse time; sample_id FKs are resolved later in attach
    for v in variants:
        assert v["source_file_id"] == "synthetic-file-id"
        assert v["tumor_sample_id"] is None
        assert v["matched_normal_sample_id"] is None
        # MAF case_id passes through
        assert v["case_id"] == "case-uuid"

    # Numeric coercion worked (string in MAF -> int/float in record)
    assert variants[0]["t_depth"] == 100
    assert variants[0]["t_alt_count"] == 50
    assert variants[0]["gnomAD_AF"] == pytest.approx(0.001)
    assert variants[0]["Start_Position"] in (100, 200)


def _patient_samples_with_aliquots() -> list[dict]:
    """A patient `samples` slot rich enough to resolve the synthetic MAF's FKs.

    Mirrors the full GDC tree: samples → portions → analytes → aliquots.
    """
    return [
        {
            "sample_id": "tumor-sample-id",
            "portions": [{"analytes": [{"aliquots": [{"aliquot_id": "tumor-aliquot-uuid"}]}]}],
        },
        {
            "sample_id": "normal-sample-id",
            "portions": [{"analytes": [{"aliquots": [{"aliquot_id": "normal-aliquot-uuid"}]}]}],
        },
    ]


def test_attach_resolves_fks_from_patient_samples(tmp_path: Path) -> None:
    project_dir = _build_mini_maf_tarball(tmp_path)
    by_case = mutations.load_for_project(project_dir)

    rows = [
        {"case_id": "case-uuid", "samples": _patient_samples_with_aliquots()},
        {"case_id": "no-mutations-here", "samples": []},
    ]
    mutations.attach(rows, by_case)

    variants = rows[0]["samples_masked_somatic_mutation"]
    # FK resolution happened against the patient's samples list
    for v in variants:
        assert v["tumor_sample_id"] == "tumor-sample-id"
        assert v["matched_normal_sample_id"] == "normal-sample-id"
    # Sort: (chr1, 100) before (chr2, 200) regardless of MAF row order
    assert [v["Chromosome"] for v in variants] == ["chr1", "chr2"]
    # Empty mutations on the second patient
    assert rows[1]["samples_masked_somatic_mutation"] == []


def test_attach_leaves_fks_null_when_aliquot_not_in_samples(tmp_path: Path) -> None:
    """If the MAF references an aliquot the patient row doesn't carry (e.g. data
    drift), the FK resolves to None rather than crashing."""
    project_dir = _build_mini_maf_tarball(tmp_path)
    by_case = mutations.load_for_project(project_dir)

    rows = [{"case_id": "case-uuid", "samples": []}]  # no aliquots known
    mutations.attach(rows, by_case)
    for v in rows[0]["samples_masked_somatic_mutation"]:
        assert v["tumor_sample_id"] is None
        assert v["matched_normal_sample_id"] is None


def test_load_for_project_missing_dir_returns_empty(tmp_path: Path) -> None:
    assert mutations.load_for_project(tmp_path / "nonexistent") == {}


def test_mutation_fields_includes_all_140_maf_columns_plus_3_fk() -> None:
    names = [f.name for f in MUTATION_FIELDS]
    assert len(names) == 143
    fk = {"tumor_sample_id", "matched_normal_sample_id", "source_file_id"}
    assert set(names) - set(_MAF_COLUMNS) == fk
