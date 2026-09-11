"""Tests for the two modalities added to reach full open-access coverage.

Methylation betas are the only source file in the dataset with **no header
row**, and the only table whose `platform` genuinely changes what a value
means. miRNA isoforms are the per-isoform companion to the mature-miRNA
table and share its `cross-mapped` column-name problem.
"""

from __future__ import annotations

import json
from pathlib import Path

import pyarrow as pa
from tcga2hf.schema import TABULAR_TABLES
from tcga2hf_pipeline import tabular

_CASE_UUID = "5fb3affa-3661-48e8-9d0f-8f0f7f6b0f11"
_SAMPLE_UUID = "e235e45d-9cbe-4c11-9a3a-1cf1d1e0aa10"
_ALIQUOT = "b3c68473-78a7-44d2-95bf-7997a7e6c0e8"

# No header line — the real files start straight at the first probe. `NA`
# is what SeSAMe writes for a masked probe.
_BETA_TXT = (
    "cg00000029\t0.176824809353763\n"
    "cg00000108\t0.94695735792137\n"
    "ch.2.30415474F\tNA\n"
    "rs9363764\t0.5121\n"
)

_ISOFORM_TSV = (
    "miRNA_ID\tisoform_coords\tread_count\treads_per_million_miRNA_mapped\t"
    "cross-mapped\tmiRNA_region\n"
    "hsa-let-7a-1\thg38:chr9:94175942-94175961:+\t1\t0.246769\tN\tprecursor\n"
    "hsa-let-7a-1\thg38:chr9:94175961-94175982:+\t4\t0.987077\tY\tmature,MIMAT0000062\n"
)


def _case() -> dict:
    return {
        "case_id": _CASE_UUID,
        "submitter_id": "TCGA-W5-AA33",
        "samples": [
            {
                "sample_id": _SAMPLE_UUID,
                "submitter_id": "TCGA-W5-AA33-01A",
                "sample_type": "Primary Tumor",
                "portions": [
                    {
                        "analytes": [
                            {
                                "aliquots": [
                                    {
                                        "aliquot_id": _ALIQUOT,
                                        "submitter_id": "TCGA-W5-AA33-01A-11D-A416-05",
                                    }
                                ]
                            }
                        ]
                    }
                ],
            }
        ],
    }


def _project(tmp_path: Path, modality: str, content: str, **extra) -> Path:
    project_dir = tmp_path / "TCGA-XYZ"
    mod_dir = project_dir / modality
    mod_dir.mkdir(parents=True)
    file_name = f"synthetic.{modality}.txt"
    (mod_dir / file_name).write_text(content)
    entry = {
        "file_id": f"synthetic-{modality}-file-id",
        "file_name": file_name,
        "cases": [{"case_id": _CASE_UUID, "submitter_id": "TCGA-W5-AA33"}],
        "associated_entities": [{"entity_id": _ALIQUOT, "entity_type": "aliquot"}],
        "_status": "downloaded",
        **extra,
    }
    (mod_dir / "manifest.json").write_text(json.dumps([entry]))
    return project_dir


def test_methylation_reads_a_headerless_two_column_file(tmp_path: Path) -> None:
    project = _project(
        tmp_path, "methylation", _BETA_TXT, platform="Illumina Human Methylation 450"
    )
    batches = list(tabular._methylation_beta_value_batches([_case()], project))
    rows = [r for b in batches for r in b.to_pylist()]
    assert len(rows) == 4
    # The first line is data, not a header: losing it would silently drop a probe.
    assert rows[0]["probe_id"] == "cg00000029"
    assert rows[0]["beta_value"] == 0.176824809353763
    assert rows[0]["aliquot_id"] == _ALIQUOT
    assert rows[0]["sample_id"] == _SAMPLE_UUID
    assert rows[0]["sample_type"] == "Primary Tumor"
    assert rows[0]["platform"] == "Illumina Human Methylation 450"


def test_masked_probes_are_null_not_zero(tmp_path: Path) -> None:
    """A masked probe means "not trustworthy"; 0.0 would mean "unmethylated"."""
    project = _project(tmp_path, "methylation", _BETA_TXT, platform="x")
    rows = [
        r
        for b in tabular._methylation_beta_value_batches([_case()], project)
        for r in b.to_pylist()
    ]
    masked = next(r for r in rows if r["probe_id"] == "ch.2.30415474F")
    assert masked["beta_value"] is None


def test_methylation_batches_conform_to_the_published_schema(tmp_path: Path) -> None:
    project = _project(tmp_path, "methylation", _BETA_TXT, platform="x")
    batches = list(tabular._methylation_beta_value_batches([_case()], project))
    assert batches
    for b in batches:
        assert b.schema.equals(TABULAR_TABLES["methylation_beta_value"])


def test_isoform_rows_rename_the_hyphenated_column(tmp_path: Path) -> None:
    project = _project(tmp_path, "mirna_isoform", _ISOFORM_TSV)
    rows = tabular._isoform_expression_quantification_rows([_case()], project)
    assert len(rows) == 2
    assert rows[0]["mirna_id"] == "hsa-let-7a-1"
    assert rows[0]["isoform_coords"] == "hg38:chr9:94175942-94175961:+"
    assert (rows[0]["read_count"], rows[0]["reads_per_million_mirna_mapped"]) == (1, 0.246769)
    # Source header is `cross-mapped`; a hyphen is not a legal SQL identifier.
    assert rows[0]["cross_mapped"] == "N"
    assert rows[1]["cross_mapped"] == "Y"
    assert rows[1]["mirna_region"] == "mature,MIMAT0000062"
    assert rows[0]["aliquot_id"] == _ALIQUOT
    pa.Table.from_pylist(rows, schema=TABULAR_TABLES["isoform_expression_quantification"])


def test_both_modalities_skip_files_that_name_no_single_aliquot(tmp_path: Path) -> None:
    project = tmp_path / "TCGA-XYZ"
    for modality, content in (("methylation", _BETA_TXT), ("mirna_isoform", _ISOFORM_TSV)):
        mod_dir = project / modality
        mod_dir.mkdir(parents=True)
        (mod_dir / "f.txt").write_text(content)
        (mod_dir / "manifest.json").write_text(
            json.dumps(
                [
                    {
                        "file_id": "fid",
                        "file_name": "f.txt",
                        "cases": [{"case_id": _CASE_UUID}],
                        "associated_entities": [],  # no aliquot to attach to
                        "_status": "downloaded",
                    }
                ]
            )
        )
    assert list(tabular._methylation_beta_value_batches([_case()], project)) == []
    assert tabular._isoform_expression_quantification_rows([_case()], project) == []
