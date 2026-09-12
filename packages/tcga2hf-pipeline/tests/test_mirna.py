"""Tests for the miRNA-Seq quantification table.

miRNA files carry no aliquot column of their own, so the FK is the file's
single associated aliquot entity — the one place these differ from the seg
files. The other thing worth pinning is the `cross-mapped` header: its
hyphen makes it both an invalid SQL identifier and unreachable as an
attribute on a pandas row, so it is renamed on the way in.
"""

from __future__ import annotations

import json
from pathlib import Path

import pyarrow as pa
from tcga2hf.schema import TABULAR_TABLES
from tcga2hf_pipeline import tabular

_CASE_UUID = "5fb3affa-3661-48e8-9d0f-8f0f7f6b0f11"
_SAMPLE_UUID = "e235e45d-9cbe-4c11-9a3a-1cf1d1e0aa10"
_ALIQUOT_UUID = "aff753ac-682c-449a-9a1e-3b2c1d0e5f22"

_MIRNA_TSV = (
    "miRNA_ID\tread_count\treads_per_million_miRNA_mapped\tcross-mapped\n"
    "hsa-let-7a-1\t32265\t10743.304123\tN\n"
    "hsa-let-7a-2\t31928\t10631.092950\tY\n"
    "hsa-mir-1302-2\t0\t0.000000\tN\n"
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
                        "portion_id": "portion-1",
                        "submitter_id": "TCGA-W5-AA33-01A-11",
                        "analytes": [
                            {
                                "aliquots": [
                                    {
                                        "aliquot_id": _ALIQUOT_UUID,
                                        "submitter_id": "TCGA-W5-AA33-01A-11R-A41D-13",
                                    }
                                ]
                            }
                        ],
                    }
                ],
            }
        ],
    }


def _build_project(
    tmp_path: Path,
    *,
    entities: list[dict] | None = None,
    file_name: str = "64cf1660.mirbase21.mirnas.quantification.txt",
    status: str = "downloaded",
    write_file: bool = True,
) -> Path:
    project_dir = tmp_path / "TCGA-XYZ"
    mod_dir = project_dir / "mirna"
    mod_dir.mkdir(parents=True)
    if write_file:
        (mod_dir / file_name).write_text(_MIRNA_TSV)
    entry = {
        "file_id": "synthetic-mirna-file-id",
        "file_name": file_name,
        "workflow_type": "BCGSC miRNA Profiling",
        "experimental_strategy": "miRNA-Seq",
        "cases": [{"case_id": _CASE_UUID, "submitter_id": "TCGA-W5-AA33"}],
        "associated_entities": (
            [
                {
                    "entity_id": _ALIQUOT_UUID,
                    "entity_type": "aliquot",
                    "entity_submitter_id": "TCGA-W5-AA33-01A-11R-A41D-13",
                }
            ]
            if entities is None
            else entities
        ),
        "_status": status,
    }
    (mod_dir / "manifest.json").write_text(json.dumps([entry]))
    return project_dir


def test_emits_one_row_per_mirna_with_resolved_fks(tmp_path: Path) -> None:
    rows = tabular._mirna_expression_quantification_rows([_case()], _build_project(tmp_path))

    assert len(rows) == 3
    assert [r["mirna_id"] for r in rows] == ["hsa-let-7a-1", "hsa-let-7a-2", "hsa-mir-1302-2"]
    for row in rows:
        assert row["case_submitter_id"] == "TCGA-W5-AA33"
        assert row["sample_id"] == _SAMPLE_UUID
        assert row["sample_type"] == "Primary Tumor"
        assert row["aliquot_id"] == _ALIQUOT_UUID
        assert row["aliquot_submitter_id"] == "TCGA-W5-AA33-01A-11R-A41D-13"


def test_cross_mapped_column_is_renamed_and_values_preserved(tmp_path: Path) -> None:
    rows = tabular._mirna_expression_quantification_rows([_case()], _build_project(tmp_path))
    assert [r["cross_mapped"] for r in rows] == ["N", "Y", "N"]


def test_counts_and_rpm_keep_their_types(tmp_path: Path) -> None:
    rows = tabular._mirna_expression_quantification_rows([_case()], _build_project(tmp_path))

    assert rows[0]["read_count"] == 32265
    assert isinstance(rows[0]["read_count"], int)
    assert rows[0]["reads_per_million_mirna_mapped"] == 10743.304123
    # A genuinely zero-count miRNA is data, not absence — it must survive.
    assert rows[2]["read_count"] == 0
    assert rows[2]["reads_per_million_mirna_mapped"] == 0.0


def test_ambiguous_aliquot_entities_are_skipped(tmp_path: Path) -> None:
    """Without exactly one aliquot there is no defensible FK."""
    two = [
        {"entity_id": _ALIQUOT_UUID, "entity_type": "aliquot", "entity_submitter_id": "a"},
        {"entity_id": "other-aliquot", "entity_type": "aliquot", "entity_submitter_id": "b"},
    ]
    assert (
        tabular._mirna_expression_quantification_rows(
            [_case()], _build_project(tmp_path, entities=two)
        )
        == []
    )
    assert (
        tabular._mirna_expression_quantification_rows(
            [_case()], _build_project(tmp_path / "none", entities=[])
        )
        == []
    )


def test_manifest_only_entries_are_skipped(tmp_path: Path) -> None:
    project = _build_project(tmp_path, status="manifest_only", write_file=False)
    assert tabular._mirna_expression_quantification_rows([_case()], project) == []


def test_missing_modality_dir_yields_no_rows(tmp_path: Path) -> None:
    project = tmp_path / "TCGA-EMPTY"
    project.mkdir()
    assert tabular._mirna_expression_quantification_rows([_case()], project) == []


def test_rows_satisfy_the_tabular_schema(tmp_path: Path) -> None:
    rows = tabular._mirna_expression_quantification_rows([_case()], _build_project(tmp_path))
    schema = TABULAR_TABLES["mirna_expression_quantification"]
    built = pa.Table.from_pylist(rows, schema=schema)
    assert built.num_rows == 3
    assert built.schema.equals(schema)
