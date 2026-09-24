"""Tests for the RPPA protein expression table.

RPPA is the one molecular modality that attaches to a *portion* rather than
an aliquot, so the sample FK has to be resolved by walking samples ->
portions. The other thing worth pinning is that the source encodes a failed
measurement as the literal string `NA`, which must land as null rather than
as the string "NA" or as a zero.
"""

from __future__ import annotations

import json
from pathlib import Path

import pyarrow as pa
from ds2hf.tcga.schema import TABULAR_TABLES
from ds2hf_pipeline.tcga import tabular

_CASE_UUID = "f47dbd99-920b-4b39-8b0e-1a2b3c4d5e6f"
_SAMPLE_UUID = "9652df1d-041f-4e14-b2a1-0c9d8e7f6a5b"
_PORTION_UUID = "e09e7127-84b5-4b17-9c3d-2e1f0a9b8c7d"

_RPPA_TSV = (
    "AGID\tlab_id\tcatalog_number\tset_id\tpeptide_target\tprotein_expression\n"
    "AGID00100\t882\tsc-628\tOld\t1433BETA\t0.057811\n"
    "AGID00111\t913\tsc-23957\tOld\t1433EPSILON\t-0.020137\n"
    "AGID00002\t3\t9456\tV1.2\t4EBP1_pS65\tNA\n"
)


def _case() -> dict:
    return {
        "case_id": _CASE_UUID,
        "submitter_id": "TCGA-W5-AA2Q",
        "samples": [
            {
                "sample_id": _SAMPLE_UUID,
                "submitter_id": "TCGA-W5-AA2Q-01A",
                "sample_type": "Primary Tumor",
                "portions": [
                    {
                        "portion_id": _PORTION_UUID,
                        "submitter_id": "TCGA-W5-AA2Q-01A-21",
                        "analytes": [],
                    }
                ],
            }
        ],
    }


def _build_project(
    tmp_path: Path,
    *,
    entities: list[dict] | None = None,
    status: str = "downloaded",
    write_file: bool = True,
) -> Path:
    project_dir = tmp_path / "TCGA-XYZ"
    mod_dir = project_dir / "protein_expression"
    mod_dir.mkdir(parents=True)
    file_name = "TCGA-W5-AA2Q-01A-21-A45N-20_RPPA_data.tsv"
    if write_file:
        (mod_dir / file_name).write_text(_RPPA_TSV)
    entry = {
        "file_id": "synthetic-rppa-file-id",
        "file_name": file_name,
        "workflow_type": None,
        "experimental_strategy": "Reverse Phase Protein Array",
        "cases": [{"case_id": _CASE_UUID, "submitter_id": "TCGA-W5-AA2Q"}],
        "associated_entities": (
            [
                {
                    "entity_id": _PORTION_UUID,
                    "entity_type": "portion",
                    "entity_submitter_id": "TCGA-W5-AA2Q-01A-21",
                }
            ]
            if entities is None
            else entities
        ),
        "_status": status,
    }
    (mod_dir / "manifest.json").write_text(json.dumps([entry]))
    return project_dir


def test_sample_is_resolved_through_the_portion(tmp_path: Path) -> None:
    rows = tabular._protein_expression_quantification_rows([_case()], _build_project(tmp_path))

    assert len(rows) == 3
    for row in rows:
        assert row["portion_id"] == _PORTION_UUID
        assert row["portion_submitter_id"] == "TCGA-W5-AA2Q-01A-21"
        assert row["sample_id"] == _SAMPLE_UUID
        assert row["sample_submitter_id"] == "TCGA-W5-AA2Q-01A"
        assert row["sample_type"] == "Primary Tumor"


def test_na_measurements_become_null_not_zero(tmp_path: Path) -> None:
    """6.5% of TCGA-CHOL's RPPA cells are the literal string `NA`."""
    rows = tabular._protein_expression_quantification_rows([_case()], _build_project(tmp_path))

    assert rows[0]["protein_expression"] == 0.057811
    assert rows[1]["protein_expression"] == -0.020137
    assert rows[2]["protein_expression"] is None
    # The row itself survives — the antibody was on the panel.
    assert rows[2]["peptide_target"] == "4EBP1_pS65"
    assert rows[2]["set_id"] == "V1.2"


def test_lab_id_stays_a_string(tmp_path: Path) -> None:
    """It is an identifier, so it must not acquire a decimal point."""
    rows = tabular._protein_expression_quantification_rows([_case()], _build_project(tmp_path))
    assert [r["lab_id"] for r in rows] == ["882", "913", "3"]


def test_unresolvable_portion_still_emits_rows_with_null_sample(tmp_path: Path) -> None:
    """A portion GDC names but the case tree doesn't is a null FK, not a dropped row.

    The measurement is real data; losing it because a biospecimen id didn't
    round-trip would be worse than reporting an unresolved sample.
    """
    unknown = [
        {
            "entity_id": "portion-not-in-case-tree",
            "entity_type": "portion",
            "entity_submitter_id": "TCGA-W5-AA2Q-01A-99",
        }
    ]
    rows = tabular._protein_expression_quantification_rows(
        [_case()], _build_project(tmp_path, entities=unknown)
    )
    assert len(rows) == 3
    assert all(r["sample_id"] is None for r in rows)
    # The submitter id falls back to what GDC reported on the entity.
    assert all(r["portion_submitter_id"] == "TCGA-W5-AA2Q-01A-99" for r in rows)


def test_non_portion_entities_are_skipped(tmp_path: Path) -> None:
    aliquot = [{"entity_id": "some-aliquot", "entity_type": "aliquot", "entity_submitter_id": "x"}]
    assert (
        tabular._protein_expression_quantification_rows(
            [_case()], _build_project(tmp_path, entities=aliquot)
        )
        == []
    )


def test_manifest_only_entries_are_skipped(tmp_path: Path) -> None:
    project = _build_project(tmp_path, status="manifest_only", write_file=False)
    assert tabular._protein_expression_quantification_rows([_case()], project) == []


def test_missing_modality_dir_yields_no_rows(tmp_path: Path) -> None:
    project = tmp_path / "TCGA-EMPTY"
    project.mkdir()
    assert tabular._protein_expression_quantification_rows([_case()], project) == []


def test_rows_satisfy_the_tabular_schema(tmp_path: Path) -> None:
    rows = tabular._protein_expression_quantification_rows([_case()], _build_project(tmp_path))
    schema = TABULAR_TABLES["protein_expression_quantification"]
    built = pa.Table.from_pylist(rows, schema=schema)
    assert built.num_rows == 3
    assert built.schema.equals(schema)
