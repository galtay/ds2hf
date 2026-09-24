"""Tests for the consolidated-row loaders of the four newest modalities.

These cover what the tabular emitters' tests don't: the struct-of-arrays
shape, and FK resolution against a *built patient row* rather than a raw GDC
case dict. The two layouts resolve sample ids from different structures, so
each needs its own coverage even though they read the same source files.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pyarrow as pa
from ds2hf.tcga.schema import PATIENTS
from ds2hf_pipeline.tcga import biospecimen_supplement, copy_number, mirna, protein_expression

_CASE_UUID = "5fb3affa-3661-48e8-9d0f-8f0f7f6b0f11"
_SAMPLE_UUID = "e235e45d-9cbe-4c11-9a3a-1cf1d1e0aa10"
_PORTION_UUID = "e09e7127-84b5-4b17-9c3d-2e1f0a9b8c7d"
_TUMOR_ALIQUOT = "b3c68473-78a7-44d2-95bf-7997a7e6c0e8"
_NORMAL_ALIQUOT = "451f8cc6-05fe-48f6-9f2e-1b8c0d3a77aa"
_MIRNA_ALIQUOT = "aff753ac-682c-449a-9a1e-3b2c1d0e5f22"


def _patient_row() -> dict[str, Any]:
    """A built patient row — samples nested the way `to_patient_rows` leaves them."""
    return {
        "case_id": _CASE_UUID,
        "case_submitter_id": "TCGA-W5-AA33",
        "samples": [
            {
                "sample_id": _SAMPLE_UUID,
                "submitter_id": "TCGA-W5-AA33-01A",
                "sample_type": "Primary Tumor",
                "portions": [
                    {
                        "portion_id": _PORTION_UUID,
                        "submitter_id": "TCGA-W5-AA33-01A-21",
                        "analytes": [
                            {
                                "aliquots": [
                                    {"aliquot_id": _TUMOR_ALIQUOT},
                                    {"aliquot_id": _MIRNA_ALIQUOT},
                                ]
                            }
                        ],
                    }
                ],
            }
        ],
    }


def _write_modality(
    tmp_path: Path,
    modality: str,
    file_name: str,
    content: str,
    entities: list[dict],
    **extra: Any,
) -> Path:
    project_dir = tmp_path / "TCGA-XYZ"
    mod_dir = project_dir / modality
    mod_dir.mkdir(parents=True, exist_ok=True)
    (mod_dir / file_name).write_text(content)
    entry = {
        "file_id": f"file-{modality}",
        "file_name": file_name,
        "cases": [{"case_id": _CASE_UUID, "submitter_id": "TCGA-W5-AA33"}],
        "associated_entities": entities,
        "_status": "downloaded",
        **extra,
    }
    (mod_dir / "manifest.json").write_text(json.dumps([entry]))
    return project_dir


def _aliquot_entity(entity_id: str) -> dict:
    return {"entity_id": entity_id, "entity_type": "aliquot", "entity_submitter_id": "x"}


# ---------------------------------------------------------------------------
# Copy number
# ---------------------------------------------------------------------------

_ASCN_TSV = (
    "GDC_Aliquot\tChromosome\tStart\tEnd\tCopy_Number\tMajor_Copy_Number\tMinor_Copy_Number\n"
    f"{_TUMOR_ALIQUOT}\tchr1\t62920\t33086177\t1\t1\t0\n"
    f"{_TUMOR_ALIQUOT}\tchr2\t100\t200\t3\t2\t1\n"
)
_MASKED_TSV = (
    "GDC_Aliquot\tChromosome\tStart\tEnd\tNum_Probes\tSegment_Mean\n"
    f"{_TUMOR_ALIQUOT}\t1\t3301765\t6028611\t1917\t-0.2829\n"
    f"{_TUMOR_ALIQUOT}\tX\t100\t200\t2\t-2.0991\n"
)


def test_allele_specific_is_struct_of_arrays_not_row_per_segment(tmp_path: Path) -> None:
    """One record per (aliquot, workflow), with segments as index-aligned arrays."""
    project = _write_modality(
        tmp_path,
        "copy_number_allele_specific",
        "a.ascat3.allelic_specific.seg.txt",
        _ASCN_TSV,
        [_aliquot_entity(_NORMAL_ALIQUOT), _aliquot_entity(_TUMOR_ALIQUOT)],
        workflow_type="ASCAT3",
        experimental_strategy="Genotyping Array",
    )
    by_case = copy_number.load_allele_specific_for_project(project)
    rows = copy_number.attach(
        [_patient_row()], by_case, "samples_allele_specific_copy_number_segment"
    )
    records = rows[0]["samples_allele_specific_copy_number_segment"]

    assert len(records) == 1
    rec = records[0]
    assert rec["aliquot_id"] == _TUMOR_ALIQUOT
    assert rec["matched_normal_aliquot_id"] == _NORMAL_ALIQUOT
    assert rec["workflow_type"] == "ASCAT3"
    assert rec["sample_id"] == _SAMPLE_UUID
    assert rec["chromosome"] == ["chr1", "chr2"]
    assert rec["copy_number"] == [1, 3]
    assert rec["major_copy_number"] == [1, 2]
    assert rec["minor_copy_number"] == [0, 1]
    # Every array is index-aligned to `chromosome`.
    n = len(rec["chromosome"])
    for key in ("start", "end", "copy_number", "major_copy_number", "minor_copy_number"):
        assert len(rec[key]) == n


def test_masked_keeps_bare_chromosome_names_as_strings(tmp_path: Path) -> None:
    project = _write_modality(
        tmp_path,
        "copy_number_masked",
        "m.nocnv_grch38.seg.v2.txt",
        _MASKED_TSV,
        [_aliquot_entity(_TUMOR_ALIQUOT)],
        workflow_type="DNAcopy",
    )
    by_case = copy_number.load_masked_for_project(project)
    rows = copy_number.attach([_patient_row()], by_case, "samples_masked_copy_number_segment")
    rec = rows[0]["samples_masked_copy_number_segment"][0]

    assert rec["chromosome"] == ["1", "X"]
    assert all(isinstance(c, str) for c in rec["chromosome"])
    assert rec["num_probes"] == [1917, 2]
    assert rec["segment_mean"] == [-0.2829, -2.0991]


def test_several_callers_for_one_aliquot_stay_separate_records(tmp_path: Path) -> None:
    """ASCAT2 and ASCAT3 disagree, so they must not be merged into one record."""
    project = tmp_path / "TCGA-XYZ"
    mod = project / "copy_number_allele_specific"
    mod.mkdir(parents=True)
    entries = []
    for wf, cn in (("ASCAT2", 2), ("ASCAT3", 4)):
        name = f"{wf}.seg.txt"
        (mod / name).write_text(
            "GDC_Aliquot\tChromosome\tStart\tEnd\tCopy_Number\t"
            "Major_Copy_Number\tMinor_Copy_Number\n"
            f"{_TUMOR_ALIQUOT}\tchr1\t1\t100\t{cn}\t{cn}\t0\n"
        )
        entries.append(
            {
                "file_id": f"file-{wf}",
                "file_name": name,
                "workflow_type": wf,
                "cases": [{"case_id": _CASE_UUID, "submitter_id": "TCGA-W5-AA33"}],
                "associated_entities": [_aliquot_entity(_TUMOR_ALIQUOT)],
                "_status": "downloaded",
            }
        )
    (mod / "manifest.json").write_text(json.dumps(entries))

    by_case = copy_number.load_allele_specific_for_project(project)
    rows = copy_number.attach(
        [_patient_row()], by_case, "samples_allele_specific_copy_number_segment"
    )
    records = rows[0]["samples_allele_specific_copy_number_segment"]

    assert len(records) == 2
    # Sorted by (aliquot_id, workflow_type) so the order is deterministic.
    assert [r["workflow_type"] for r in records] == ["ASCAT2", "ASCAT3"]
    assert [r["copy_number"][0] for r in records] == [2, 4]


def test_cases_without_copy_number_get_empty_lists(tmp_path: Path) -> None:
    project = tmp_path / "TCGA-EMPTY"
    project.mkdir()
    assert copy_number.load_allele_specific_for_project(project) == {}
    rows = copy_number.attach([_patient_row()], {}, "samples_masked_copy_number_segment")
    assert rows[0]["samples_masked_copy_number_segment"] == []


# ---------------------------------------------------------------------------
# miRNA
# ---------------------------------------------------------------------------

_MIRNA_TSV = (
    "miRNA_ID\tread_count\treads_per_million_miRNA_mapped\tcross-mapped\n"
    "hsa-let-7a-1\t32265\t10743.304123\tN\n"
    "hsa-let-7a-2\t31928\t10631.092950\tY\n"
)


def test_mirna_record_shape_and_fk(tmp_path: Path) -> None:
    project = _write_modality(
        tmp_path,
        "mirna",
        "x.mirbase21.mirnas.quantification.txt",
        _MIRNA_TSV,
        [_aliquot_entity(_MIRNA_ALIQUOT)],
    )
    rows = mirna.attach([_patient_row()], mirna.load_for_project(project))
    rec = rows[0]["samples_mirna_expression_quantification"][0]

    assert rec["sample_id"] == _SAMPLE_UUID
    assert rec["aliquot_id"] == _MIRNA_ALIQUOT
    assert rec["mirna_id"] == ["hsa-let-7a-1", "hsa-let-7a-2"]
    assert rec["read_count"] == [32265, 31928]
    # Source header `cross-mapped` is renamed; values untouched.
    assert rec["cross_mapped"] == ["N", "Y"]


# ---------------------------------------------------------------------------
# RPPA
# ---------------------------------------------------------------------------

_RPPA_TSV = (
    "AGID\tlab_id\tcatalog_number\tset_id\tpeptide_target\tprotein_expression\n"
    "AGID00100\t882\tsc-628\tOld\t1433BETA\t0.057811\n"
    "AGID00002\t3\t9456\tV1.2\t4EBP1_pS65\tNA\n"
)


def test_rppa_resolves_sample_through_the_portion(tmp_path: Path) -> None:
    project = _write_modality(
        tmp_path,
        "protein_expression",
        "rppa.tsv",
        _RPPA_TSV,
        [
            {
                "entity_id": _PORTION_UUID,
                "entity_type": "portion",
                "entity_submitter_id": "TCGA-W5-AA33-01A-21",
            }
        ],
    )
    rows = protein_expression.attach([_patient_row()], protein_expression.load_for_project(project))
    rec = rows[0]["samples_protein_expression_quantification"][0]

    assert rec["portion_id"] == _PORTION_UUID
    assert rec["sample_id"] == _SAMPLE_UUID
    assert rec["peptide_target"] == ["1433BETA", "4EBP1_pS65"]
    # `NA` is a failed measurement -> null, and the entry is kept in place so
    # the arrays stay index-aligned.
    assert rec["protein_expression"] == [0.057811, None]
    # Identifier columns stay strings.
    assert rec["lab_id"] == ["882", "3"]


# ---------------------------------------------------------------------------
# Biospecimen supplements
# ---------------------------------------------------------------------------


def _biotab(header: list[str], rows: list[list[str]]) -> str:
    lines = [
        "\t".join(header),
        "\t".join(f"CDE_{h}" for h in header),
        "\t".join(f"CDE_ID:{i}" for i in range(len(header))),
    ]
    return "\n".join(lines + ["\t".join(r) for r in rows]) + "\n"


def test_biospecimen_supplements_group_by_patient_barcode(tmp_path: Path) -> None:
    supp = tmp_path / "biospecimen_supplement"
    supp.mkdir()
    (supp / "nationwidechildrens.org_biospecimen_slide_chol.txt").write_text(
        _biotab(
            ["bcr_slide_barcode", "percent_tumor_nuclei"],
            [["TCGA-W5-AA33-01A-01-TS1", "80"], ["TCGA-W5-AA33-01A-02-TS2", "70"]],
        )
    )
    (supp / "nationwidechildrens.org_biospecimen_sample_chol.txt").write_text(
        _biotab(["bcr_patient_barcode", "composition"], [["TCGA-W5-AA33", "Solid Tissue"]])
    )

    by_case = biospecimen_supplement.load_supplements_for_project(supp)
    rows = biospecimen_supplement.attach_supplements([_patient_row()], by_case)
    supp_col = rows[0]["biospecimen_supplement"]

    assert len(supp_col["slide"]) == 2
    assert supp_col["slide"][0]["percent_tumor_nuclei"] == "80"
    assert len(supp_col["sample"]) == 1
    # Every declared form is present as a key, empty where the project has none.
    assert set(supp_col) == set(biospecimen_supplement.TABULAR_FORM_KINDS)
    assert supp_col["cqcf"] == []


def test_patient_without_biospecimen_data_gets_none(tmp_path: Path) -> None:
    rows = biospecimen_supplement.attach_supplements([_patient_row()], {})
    assert rows[0]["biospecimen_supplement"] is None


# ---------------------------------------------------------------------------
# Schema conformance
# ---------------------------------------------------------------------------


def test_all_four_modalities_satisfy_the_patients_schema(tmp_path: Path) -> None:
    """A row carrying all four new columns must load under the strict schema."""
    row = _patient_row()
    ascn = _write_modality(
        tmp_path / "a",
        "copy_number_allele_specific",
        "a.seg.txt",
        _ASCN_TSV,
        [_aliquot_entity(_TUMOR_ALIQUOT)],
        workflow_type="ASCAT3",
    )
    masked = _write_modality(
        tmp_path / "m",
        "copy_number_masked",
        "m.seg.txt",
        _MASKED_TSV,
        [_aliquot_entity(_TUMOR_ALIQUOT)],
        workflow_type="DNAcopy",
    )
    mi = _write_modality(
        tmp_path / "n", "mirna", "n.txt", _MIRNA_TSV, [_aliquot_entity(_MIRNA_ALIQUOT)]
    )
    rp = _write_modality(
        tmp_path / "p",
        "protein_expression",
        "p.tsv",
        _RPPA_TSV,
        [{"entity_id": _PORTION_UUID, "entity_type": "portion", "entity_submitter_id": "x"}],
    )
    copy_number.attach(
        [row],
        copy_number.load_allele_specific_for_project(ascn),
        "samples_allele_specific_copy_number_segment",
    )
    copy_number.attach(
        [row], copy_number.load_masked_for_project(masked), "samples_masked_copy_number_segment"
    )
    mirna.attach([row], mirna.load_for_project(mi))
    protein_expression.attach([row], protein_expression.load_for_project(rp))

    # Fill the columns the schema requires that this fixture doesn't exercise.
    for col in PATIENTS.names:
        row.setdefault(col, [] if str(PATIENTS.field(col).type).startswith("list") else None)

    table = pa.Table.from_pylist([row], schema=PATIENTS)
    assert table.num_rows == 1
    assert table.schema.equals(PATIENTS)
    assert len(table["samples_allele_specific_copy_number_segment"][0]) == 1
    assert len(table["samples_mirna_expression_quantification"][0]) == 1
    assert len(table["samples_protein_expression_quantification"][0]) == 1


def test_every_molecular_column_is_initialized_by_to_patient_rows() -> None:
    """A project with zero files for a modality must still get the column.

    `attach` is only called when a modality has data, so `to_patient_rows`
    is the sole thing putting an empty column on the row otherwise. Missing
    an entry there means the build crashes on the first project lacking that
    modality — TCGA-LAML has no RPPA and no pathology reports, and is 14th
    alphabetically, so the failure lands 40 minutes into a full rebuild.
    """
    from ds2hf_pipeline.tcga.clinical import to_patient_rows

    rows = to_patient_rows([{"case_id": "c", "submitter_id": "TCGA-XX-0000"}])
    molecular = {n for n in PATIENTS.names if n.startswith("samples_")}
    missing = molecular - set(rows[0])
    assert missing == set(), f"columns never initialized: {sorted(missing)}"
    assert all(rows[0][n] == [] for n in molecular)


def test_cases_table_excludes_every_molecular_column() -> None:
    """A new molecular column must not leak into the tabular `cases` table."""
    from ds2hf.tcga.schema import TABULAR_CASES_FIELDS

    names = {f.name for f in TABULAR_CASES_FIELDS}
    molecular = {n for n in PATIENTS.names if n.startswith("samples_")}
    assert names & molecular == set(), f"leaked into cases: {sorted(names & molecular)}"
