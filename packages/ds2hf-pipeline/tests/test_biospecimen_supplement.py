"""Tests for the Biospecimen Supplement biotab parser.

Two things carry real risk here. First, form classification: the filenames
nest (`biospecimen_slide` is not a substring of `biospecimen_diagnostic_slides`
only by luck of GDC naming, and two different submitters ship these files).
Second, patient-barcode recovery: the specimen-level forms are keyed on their
own entity and several omit `bcr_patient_barcode` entirely, so the patient FK
has to come out of whichever barcode column the form does carry.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from ds2hf_pipeline.tcga import biospecimen_supplement as bs


def _biotab(header: list[str], rows: list[list[str]]) -> str:
    """A BCR biotab: names, aliased names, CDE tags, then data."""
    lines = [
        "\t".join(header),
        "\t".join(f"CDE_{h}" for h in header),
        "\t".join(f"CDE_ID:{i}" for i in range(len(header))),
    ]
    lines += ["\t".join(r) for r in rows]
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# Form classification
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("file_name", "expected"),
    [
        ("nationwidechildrens.org_biospecimen_sample_chol.txt", "biospecimen_sample"),
        ("nationwidechildrens.org_biospecimen_slide_chol.txt", "biospecimen_slide"),
        # The nesting cases: each must reach its own form, not the shorter one.
        (
            "nationwidechildrens.org_biospecimen_diagnostic_slides_chol.txt",
            "biospecimen_diagnostic_slides",
        ),
        (
            "nationwidechildrens.org_biospecimen_shipment_portion_chol.txt",
            "biospecimen_shipment_portion",
        ),
        ("nationwidechildrens.org_biospecimen_portion_chol.txt", "biospecimen_portion"),
        ("nationwidechildrens.org_ssf_tumor_samples_chol.txt", "ssf_tumor_samples"),
        ("nationwidechildrens.org_auxiliary_chol.txt", "auxiliary"),
        # A second submitter ships TCGA-LUAD; classification must not depend
        # on the `nationwidechildrens.org` prefix.
        ("genome.wustl.edu_biospecimen_cqcf_luad.txt", "biospecimen_cqcf"),
        ("genome.wustl.edu_biospecimen_sample_luad.txt", "biospecimen_sample"),
        # Clinical biotabs live in a sibling directory and are not ours.
        ("nationwidechildrens.org_clinical_patient_chol.txt", None),
        ("some_unrelated_file.txt", None),
    ],
)
def test_form_kind_classification(file_name: str, expected: str | None) -> None:
    assert bs._form_kind(file_name) == expected


def test_every_wanted_form_maps_to_a_declared_tabular_suffix() -> None:
    """A form we download but can't name a table for would be silently dropped."""
    suffixes = {bs._suffix_for(kind) for kind in bs.WANTED_FORMS}
    assert suffixes == set(bs.TABULAR_FORM_KINDS)


# ---------------------------------------------------------------------------
# Patient barcode recovery
# ---------------------------------------------------------------------------


def test_patient_barcode_recovered_from_entity_barcodes() -> None:
    """Forms without `bcr_patient_barcode` still resolve a patient FK."""
    assert (
        bs._case_submitter_id({"bcr_aliquot_barcode": "TCGA-3X-AAV9-01A-11D-A42S-01"})
        == "TCGA-3X-AAV9"
    )
    assert bs._case_submitter_id({"bcr_sample_barcode": "TCGA-W5-AA33-10A"}) == "TCGA-W5-AA33"
    assert bs._case_submitter_id({"bcr_patient_barcode": "TCGA-W5-AA33"}) == "TCGA-W5-AA33"


def test_patient_barcode_column_wins_over_entity_columns() -> None:
    rec = {
        "bcr_patient_barcode": "TCGA-AA-1111",
        "bcr_sample_barcode": "TCGA-BB-2222-01A",
    }
    assert bs._case_submitter_id(rec) == "TCGA-AA-1111"


@pytest.mark.parametrize(
    "rec",
    [
        {},
        {"bcr_patient_barcode": ""},
        {"bcr_patient_barcode": "[Not Available]"},
        # BCR leaves literal header echoes in some forms.
        {"bcr_patient_barcode": "bcr_patient_barcode"},
        {"bcr_sample_barcode": "TCGA-AA"},
    ],
)
def test_unresolvable_barcodes_return_none(rec: dict) -> None:
    assert bs._case_submitter_id(rec) is None


# ---------------------------------------------------------------------------
# Tabular emission
# ---------------------------------------------------------------------------


def test_build_tabular_rows_routes_forms_and_prefixes_the_fk(tmp_path: Path) -> None:
    supp = tmp_path / "biospecimen_supplement"
    supp.mkdir()
    (supp / "nationwidechildrens.org_biospecimen_slide_chol.txt").write_text(
        _biotab(
            ["bcr_patient_uuid", "bcr_slide_barcode", "percent_tumor_nuclei", "percent_necrosis"],
            [
                ["uuid-1", "TCGA-W5-AA33-01A-01-TS1", "80", "5"],
                ["uuid-2", "TCGA-W5-AA34-01A-01-TS1", "60", "10"],
            ],
        )
    )
    (supp / "nationwidechildrens.org_ssf_tumor_samples_chol.txt").write_text(
        _biotab(
            ["bcr_patient_barcode", "tumor_nuclei_requirements_indicator"],
            [["TCGA-W5-AA33", "YES"]],
        )
    )

    rows = bs.build_tabular_rows(supp)

    assert [r["case_submitter_id"] for r in rows["slide"]] == ["TCGA-W5-AA33", "TCGA-W5-AA34"]
    assert rows["slide"][0]["percent_tumor_nuclei"] == "80"
    assert len(rows["ssf_tumor_samples"]) == 1
    # Forms with no file for this project are present but empty, so the
    # caller sees a stable key set.
    assert rows["cqcf"] == []
    assert rows["auxiliary"] == []


def test_rows_without_a_resolvable_patient_are_dropped(tmp_path: Path) -> None:
    supp = tmp_path / "biospecimen_supplement"
    supp.mkdir()
    (supp / "nationwidechildrens.org_biospecimen_sample_chol.txt").write_text(
        _biotab(
            ["bcr_patient_barcode", "sample_type"],
            [["TCGA-W5-AA33", "Primary Tumor"], ["", "Primary Tumor"]],
        )
    )
    rows = bs.build_tabular_rows(supp)
    assert len(rows["sample"]) == 1
    assert rows["sample"][0]["case_submitter_id"] == "TCGA-W5-AA33"


def test_two_submitters_of_one_form_are_concatenated(tmp_path: Path) -> None:
    """TCGA-LUAD ships both nationwidechildrens.org and genome.wustl.edu forms."""
    supp = tmp_path / "biospecimen_supplement"
    supp.mkdir()
    (supp / "nationwidechildrens.org_biospecimen_sample_luad.txt").write_text(
        _biotab(["bcr_patient_barcode", "shared_col"], [["TCGA-AA-1111", "a"]])
    )
    (supp / "genome.wustl.edu_biospecimen_sample_luad.txt").write_text(
        _biotab(["bcr_patient_barcode", "wustl_only_col"], [["TCGA-BB-2222", "b"]])
    )

    rows = bs.build_tabular_rows(supp)
    assert len(rows["sample"]) == 2
    by_patient = {r["case_submitter_id"]: r for r in rows["sample"]}
    assert set(by_patient) == {"TCGA-AA-1111", "TCGA-BB-2222"}
    # Each row keeps only its own submitter's columns; parquet infers the
    # union and pads the rest with null.
    assert "wustl_only_col" not in by_patient["TCGA-AA-1111"]
    assert "shared_col" not in by_patient["TCGA-BB-2222"]


def test_missing_dir_returns_stable_empty_shape(tmp_path: Path) -> None:
    rows = bs.build_tabular_rows(tmp_path / "does-not-exist")
    assert set(rows) == set(bs.TABULAR_FORM_KINDS)
    assert all(v == [] for v in rows.values())
