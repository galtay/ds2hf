from __future__ import annotations

import json
from pathlib import Path

import pyarrow.parquet as pq
import pytest
from tcga2hf.schema import (
    ALIQUOT_FIELDS,
    ANALYTE_FIELDS,
    DEMOGRAPHIC_FIELDS,
    DIAGNOSIS_FIELDS,
    EXPOSURE_FIELDS,
    FAMILY_HISTORY_FIELDS,
    FOLLOW_UP_FIELDS,
    PATIENTS,
    PORTION_FIELDS,
    SAMPLE_FIELDS,
    TREATMENT_FIELDS,
)
from tcga2hf_pipeline import clinical, dataset_card

FIXTURE = Path(__file__).parent / "fixtures" / "case_chol_one.json"


@pytest.fixture
def one_case() -> list[dict]:
    return [json.loads(FIXTURE.read_text())]


def test_to_patient_rows_one_row_per_case(one_case: list[dict]) -> None:
    rows = clinical.to_patient_rows(one_case)
    assert len(rows) == 1
    row = rows[0]
    assert row["case_id"] == one_case[0]["case_id"]
    assert row["case_submitter_id"].startswith("TCGA-")
    assert row["project_id"] == "TCGA-CHOL"
    assert row["gdc_portal_url"] == f"https://portal.gdc.cancer.gov/cases/{row['case_id']}"


def test_demographic_is_a_struct_with_expected_keys(one_case: list[dict]) -> None:
    row = clinical.to_patient_rows(one_case)[0]
    demo = row["demographic"]
    assert demo is not None
    assert set(demo.keys()) == {f.name for f in DEMOGRAPHIC_FIELDS}


def test_diagnoses_list_length_and_treatments_nested(one_case: list[dict]) -> None:
    row = clinical.to_patient_rows(one_case)[0]
    src_dx = one_case[0].get("diagnoses") or []
    assert len(row["diagnoses"]) == len(src_dx)

    expected_dx_keys = {f.name for f in DIAGNOSIS_FIELDS}
    expected_tx_keys = {f.name for f in TREATMENT_FIELDS}
    for src, out in zip(src_dx, row["diagnoses"], strict=True):
        assert set(out.keys()) == expected_dx_keys
        src_tx = src.get("treatments") or []
        assert len(out["treatments"]) == len(src_tx)
        for tx in out["treatments"]:
            assert set(tx.keys()) == expected_tx_keys


def test_other_lists_match_source_lengths(one_case: list[dict]) -> None:
    row = clinical.to_patient_rows(one_case)[0]
    src = one_case[0]
    assert len(row["follow_ups"]) == len(src.get("follow_ups") or [])
    assert len(row["exposures"]) == len(src.get("exposures") or [])
    assert len(row["family_histories"]) == len(src.get("family_histories") or [])

    if row["follow_ups"]:
        assert set(row["follow_ups"][0].keys()) == {f.name for f in FOLLOW_UP_FIELDS}
    if row["family_histories"]:
        assert set(row["family_histories"][0].keys()) == {f.name for f in FAMILY_HISTORY_FIELDS}
    # exposures may be empty for CHOL; if present, key set should match
    if row["exposures"]:
        assert set(row["exposures"][0].keys()) == {f.name for f in EXPOSURE_FIELDS}


def test_samples_preserve_full_gdc_hierarchy(one_case: list[dict]) -> None:
    """Each sample carries portions → analytes → aliquots verbatim; counts at
    each level match the source GDC response."""
    row = clinical.to_patient_rows(one_case)[0]
    src_samples = one_case[0].get("samples") or []
    assert len(row["samples"]) == len(src_samples)

    sample_keys = {f.name for f in SAMPLE_FIELDS}
    portion_keys = {f.name for f in PORTION_FIELDS}
    analyte_keys = {f.name for f in ANALYTE_FIELDS}
    aliquot_keys = {f.name for f in ALIQUOT_FIELDS}
    for sample in row["samples"]:
        assert set(sample.keys()) == sample_keys
        for portion in sample["portions"]:
            assert set(portion.keys()) == portion_keys
            for analyte in portion["analytes"]:
                assert set(analyte.keys()) == analyte_keys
                for aliquot in analyte["aliquots"]:
                    assert set(aliquot.keys()) == aliquot_keys

    # Counts at each tree level must match the source.
    src_portion_count = sum(len(s.get("portions") or []) for s in src_samples)
    src_analyte_count = sum(
        len(p.get("analytes") or []) for s in src_samples for p in (s.get("portions") or [])
    )
    src_aliquot_count = sum(
        len(a.get("aliquots") or [])
        for s in src_samples
        for p in (s.get("portions") or [])
        for a in (p.get("analytes") or [])
    )
    out_portion_count = sum(len(s["portions"]) for s in row["samples"])
    out_analyte_count = sum(len(p["analytes"]) for s in row["samples"] for p in s["portions"])
    out_aliquot_count = sum(
        len(a["aliquots"]) for s in row["samples"] for p in s["portions"] for a in p["analytes"]
    )
    assert out_portion_count == src_portion_count
    assert out_analyte_count == src_analyte_count
    assert out_aliquot_count == src_aliquot_count


def test_analyte_type_lives_on_analyte_not_aliquot(one_case: list[dict]) -> None:
    """analyte_type now lives on the analyte (where GDC puts it), not hoisted
    onto each aliquot. This patient has DNA aliquots."""
    row = clinical.to_patient_rows(one_case)[0]
    analyte_types = {
        a.get("analyte_type") for s in row["samples"] for p in s["portions"] for a in p["analytes"]
    }
    analyte_types.discard(None)
    assert analyte_types & {"DNA", "RNA"}, f"expected DNA/RNA in analytes, got {analyte_types}"


def test_samples_sorted_by_days_to_collection_then_id() -> None:
    case = {
        "case_id": "c1",
        "submitter_id": "TCGA-XX-1",
        "project": {"project_id": "TCGA-CHOL"},
        "samples": [
            {"sample_id": "s2", "days_to_collection": 100, "portions": []},
            {"sample_id": "s3", "days_to_collection": None, "portions": []},
            {"sample_id": "s1", "days_to_collection": 50, "portions": []},
        ],
    }
    row = clinical.to_patient_rows([case])[0]
    assert [s["sample_id"] for s in row["samples"]] == ["s1", "s2", "s3"]


def test_write_patients_parquet_round_trip(one_case: list[dict], tmp_path: Path) -> None:
    rows = clinical.to_patient_rows(one_case)
    out = clinical.write_patients(rows, tmp_path, project_id="TCGA-CHOL")
    assert out == tmp_path / "TCGA-CHOL" / "data.parquet"
    assert out.exists()

    table = pq.read_table(out)
    assert table.schema.equals(PATIENTS, check_metadata=False)
    assert table.num_rows == 1

    record = table.to_pylist()[0]
    assert record["case_id"] == one_case[0]["case_id"]
    assert isinstance(record["diagnoses"], list)
    if record["diagnoses"]:
        assert isinstance(record["diagnoses"][0]["treatments"], list)


def test_write_patients_uses_row_groups_and_page_index(
    one_case: list[dict], tmp_path: Path
) -> None:
    """The HF Dataset Viewer needs bounded row groups + a page index to scan
    large parquets without hitting the 300 MB limit. Confirm both are written."""
    # Multiply our one fixture row to span more than one row group.
    n = clinical._ROW_GROUP_SIZE + 5
    rows = clinical.to_patient_rows(one_case * n)
    out = clinical.write_patients(rows, tmp_path, project_id="TCGA-CHOL")

    pf = pq.ParquetFile(out)
    assert pf.metadata.num_rows == n
    # n rows split into row groups of _ROW_GROUP_SIZE → ceil(n / size) groups.
    assert pf.metadata.num_row_groups == 2
    assert pf.metadata.row_group(0).num_rows == clinical._ROW_GROUP_SIZE
    assert pf.metadata.row_group(1).num_rows == 5

    # Page index presence: with `write_page_index=True`, every column should
    # carry both a column index (per-page min/max) and an offset index
    # (per-page byte offsets) the Viewer can use for random access.
    rg0 = pf.metadata.row_group(0)
    assert all(rg0.column(i).has_column_index for i in range(rg0.num_columns))
    assert all(rg0.column(i).has_offset_index for i in range(rg0.num_columns))


def test_diagnoses_sorted_by_days_then_id() -> None:
    case = {
        "case_id": "c1",
        "submitter_id": "TCGA-XX-1",
        "project": {"project_id": "TCGA-CHOL"},
        "diagnoses": [
            {"diagnosis_id": "d2", "days_to_diagnosis": 100.0, "treatments": []},
            {"diagnosis_id": "d3", "days_to_diagnosis": None, "treatments": []},
            {"diagnosis_id": "d1", "days_to_diagnosis": 50.0, "treatments": []},
            {"diagnosis_id": "d0", "days_to_diagnosis": None, "treatments": []},
        ],
    }
    row = clinical.to_patient_rows([case])[0]
    ids = [d["diagnosis_id"] for d in row["diagnoses"]]
    # Temporal first (50, 100), then nulls last sorted by id (d0, d3).
    assert ids == ["d1", "d2", "d0", "d3"]


def test_treatments_sorted_within_each_diagnosis() -> None:
    case = {
        "case_id": "c1",
        "submitter_id": "TCGA-XX-1",
        "project": {"project_id": "TCGA-CHOL"},
        "diagnoses": [
            {
                "diagnosis_id": "d1",
                "days_to_diagnosis": 0.0,
                "treatments": [
                    {"treatment_id": "t-c", "days_to_treatment_start": None},
                    {"treatment_id": "t-a", "days_to_treatment_start": 30.0},
                    {"treatment_id": "t-b", "days_to_treatment_start": 10.0},
                ],
            }
        ],
    }
    row = clinical.to_patient_rows([case])[0]
    tx_ids = [t["treatment_id"] for t in row["diagnoses"][0]["treatments"]]
    assert tx_ids == ["t-b", "t-a", "t-c"]


def test_other_lists_have_deterministic_order() -> None:
    case = {
        "case_id": "c1",
        "submitter_id": "TCGA-XX-1",
        "project": {"project_id": "TCGA-CHOL"},
        "follow_ups": [
            {"follow_up_id": "f3", "days_to_follow_up": 100.0},
            {"follow_up_id": "f1", "days_to_follow_up": 10.0},
            {"follow_up_id": "f2", "days_to_follow_up": None},
        ],
        "exposures": [
            {"exposure_id": "e2"},
            {"exposure_id": "e1"},
        ],
        "family_histories": [
            {"family_history_id": "h2"},
            {"family_history_id": "h1"},
        ],
    }
    row = clinical.to_patient_rows([case])[0]
    assert [f["follow_up_id"] for f in row["follow_ups"]] == ["f1", "f3", "f2"]
    assert [e["exposure_id"] for e in row["exposures"]] == ["e1", "e2"]
    assert [fh["family_history_id"] for fh in row["family_histories"]] == ["h1", "h2"]


def test_to_patient_rows_is_byte_reproducible() -> None:
    """Same source -> same Parquet bytes (modulo order of input cases, which we
    don't sort). With fixed input ordering, repeated builds should be byte-equal."""
    import json

    case = json.loads(FIXTURE.read_text())
    a = clinical.to_patient_rows([case])
    b = clinical.to_patient_rows([case])
    assert a == b


def test_write_card_release_line_single(tmp_path: Path) -> None:
    card = dataset_card.write_card(
        tmp_path,
        projects=["TCGA-CHOL", "TCGA-DLBC"],
        gdc_releases={"TCGA-CHOL": "Data Release 45.0", "TCGA-DLBC": "Data Release 45.0"},
    )
    text = card.read_text()
    assert "**GDC data release:** Data Release 45.0" in text


def test_write_card_release_line_per_project_when_mixed(tmp_path: Path) -> None:
    card = dataset_card.write_card(
        tmp_path,
        projects=["TCGA-CHOL", "TCGA-DLBC"],
        gdc_releases={"TCGA-CHOL": "Data Release 45.0", "TCGA-DLBC": "Data Release 44.0"},
    )
    text = card.read_text()
    assert "**GDC data releases (per project):**" in text
    assert "`TCGA-CHOL`: Data Release 45.0" in text
    assert "`TCGA-DLBC`: Data Release 44.0" in text


def test_write_card_release_line_unknown_when_missing(tmp_path: Path) -> None:
    card = dataset_card.write_card(tmp_path, projects=["TCGA-CHOL"])
    text = card.read_text()
    assert "**GDC data release:** unknown" in text


def _touch_patient_parquets(root: Path, projects: list[str]) -> None:
    """Create the per-project parquet paths `_configs_yaml` checks for.

    Contents don't matter — the card generator only tests existence, so a
    config never points at a file the upload wouldn't carry.
    """
    for project in projects:
        path = root / project / "data.parquet"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()


def test_write_card_emits_one_config_per_project(tmp_path: Path) -> None:
    _touch_patient_parquets(tmp_path, ["TCGA-DLBC", "TCGA-CHOL"])
    card = dataset_card.write_card(tmp_path, projects=["TCGA-DLBC", "TCGA-CHOL"])
    text = card.read_text()
    assert text.startswith("---\n")
    parts = text.split("---\n", 2)
    assert len(parts) == 3
    front = parts[1]

    # One config block per project, sorted alphabetically.
    assert "config_name: TCGA-CHOL" in front
    assert "config_name: TCGA-DLBC" in front
    assert front.index("config_name: TCGA-CHOL") < front.index("config_name: TCGA-DLBC")

    # Per-project data_files path with explicit train split.
    assert "path: TCGA-CHOL/data.parquet" in front
    assert "path: TCGA-DLBC/data.parquet" in front
    assert "split: train" in front


def test_write_card_skips_projects_without_a_parquet(tmp_path: Path) -> None:
    """A project with no parquet on disk gets no config.

    Guards the incremental path: a `--project` build lists every project in
    the tree, and any whose output is absent must not be advertised — HF
    Data Studio surfaces a dangling path as an error rather than an
    omission.
    """
    _touch_patient_parquets(tmp_path, ["TCGA-CHOL"])
    card = dataset_card.write_card(tmp_path, projects=["TCGA-CHOL", "TCGA-DLBC"])
    front = card.read_text().split("---\n", 2)[1]

    assert "config_name: TCGA-CHOL" in front
    assert "TCGA-DLBC" not in front


# ---------------------------------------------------------------------------
# Full `/cases` expansion coverage
#
# GDC exposes 41 expandable groups; we request all of them except the
# `files.*` subtree. That does not fit in one request — the API silently
# truncates a long `expand` list, returning 200 with an empty `hits` — so the
# fetch is two calls merged on `case_id`. These tests pin the merge and, more
# importantly, the guard against that silent failure.
# ---------------------------------------------------------------------------


class _FakeClient:
    """Records each `cases()` call and replays canned halves."""

    def __init__(self, clinical_hits: list[dict], biospecimen_hits: list[dict]) -> None:
        self._clinical = clinical_hits
        self._biospecimen = biospecimen_hits
        self.calls: list[list[str]] = []

    def cases(self, filters, fields, expand, page_size):  # noqa: ANN001
        self.calls.append(expand)
        return self._clinical if "demographic" in expand else self._biospecimen


def test_expansions_cover_every_non_file_group() -> None:
    """The two lists must partition the non-`files` groups, with no overlap."""
    both = clinical.CLINICAL_EXPANSIONS + clinical.BIOSPECIMEN_EXPANSIONS
    assert len(both) == len(set(both)), "an expansion is listed twice"
    assert clinical.EXPANSIONS == both
    # Biospecimen groups are exactly the `samples*` ones.
    assert all(e.startswith("samples") for e in clinical.BIOSPECIMEN_EXPANSIONS)
    assert not any(e.startswith("samples") for e in clinical.CLINICAL_EXPANSIONS)
    # Neither half may exceed the measured API ceiling of 21 groups.
    assert len(clinical.CLINICAL_EXPANSIONS) <= 21
    assert len(clinical.BIOSPECIMEN_EXPANSIONS) <= 21
    # `files.*` stays out: it is ~12x the payload and duplicates `files`.
    assert not any(e.startswith("files") for e in both)


def test_fetch_clinical_merges_the_two_halves() -> None:
    client = _FakeClient(
        clinical_hits=[
            {
                "case_id": "c1",
                "submitter_id": "TCGA-XX-0001",
                "demographic": {"gender": "female"},
            }
        ],
        biospecimen_hits=[{"case_id": "c1", "samples": [{"sample_id": "s1"}]}],
    )
    merged = clinical.fetch_clinical(["TCGA-XYZ"], client)
    assert len(client.calls) == 2
    assert len(merged) == 1
    row = merged[0]
    # Both halves land on the same record.
    assert row["demographic"] == {"gender": "female"}
    assert row["samples"] == [{"sample_id": "s1"}]


def test_fetch_clinical_raises_when_one_half_comes_back_empty() -> None:
    """The signature of GDC silently truncating an over-long `expand`.

    Returning the clinical half alone would ship a dataset whose every case
    has no biospecimen tree, with nothing in the logs to say so.
    """
    client = _FakeClient(
        clinical_hits=[{"case_id": "c1", "demographic": {}}],
        biospecimen_hits=[],
    )
    with pytest.raises(RuntimeError, match="silently truncated"):
        clinical.fetch_clinical(["TCGA-XYZ"], client)


def test_no_cases_at_all_is_not_an_error() -> None:
    """A project with genuinely zero cases returns empty, not a raise."""
    client = _FakeClient(clinical_hits=[], biospecimen_hits=[])
    assert clinical.fetch_clinical(["TCGA-XYZ"], client) == []


def test_newly_expanded_entities_land_on_the_row() -> None:
    """Each new expansion is picked into the row at its own level."""
    case = {
        "case_id": "c1",
        "submitter_id": "TCGA-XX-0001",
        "project": {"project_id": "TCGA-XYZ", "program": {"name": "TCGA", "program_id": "p1"}},
        "annotations": [{"annotation_id": "a1", "category": "Item is noncanonical"}],
        "tissue_source_site": {"tissue_source_site_id": "t1", "code": "4G", "name": "Sapienza"},
        "summary": {
            "file_count": 55,
            "file_size": 123,
            "data_categories": [{"data_category": "Clinical", "file_count": 10}],
            "experimental_strategies": [{"experimental_strategy": "WXS", "file_count": 16}],
        },
        "diagnoses": [
            {
                "diagnosis_id": "d1",
                "annotations": [{"annotation_id": "a2"}],
                "pathology_details": [{"pathology_detail_id": "pd1", "percent_tumor_nuclei": 80.0}],
            }
        ],
        "follow_ups": [
            {
                "follow_up_id": "f1",
                "molecular_tests": [{"molecular_test_id": "m1", "gene_symbol": "IDH1"}],
                "other_clinical_attributes": [
                    {"other_clinical_attribute_id": "o1", "risk_factors": ["Alcohol"]}
                ],
            }
        ],
        "samples": [
            {
                "sample_id": "s1",
                "annotations": [{"annotation_id": "a3"}],
                "portions": [
                    {
                        "portion_id": "p1",
                        "center": {"center_id": "ce1", "name": "BCR"},
                        "slides": [{"slide_id": "sl1", "percent_tumor_cells": 20.0}],
                        "analytes": [
                            {
                                "analyte_id": "an1",
                                "aliquots": [
                                    {
                                        "aliquot_id": "al1",
                                        "center": {"center_id": "ce2", "name": "Broad"},
                                    }
                                ],
                            }
                        ],
                    }
                ],
            }
        ],
    }
    row = clinical.to_patient_rows([case])[0]

    assert row["annotations"][0]["category"] == "Item is noncanonical"
    assert row["tissue_source_site"]["code"] == "4G"
    assert row["program"]["name"] == "TCGA"
    assert row["summary"]["file_count"] == 55
    assert row["summary"]["data_categories"][0]["data_category"] == "Clinical"

    dx = row["diagnoses"][0]
    assert dx["annotations"][0]["annotation_id"] == "a2"
    assert dx["pathology_details"][0]["percent_tumor_nuclei"] == 80.0

    fu = row["follow_ups"][0]
    assert fu["molecular_tests"][0]["gene_symbol"] == "IDH1"
    # GDC types this `keyword` but returns an array; the schema must match.
    assert fu["other_clinical_attributes"][0]["risk_factors"] == ["Alcohol"]

    sample = row["samples"][0]
    assert sample["annotations"][0]["annotation_id"] == "a3"
    portion = sample["portions"][0]
    assert portion["center"]["name"] == "BCR"
    assert portion["slides"][0]["percent_tumor_cells"] == 20.0
    assert portion["analytes"][0]["aliquots"][0]["center"]["name"] == "Broad"


def test_expanded_row_conforms_to_the_published_cases_schema() -> None:
    """A row with every new entity populated must fit TABULAR_TABLES['cases']."""
    import pyarrow as pa
    from tcga2hf.schema import TABULAR_TABLES
    from tcga2hf_pipeline import tabular

    case = {
        "case_id": "c1",
        "submitter_id": "TCGA-XX-0001",
        "project": {"project_id": "TCGA-XYZ", "program": {"name": "TCGA"}},
        "annotations": [{"annotation_id": "a1"}],
        "tissue_source_site": {"tissue_source_site_id": "t1", "code": "4G"},
        "summary": {
            "file_count": 1,
            "file_size": 2,
            "data_categories": [],
            "experimental_strategies": [],
        },
        "diagnoses": [
            {
                "diagnosis_id": "d1",
                "pathology_details": [{"pathology_detail_id": "pd1"}],
            }
        ],
        "follow_ups": [
            {
                "follow_up_id": "f1",
                "molecular_tests": [{"molecular_test_id": "m1"}],
                "other_clinical_attributes": [
                    {
                        "other_clinical_attribute_id": "o1",
                        "risk_factors": ["Alcohol"],
                        "comorbidities": ["Diabetes"],
                        "viral_hepatitis_serology_tests": ["HBV"],
                    }
                ],
            }
        ],
        "samples": [],
    }
    rows = tabular._cases_rows([case])
    table = pa.Table.from_pylist(rows, schema=TABULAR_TABLES["cases"])
    assert table.num_rows == 1
