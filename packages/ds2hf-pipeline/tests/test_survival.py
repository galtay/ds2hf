from __future__ import annotations

import json
from pathlib import Path

import pytest
from ds2hf_pipeline.tcga import clinical, survival

FIXTURE = Path(__file__).parent / "fixtures" / "case_chol_one.json"


@pytest.fixture
def chol_row() -> dict:
    """One CHOL patient row, the same fixture used by `test_clinical`.

    `TCGA-ZH-A8Y3` — Dead at day 561, with a recorded recurrence and
    margin-status follow-ups. Liu's CDR records:
        OS=(1, 561), DSS=(1, 561), PFI=(1, 561), DFI=NA  (CHOL with no R0)
    """
    case = json.loads(FIXTURE.read_text())
    return clinical.to_patient_rows([case])[0]


def test_attach_survival_populates_struct(chol_row: dict) -> None:
    """Sanity: attach_survival writes a `survival_derived` struct with the 8 fields."""
    rows = [chol_row]
    survival.attach_survival(rows)
    assert "survival_derived" in rows[0]
    sd = rows[0]["survival_derived"]
    expected = {
        "os_event",
        "os_time",
        "dss_event",
        "dss_time",
        "pfi_event",
        "pfi_time",
        "dfi_event",
        "dfi_time",
    }
    assert expected <= sd.keys()


def test_canary_case_matches_cdr(chol_row: dict) -> None:
    """TCGA-ZH-A8Y3 should produce the same OS/DSS/PFI/DFI values as Liu's CDR.

    Locks in the algorithm against the canary patient: any future change
    to `survival.py` that breaks this is a regression. CHOL falls into
    the residual_disease fallback chain; this patient lacks R0 so DFI
    is not derivable (matches CDR's NA).
    """
    rows = [chol_row]
    survival.attach_survival(rows)
    sd = rows[0]["survival_derived"]
    assert (sd["os_event"], sd["os_time"]) == (1, 561.0)
    assert (sd["dss_event"], sd["dss_time"]) == (1, 561.0)
    assert (sd["pfi_event"], sd["pfi_time"]) == (1, 561.0)
    assert (sd["dfi_event"], sd["dfi_time"]) == (None, None)


def test_short_project_strips_tcga_prefix() -> None:
    assert survival._short_project("TCGA-CHOL") == "CHOL"
    assert survival._short_project(None) is None
    assert survival._short_project("CHOL") == "CHOL"  # already short


def test_is_stage_iv_matches_only_iv_stages() -> None:
    """Stage IV variants match; III/IIIA do not (no Roman-numeral collision)."""

    def case_with_stage(stage: str | None) -> dict:
        return {
            "diagnoses": [{"diagnosis_is_primary_disease": True, "ajcc_pathologic_stage": stage}]
        }

    assert survival._is_stage_iv(case_with_stage("Stage IV"))
    assert survival._is_stage_iv(case_with_stage("Stage IVA"))
    assert survival._is_stage_iv(case_with_stage("Stage IVB"))
    assert not survival._is_stage_iv(case_with_stage("Stage III"))
    assert not survival._is_stage_iv(case_with_stage("Stage IIIA"))
    assert not survival._is_stage_iv(case_with_stage("Stage I"))
    assert not survival._is_stage_iv(case_with_stage(None))
    assert not survival._is_stage_iv({"diagnoses": []})


def test_dfi_excluded_for_no_field_tumor_types() -> None:
    """SKCM/THYM/UVM/LAML always return DFI=NA regardless of patient state."""
    case = {
        "demographic": {"vital_status": "Alive", "days_to_death": None},
        "diagnoses": [{"diagnosis_is_primary_disease": True, "days_to_last_follow_up": 365}],
        "follow_ups": [],
    }
    for proj in ("TCGA-SKCM", "TCGA-THYM", "TCGA-UVM", "TCGA-LAML"):
        assert survival.derive_dfi(case, proj) == (None, None), proj


def test_dfi_excluded_for_stage_iv() -> None:
    case = {
        "demographic": {"vital_status": "Alive", "days_to_death": None},
        "diagnoses": [
            {
                "diagnosis_is_primary_disease": True,
                "ajcc_pathologic_stage": "Stage IV",
                "days_to_last_follow_up": 365,
                "residual_disease": "R0",
            }
        ],
        "follow_ups": [],
    }
    assert survival.derive_dfi(case, "TCGA-CHOL") == (None, None)


def test_pfi_uses_progression_or_recurrence_yes_at_follow_up_day() -> None:
    """Cases where the GDC records progression_or_recurrence=Yes without
    a dedicated days_to_recurrence/progression should still produce a PFI
    event at the follow-up's days_to_follow_up. Mirrors TCGA-UY-A78L."""
    case = {
        "demographic": {"vital_status": "Alive", "days_to_death": None},
        "diagnoses": [{"diagnosis_is_primary_disease": True, "days_to_last_follow_up": 1127}],
        "follow_ups": [
            {
                "days_to_follow_up": 355,
                "progression_or_recurrence": "Yes",
                "progression_or_recurrence_type": "Unknown",
            },
            {"days_to_follow_up": 1127, "disease_response": "WT-With Tumor"},
        ],
    }
    assert survival.derive_pfi(case) == (1, 355.0)


def test_dfi_subsequent_primary_same_organ_is_event() -> None:
    """Subsequent primary in the same organ as the original tumor is treated
    as a same-disease event by Liu (DFI=1), not censored."""
    case = {
        "demographic": {"vital_status": "Alive", "days_to_death": None},
        "diagnoses": [
            {
                "diagnosis_is_primary_disease": True,
                "tissue_or_organ_of_origin": "Breast, NOS",
                "residual_disease": "R0",
                "days_to_last_follow_up": 5000,
            },
            {
                "diagnosis_is_primary_disease": False,
                "classification_of_tumor": "Subsequent Primary",
                "tissue_or_organ_of_origin": "Breast, NOS",
                "days_to_diagnosis": 3076,
            },
        ],
        "follow_ups": [],
    }
    assert survival.derive_dfi(case, "TCGA-LIHC") == (1, 3076.0)


def test_dfi_subsequent_primary_different_organ_is_censored() -> None:
    """Subsequent primary in a different organ is Liu's 'new primary in other
    organ' censoring case — DFI stays at the censoring time, not the event."""
    case = {
        "demographic": {"vital_status": "Alive", "days_to_death": None},
        "diagnoses": [
            {
                "diagnosis_is_primary_disease": True,
                "tissue_or_organ_of_origin": "Liver",
                "residual_disease": "R0",
                "days_to_last_follow_up": 854,
            },
            {
                "diagnosis_is_primary_disease": False,
                "classification_of_tumor": "Subsequent Primary",
                "tissue_or_organ_of_origin": "Kidney, NOS",
                "days_to_diagnosis": 110,
            },
        ],
        "follow_ups": [],
    }
    # DFI=0, time=854 (last follow-up); the day-110 subsequent primary was censored.
    assert survival.derive_dfi(case, "TCGA-LIHC") == (0, 854.0)


def test_dss_dead_with_tumor_free_at_death_returns_zero() -> None:
    """Liu's spec: Dead AND tumor_status=TUMOR FREE -> DSS=0, regardless of cause."""
    case = {
        "demographic": {
            "vital_status": "Dead",
            "days_to_death": 500,
            "cause_of_death": "Cancer Related",
        },
        "diagnoses": [{"diagnosis_is_primary_disease": True, "days_to_last_follow_up": 500}],
        "follow_ups": [
            {"days_to_follow_up": 500, "disease_response": "TF-Tumor Free"},
        ],
    }
    assert survival.derive_dss(case) == (0, 500.0)


def test_dss_dead_with_unknown_status_defaults_to_event() -> None:
    """Dead patients with no TF signal default to DSS=1 (the conservative
    survival-analysis default that matches Liu's CDR pattern for Dead+unknown)."""
    case = {
        "demographic": {"vital_status": "Dead", "days_to_death": 365, "cause_of_death": None},
        "diagnoses": [{"diagnosis_is_primary_disease": True, "days_to_last_follow_up": 365}],
        "follow_ups": [],
    }
    assert survival.derive_dss(case) == (1, 365.0)
