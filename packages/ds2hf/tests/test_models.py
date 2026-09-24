"""Tests for the TcgaHfPatient pydantic reference implementation.

Two flavors:
  - synthetic: tiny in-memory rows that exercise specific transformations.
  - live: round-trip the actual rebuilt CHOL parquet through TcgaHfPatient and
    confirm the rich convenience methods work end-to-end against real data.
"""

from __future__ import annotations

from pathlib import Path

import pyarrow.parquet as pq
import pytest
from ds2hf.tcga import schema
from ds2hf.tcga.models import (
    Aliquot,
    Analyte,
    Demographic,
    Diagnosis,
    Exposure,
    FamilyHistory,
    FollowUp,
    GeneExpression,
    Mutation,
    PathologyReport,
    Portion,
    Sample,
    TcgaHfPatient,
    Treatment,
)
from pydantic import ValidationError

PROCESSED = Path.home() / "data/tcga2hf/processed"
LIVE_REQUIRED = pytest.mark.skipif(
    not (PROCESSED / "TCGA-CHOL/data.parquet").exists(),
    reason="run `ds2hf-pipeline tcga build` first to populate $HOME/data/tcga2hf/processed",
)


# ---------------------------------------------------------------------------
# Synthetic-data unit tests
# ---------------------------------------------------------------------------


def _mk_patient_with_pair() -> TcgaHfPatient:
    """A minimal patient with one tumor sample, one normal sample, and one
    mutation linking them — enough to exercise tumor_normal_pairs +
    aliquot_to_sample."""
    return TcgaHfPatient(
        case_id="case-1",
        case_submitter_id="TCGA-XX-1",
        project_id="TCGA-CHOL",
        samples=[
            Sample(
                sample_id="tumor-s",
                tissue_type="Tumor",
                portions=[
                    Portion(
                        portion_id="tumor-p",
                        is_ffpe=False,
                        analytes=[
                            Analyte(
                                analyte_id="tumor-a",
                                analyte_type="DNA",
                                aliquots=[Aliquot(aliquot_id="tumor-aq")],
                            )
                        ],
                    )
                ],
            ),
            Sample(
                sample_id="normal-s",
                tissue_type="Normal",
                portions=[
                    Portion(
                        portion_id="normal-p",
                        analytes=[
                            Analyte(
                                analyte_id="normal-a",
                                analyte_type="DNA",
                                aliquots=[Aliquot(aliquot_id="normal-aq")],
                            )
                        ],
                    )
                ],
            ),
        ],
        samples_masked_somatic_mutation=[
            Mutation.model_validate(
                {
                    "tumor_sample_id": "tumor-s",
                    "matched_normal_sample_id": "normal-s",
                    "Hugo_Symbol": "TP53",
                    "Variant_Classification": "Missense_Mutation",
                    "Tumor_Sample_UUID": "tumor-aq",
                    "Matched_Norm_Sample_UUID": "normal-aq",
                }
            )
        ],
    )


def test_strict_mode_rejects_unknown_fields() -> None:
    """`extra="forbid"` is the contract: a row with a key not in the GDC
    dictionary must fail to validate, not be silently dropped."""
    with pytest.raises(ValidationError):
        TcgaHfPatient(
            case_id="case-1",
            case_submitter_id="TCGA-XX-1",
            project_id="TCGA-CHOL",
            not_a_real_field="oops",  # type: ignore[call-arg]
        )
    with pytest.raises(ValidationError):
        Sample.model_validate(
            {"sample_id": "s-1", "made_up_attribute": "oops"},
        )


def test_pydantic_fields_match_pa_fields_exactly() -> None:
    """The pydantic field sets must equal the GDC dictionary (`*_FIELDS`)
    field sets. If gdcdictionary adds/removes a field upstream, the pyarrow
    schema is regenerated; this test fails until pydantic catches up.

    Designed-in exception: TcgaHfPatient declares `clinical_supplement` as
    a flex `dict` field outside PATIENT_FIELDS because BCR biotab columns
    vary by cancer type — they're inferred per project at parquet write
    time rather than enumerated in the global schema.
    """
    # (cls, fields, extra_pyd_fields_intentionally_outside_schema)
    pairs: list[tuple[type, list, set[str]]] = [
        (Aliquot, schema.ALIQUOT_FIELDS, set()),
        (Analyte, schema.ANALYTE_FIELDS, set()),
        (Portion, schema.PORTION_FIELDS, set()),
        (Sample, schema.SAMPLE_FIELDS, set()),
        (Demographic, schema.DEMOGRAPHIC_FIELDS, set()),
        (Treatment, schema.TREATMENT_FIELDS, set()),
        (Diagnosis, schema.DIAGNOSIS_FIELDS, set()),
        (FollowUp, schema.FOLLOW_UP_FIELDS, set()),
        (Exposure, schema.EXPOSURE_FIELDS, set()),
        (FamilyHistory, schema.FAMILY_HISTORY_FIELDS, set()),
        (Mutation, schema.MUTATION_FIELDS, set()),
        (GeneExpression, schema.EXPRESSION_FIELDS, set()),
        (PathologyReport, schema.PATHOLOGY_REPORT_FIELDS, set()),
        (TcgaHfPatient, schema.PATIENT_FIELDS, {"clinical_supplement"}),
    ]
    for cls, fields, extras in pairs:
        pyd_names = set(cls.model_fields) - extras
        pa_names = {f.name for f in fields}
        assert pyd_names == pa_names, (
            f"{cls.__name__}: pydantic - pa = {pyd_names - pa_names}; "
            f"pa - pydantic = {pa_names - pyd_names}"
        )


def _mk_patient_with_report(*, pdf: bytes = b"%PDF-1.4 fake") -> TcgaHfPatient:
    """A patient with one tumor sample and its pathology report attached."""
    return TcgaHfPatient(
        case_id="case-1",
        case_submitter_id="TCGA-XX-1",
        project_id="TCGA-CHOL",
        samples=[
            Sample(
                sample_id="tumor-s",
                submitter_id="TCGA-XX-1-01A",
                tissue_type="Tumor",
                pathology_report_uuid="REPORT-UUID",
            )
        ],
        samples_pathology_report=[
            PathologyReport(
                sample_id="tumor-s",
                sample_submitter_id="TCGA-XX-1-01A",
                pathology_report_uuid="REPORT-UUID",
                source_file_id="file-1",
                file_name="TCGA-XX-1.REPORT-UUID.PDF",
                file_size=len(pdf),
                pdf_bytes=pdf,
            )
        ],
    )


def test_pathology_report_bytes_survive_the_model() -> None:
    """pa.binary() -> pydantic `bytes` with no re-encoding."""
    p = _mk_patient_with_report()
    (report,) = p.samples_pathology_report
    assert report.pdf_bytes == b"%PDF-1.4 fake"
    assert p.pathology_reports_by_sample() == {"tumor-s": [report]}


def test_write_pdf_round_trips_to_disk(tmp_path: Path) -> None:
    p = _mk_patient_with_report()
    out = p.samples_pathology_report[0].write_pdf(tmp_path / "nested" / "report.pdf")
    assert out.read_bytes() == b"%PDF-1.4 fake"


def test_write_pdf_raises_when_bytes_were_projected_away() -> None:
    """Column projection is the expected way to read these rows cheaply, so
    the failure has to name that cause rather than raise AttributeError."""
    report = PathologyReport(source_file_id="file-1")
    with pytest.raises(ValueError, match="pdf_bytes"):
        report.write_pdf("unused.pdf")


def test_samples_missing_pathology_report_flags_unfetched_reports() -> None:
    """GDC populates sample.pathology_report_uuid whether or not we fetched
    the PDF, so the gap has to be visible to a consumer."""
    p = _mk_patient_with_report()
    assert p.samples_missing_pathology_report() == []

    p.samples_pathology_report = []
    assert [s.sample_id for s in p.samples_missing_pathology_report()] == ["tumor-s"]


def test_aliquot_to_sample_walks_full_tree() -> None:
    p = _mk_patient_with_pair()
    assert p.aliquot_to_sample() == {"tumor-aq": "tumor-s", "normal-aq": "normal-s"}


def test_aliquot_lookup_returns_full_provenance_chain() -> None:
    p = _mk_patient_with_pair()
    chain = p.aliquot_lookup()["tumor-aq"]
    sample, portion, analyte, aliquot = chain
    assert sample.sample_id == "tumor-s"
    assert portion.portion_id == "tumor-p"
    assert analyte.analyte_id == "tumor-a"
    assert analyte.analyte_type == "DNA"
    assert aliquot.aliquot_id == "tumor-aq"


def test_tumor_and_normal_sample_filters() -> None:
    p = _mk_patient_with_pair()
    assert [s.sample_id for s in p.tumor_samples()] == ["tumor-s"]
    assert [s.sample_id for s in p.normal_samples()] == ["normal-s"]


def test_tumor_normal_pairs_dedupes() -> None:
    p = _mk_patient_with_pair()
    # Add a second mutation pointing at the same pair
    p.samples_masked_somatic_mutation.append(
        Mutation.model_validate(
            {
                "tumor_sample_id": "tumor-s",
                "matched_normal_sample_id": "normal-s",
                "Hugo_Symbol": "KRAS",
                "Tumor_Sample_UUID": "tumor-aq",
                "Matched_Norm_Sample_UUID": "normal-aq",
            }
        )
    )
    pairs = p.tumor_normal_pairs()
    assert len(pairs) == 1
    assert pairs[0][0].sample_id == "tumor-s"
    assert pairs[0][1].sample_id == "normal-s"


def test_mutations_grouped_by_gene_and_consequence() -> None:
    p = _mk_patient_with_pair()
    p.samples_masked_somatic_mutation.append(
        Mutation.model_validate(
            {
                "tumor_sample_id": "tumor-s",
                "matched_normal_sample_id": "normal-s",
                "Hugo_Symbol": "TP53",
                "Variant_Classification": "Silent",
                "Tumor_Sample_UUID": "tumor-aq",
                "Matched_Norm_Sample_UUID": "normal-aq",
            }
        )
    )
    by_gene = p.mutations_by_gene()
    assert set(by_gene) == {"TP53"}
    assert len(by_gene["TP53"]) == 2
    assert p.mutations_by_consequence() == {"Missense_Mutation": 1, "Silent": 1}


def test_timeline_includes_clinical_and_biospecimen_on_same_anchor() -> None:
    p = _mk_patient_with_pair()
    p.days_to_consent = -10  # consent 10 days before diagnosis
    p.diagnoses = [
        Diagnosis(
            diagnosis_id="d1",
            days_to_diagnosis=0,
            primary_diagnosis="Cholangiocarcinoma",
            treatments=[
                Treatment(treatment_id="t1", days_to_treatment_start=30, treatment_type="Surgery"),
                Treatment(
                    treatment_id="t2", days_to_treatment_start=10, treatment_type="Pharmaceutical"
                ),
            ],
        )
    ]
    p.samples[0].days_to_sample_procurement = 5  # tumor surgery day 5
    p.samples[0].days_to_collection = 1500  # BCR received the case package on day 1500
    p.samples[1].days_to_collection = 1500  # same per-case batch event

    events = p.timeline()
    days = [e.day for e in events]
    assert days == sorted(days)

    # All categories are present
    cats = [e.category for e in events]
    assert "consent" in cats
    assert "diagnosis" in cats
    assert "treatment_start" in cats
    assert "sample_procurement" in cats
    assert "bcr_receipt" in cats

    # Two BCR receipt events at the same day (one per sample)
    assert cats.count("bcr_receipt") == 2

    # Consent (day -10) before diagnosis (day 0)
    assert events[0].category == "consent"
    assert events[1].category == "diagnosis"


# ---------------------------------------------------------------------------
# Live round-trip against the real CHOL parquet
# ---------------------------------------------------------------------------


@LIVE_REQUIRED
def test_round_trip_real_chol_patient() -> None:
    rows = pq.read_table(PROCESSED / "TCGA-CHOL/data.parquet").to_pylist()
    target = next(r for r in rows if r["case_submitter_id"] == "TCGA-W5-AA39")
    p = TcgaHfPatient.model_validate(target)

    # Basic identity
    assert p.case_submitter_id == "TCGA-W5-AA39"
    assert p.project_id == "TCGA-CHOL"
    assert p.demographic is not None
    assert p.demographic.sex_at_birth == "male"

    # Biospecimen tree
    assert len(p.samples) > 0
    assert len(p.all_samples_by_id()) == len(p.samples)
    aliquot_map = p.aliquot_to_sample()
    assert all(sid in p.all_samples_by_id() for sid in aliquot_map.values())

    # Mutations
    by_gene = p.mutations_by_gene()
    assert len(by_gene) > 0  # CHOL patients all have mutations

    # Tumor/normal pairs from real mutations
    pairs = p.tumor_normal_pairs()
    assert len(pairs) >= 1
    for tumor, normal in pairs:
        assert tumor.is_tumor
        assert normal.is_normal

    # Expression — TCGA-W5-AA39 has 1 RNA-Seq aliquot
    assert len(p.samples_gene_expression_quantification) == 1
    er = p.samples_gene_expression_quantification[0]
    assert len(er.gene_id) == 60660
    # ALB is a known top-expressed gene in liver/biliary tissue
    alb = p.expression_for_gene("ALB")
    assert er.aliquot_id in alb
    assert alb[er.aliquot_id]["tpm_unstranded"] is not None

    # Timeline reconstructs without crash, sorted ascending, and includes
    # both clinical and biospecimen events on the same anchor.
    events = p.timeline()
    assert events == sorted(events, key=lambda e: e.day)
    categories = {e.category for e in events}
    assert "diagnosis" in categories or "follow_up" in categories
    # bcr_receipt should appear since CHOL data has days_to_collection
    assert "bcr_receipt" in categories
