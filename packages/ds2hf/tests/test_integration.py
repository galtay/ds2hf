"""Integration tests against locally built parquets (skip if absent).

Stress-tests `TcgaHfPatient.model_validate` over every patient row in every
project parquet, then exercises the timeline + consistency_check methods
across the cohort. Several of these tests assert the *presence* of a known
GDC quirk (sample collections post-dating death), which means if GDC ever
fixes the underlying anchor inconsistency, those tests will fail and we'll
know to update our docs.

Run only these:                  uv run pytest -m integration
Skip these in regular runs:      uv run pytest -m "not integration"
"""

from __future__ import annotations

import time
from collections import Counter
from pathlib import Path

import pyarrow.parquet as pq
import pytest
from ds2hf.tcga.models import TcgaHfPatient

PROCESSED = Path.home() / "data/tcga2hf/processed"

pytestmark = [
    pytest.mark.integration,
    pytest.mark.skipif(
        not (PROCESSED / "TCGA-CHOL/data.parquet").exists(),
        reason="run `ds2hf-pipeline tcga build` first to populate $HOME/data/tcga2hf/processed",
    ),
]


@pytest.fixture(scope="module")
def all_patients() -> list[TcgaHfPatient]:
    """Validate every row from every project parquet through TcgaHfPatient."""
    patients: list[TcgaHfPatient] = []
    for project_dir in sorted(PROCESSED.glob("TCGA-*")):
        rows = pq.read_table(project_dir / "data.parquet").to_pylist()
        for row in rows:
            patients.append(TcgaHfPatient.model_validate(row))
    return patients


def test_validate_every_row(all_patients: list[TcgaHfPatient]) -> None:
    assert len(all_patients) > 0
    # case_id is unique per patient
    case_ids = {p.case_id for p in all_patients}
    assert len(case_ids) == len(all_patients)
    # case_submitter_id follows TCGA barcode prefix
    assert all(p.case_submitter_id.startswith("TCGA-") for p in all_patients)


def test_validation_throughput_is_reasonable(all_patients: list[TcgaHfPatient]) -> None:
    """Re-validate a slice and confirm pydantic isn't slowing to a crawl with
    nested 60k-element gene vectors. Loose ceiling — we just want to catch
    20x regressions, not benchmark."""
    sample = all_patients[: min(10, len(all_patients))]
    payloads = [p.model_dump() for p in sample]
    start = time.time()
    for payload in payloads:
        TcgaHfPatient.model_validate(payload)
    avg_ms = (time.time() - start) / len(payloads) * 1000
    assert avg_ms < 500, f"validation slowed to {avg_ms:.0f} ms/patient"


# ---------------------------------------------------------------------------
# Anchor invariants we expect to hold for every patient
# ---------------------------------------------------------------------------


def test_days_to_birth_is_always_negative(all_patients: list[TcgaHfPatient]) -> None:
    """days_to_birth must be negative when present — birth precedes index."""
    for p in all_patients:
        if p.demographic and p.demographic.days_to_birth is not None:
            assert p.demographic.days_to_birth < 0, p.case_submitter_id


def test_primary_diagnosis_days_to_diagnosis_is_mostly_zero(
    all_patients: list[TcgaHfPatient],
) -> None:
    """The primary diagnosis IS the index event per the GDC dictionary, so
    `days_to_diagnosis == 0` should hold for almost every primary
    diagnosis. At full TCGA scale a small fraction of cases violate this
    (e.g. patients with synchronous malignancies whose second 'primary
    disease' diagnosis is dated at non-zero days). Assert the rule holds
    on the vast majority of the cohort rather than 100%."""
    total_primary = 0
    nonzero = 0
    for p in all_patients:
        for dx in p.diagnoses:
            if dx.diagnosis_is_primary_disease and dx.days_to_diagnosis is not None:
                total_primary += 1
                if dx.days_to_diagnosis != 0:
                    nonzero += 1
    if total_primary == 0:
        pytest.skip("no primary diagnoses with days_to_diagnosis in cohort")
    rate = nonzero / total_primary
    assert rate < 0.05, (
        f"{nonzero}/{total_primary} primary diagnoses ({rate:.1%}) have non-zero "
        "days_to_diagnosis — exceeds the 5% tolerance for this GDC quirk"
    )


def test_age_at_diagnosis_uses_diagnosis_anchor_not_index(
    all_patients: list[TcgaHfPatient],
) -> None:
    """For non-primary diagnoses, age_at_diagnosis should equal
    -days_to_birth + days_to_diagnosis. Confirms age_at_diagnosis is anchored to
    the diagnosis date, not the index date."""
    for p in all_patients:
        if not (p.demographic and p.demographic.days_to_birth is not None):
            continue
        age_at_index = -p.demographic.days_to_birth
        for dx in p.diagnoses:
            if dx.age_at_diagnosis is None or dx.days_to_diagnosis is None:
                continue
            expected = age_at_index + int(dx.days_to_diagnosis)
            assert dx.age_at_diagnosis == expected, (
                f"{p.case_submitter_id} age_at_diagnosis={dx.age_at_diagnosis} "
                f"!= -days_to_birth ({age_at_index}) + days_to_diagnosis "
                f"({dx.days_to_diagnosis}) = {expected}"
            )


# ---------------------------------------------------------------------------
# Known GDC quirks — assert their *presence* so we notice if upstream changes
# ---------------------------------------------------------------------------


def test_sample_collection_uses_different_anchor_than_other_events(
    all_patients: list[TcgaHfPatient],
) -> None:
    """Empirical: many TCGA samples have days_to_collection > days_to_death,
    proving days_to_collection is not anchored to the same index date as the
    other days_to_* fields. If this stops being true (GDC fixes the anchor),
    revisit our `TimelineEvent.anchor` documentation."""
    post_death = []
    for p in all_patients:
        dod = p.demographic.days_to_death if p.demographic else None
        if dod is None:
            continue
        for s in p.samples:
            dtc = s.days_to_collection
            if dtc is not None and dtc > dod:
                post_death.append((p.case_submitter_id, dtc - dod))
    assert len(post_death) > 0, (
        "expected GDC's days_to_collection anchor quirk to surface in our cohort"
    )


def test_consistency_check_method_works_per_patient(
    all_patients: list[TcgaHfPatient],
) -> None:
    """consistency_check returns the expected keys and known cohort signals."""
    expected_keys = {
        "bcr_receipts_after_death",
        "samples_with_no_temporal_data",
        "pre_index_treatments",
        "non_diagnosis_index",
    }
    cohort_total: Counter[str] = Counter()
    for p in all_patients:
        report = p.consistency_check()
        assert set(report) == expected_keys
        cohort_total.update(report)
    # BCR receipts post-dating death are common in TCGA — preserved tissue
    # shipped years after death, not an anchor inconsistency.
    assert cohort_total["bcr_receipts_after_death"] > 0
    # Sample-only patient registrations exist.
    assert cohort_total["samples_with_no_temporal_data"] > 0


def test_index_date_is_diagnosis_for_tcga(all_patients: list[TcgaHfPatient]) -> None:
    """Empirical: every TCGA case in our cohort uses 'Diagnosis' as the index_date
    (or has it null). Confirms the GDC-documented TCGA convention; if this ever
    changes, downstream timeline math may need re-anchoring per case."""
    for p in all_patients:
        if p.index_date is not None:
            assert p.index_date == "Diagnosis", (p.case_submitter_id, p.index_date)


def test_timeline_includes_every_dated_category(
    all_patients: list[TcgaHfPatient],
) -> None:
    """timeline() now returns clinical AND biospecimen events on the unified
    index anchor. Confirm the cohort actually surfaces every category we
    document, and that all events are sorted ascending."""
    seen: set[str] = set()
    for p in all_patients:
        events = p.timeline()
        days = [e.day for e in events]
        assert days == sorted(days), p.case_submitter_id
        for ev in events:
            seen.add(ev.category)

    # Categories that should appear in CHOL+DLBC
    for required in {"diagnosis", "treatment_start", "follow_up", "bcr_receipt", "death"}:
        assert required in seen, f"missing category {required}; saw {seen}"


# ---------------------------------------------------------------------------
# Cohort-level timeline shape
# ---------------------------------------------------------------------------


def test_timelines_sorted_ascending_when_present(all_patients: list[TcgaHfPatient]) -> None:
    """`timeline()` always returns events sorted ascending by day. Some
    patients have empty timelines (see test below) — that's a separate
    concern."""
    for p in all_patients:
        days = [e.day for e in p.timeline()]
        assert days == sorted(days), p.case_submitter_id


def test_some_patients_have_no_clinical_metadata(
    all_patients: list[TcgaHfPatient],
) -> None:
    """Real GDC data quirk: some TCGA cases are registered with biospecimens
    but no clinical metadata — empty `diagnoses`, empty `follow_ups`,
    `demographic is None`, and `days_to_collection` is null on every
    sample (so they produce empty timelines). Surface this so dataset
    consumers don't assume every row has clinical context.

    The `no_timeline` and `no_clinical` sets *mostly* coincide (an empty
    timeline almost always means no clinical metadata), but at full TCGA
    scale a handful of patients diverge — e.g. samples with
    `days_to_collection` populated despite missing clinical metadata
    (counts toward a timeline event), or clinical metadata where every
    `days_to_*` happens to be null (no timeline). Both directions are
    rare; we assert the overlap is large rather than exact.
    """
    no_timeline = {p.case_submitter_id for p in all_patients if not p.timeline()}
    no_clinical = {
        p.case_submitter_id
        for p in all_patients
        if not p.diagnoses and not p.follow_ups and p.demographic is None
    }
    sym_diff = no_timeline.symmetric_difference(no_clinical)
    union = no_timeline | no_clinical
    if union:
        # ≥ 90% of the patients in either set are in both sets.
        assert len(sym_diff) <= 0.1 * len(union), (
            f"no_timeline ⊕ no_clinical = {len(sym_diff)} / {len(union)} — "
            "divergence between empty-timeline and no-clinical-metadata sets is "
            "larger than expected"
        )
    # At least one such patient exists, and they're a small minority of the cohort.
    assert 0 < len(no_clinical) < 0.25 * len(all_patients), (
        f"{len(no_clinical)}/{len(all_patients)} no-clinical patients — unexpected proportion"
    )
