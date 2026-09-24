"""Re-derive Liu et al. 2018 CDR survival endpoints from the live GDC data.

Companion to `cdr.py`: where `cdr.attach_cdr` lifts curated values from
Liu's frozen 2018 workbook, this module re-computes the same four
endpoints (OS, DSS, DFI, PFI) from the patient row's clinical structure
using Liu's documented algorithm. The two streams together let us:

  1. Validate the pipeline by comparing re-derived values to CDR for the
     ~11,160 patients Liu covers, and
  2. Extend coverage to the ~268 post-2018-freeze patients CDR can't reach.

## Liu's algorithm (verbatim, from the TCGA-CDR_Notes sheet)

  - **OS**: 1 for death from any cause, 0 for alive.
    OS.time: last_contact_days_to or death_days_to, whichever is larger.

  - **DSS**: 1 for `vital_status == Dead AND tumor_status == WITH TUMOR`,
    or `cause_of_death` indicates the cancer; 0 for Alive or
    `Dead AND tumor_status == TUMOR FREE`. Else NA.
    DSS.time: same as OS.time.

  - **DFI**: 1 for new tumor event (recurrence/metastasis/new primary,
    including type N/A). Disease-free at end of first course is required:
    `treatment_outcome_first_course == "Complete Remission/Response"`,
    falling back to `residual_tumor == "R0"` for tumor types without
    that field, falling back to `margin_status == "negative"` for those
    without either. Tumor types with none of the three: DFI is NA.
    "New primary in other organ" is censored, not an event. Stage IV
    excluded; dead-with-tumor-no-event excluded. 0 for censored otherwise.
    DFI.time: new_tumor_event_dx_days_to for events, else
    last_contact_days_to or death_days_to (whichever applies).

  - **PFI**: 1 for new tumor event (any kind, including type N/A) OR
    Dead-with-cancer-without-new-tumor-event. 0 for censored otherwise.
    PFI.time: new_tumor_event_dx_days_to or death_days_to for events;
    last_contact_days_to or death_days_to for censored.

## Field mapping (Liu's old TCGA name -> modern GDC path in our schema)

  - vital_status                       demographic.vital_status
  - death_days_to                      demographic.days_to_death
  - cause_of_death                     demographic.cause_of_death
  - last_contact_days_to               max(diagnoses[primary].days_to_last_follow_up,
                                            max(follow_ups[].days_to_follow_up))
  - tumor_status                       latest non-null follow_ups[].disease_response,
                                            mapped 'WT-With Tumor' -> WITH TUMOR,
                                            'TF-Tumor Free' -> TUMOR FREE
  - new_tumor_event_dx_days_to         min(populated follow_ups[].days_to_recurrence
                                            or .days_to_progression)
  - new_tumor_event_type               follow_ups[].progression_or_recurrence_type
                                            at the earliest event
  - treatment_outcome_first_course     diagnoses[primary].treatments[0].treatment_outcome
                                            (treatments pre-sorted by days_to_treatment_start)
  - residual_tumor (Liu) ->
      residual_disease (modern GDC)    diagnoses[primary].residual_disease
  - margin_status                      diagnoses[primary].treatments[].margin_status
                                            (lives on TREATMENT_FIELDS, not DIAGNOSIS_FIELDS)
  - Stage IV check                     diagnoses[primary].ajcc_pathologic_stage starts "Stage IV"

The "primary diagnosis" is the one with `diagnosis_is_primary_disease == True`,
falling back to the earliest by `days_to_diagnosis` if none is flagged.

## DFI tumor-type special cases

  - no DFI available:          SKCM, THYM, UVM, LAML
  - SARC: residual_disease only (Liu chose it over margin_status when both
                                  populated; clinically the right field for
                                  end-of-first-course assessment)
  - all other tumor types: disease-free is the OR of `treatment_outcome ==
    "Complete Response"`, `residual_disease == "R0"`, or
    `margin_status == "Uninvolved"`. Liu's STAR Methods is explicit about
    this being an OR — three independent signals, any one suffices. Earlier
    drafts of this module implemented a per-tumor-type chain instead, which
    over-excluded patients with a positive signal in a non-primary field.
"""

from __future__ import annotations

from typing import Any

# Tumor types where Liu defined no DFI. SKCM/THYM/UVM lack any usable
# disease-free signal in the three fields; LAML is a liquid tumor where
# "disease-free interval" doesn't apply at all (CDR records DFI=NA for
# 100% of LAML cases). Matched against the TCGA project_id suffix.
_DFI_NO_FIELD: set[str] = {"LAML", "SKCM", "THYM", "UVM"}

# SARC is the one tumor type where both `residual_disease` and `margin_status`
# are populated; Liu's STAR Methods says they chose `residual_tumor` (modern
# `residual_disease`) for SARC because the values are highly consistent and
# clinically that field reflects end-of-first-course state. For every other
# tumor type, all three disease-free signals are pooled with OR.
_DFI_RESIDUAL_ONLY: set[str] = {"SARC"}


# ---------------------------------------------------------------------------
# Field extractors
# ---------------------------------------------------------------------------


def _short_project(project_id: str | None) -> str | None:
    """`TCGA-CHOL` -> `CHOL`. Liu's per-type tables key on the suffix."""
    if not project_id or "-" not in project_id:
        return project_id
    return project_id.split("-", 1)[1]


def _primary_diagnosis(case: dict[str, Any]) -> dict[str, Any] | None:
    """Return the diagnosis flagged `diagnosis_is_primary_disease`, else the earliest.

    Some TCGA cases have multiple diagnosis records (prior primary, recurrence,
    etc.); Liu's algorithm always anchors on the *primary* one. We trust the
    GDC's own flag first and fall back to "earliest by days_to_diagnosis"
    only when no diagnosis is marked primary.
    """
    diagnoses = case.get("diagnoses") or []
    for dx in diagnoses:
        if dx.get("diagnosis_is_primary_disease"):
            return dx
    # `clinical._patient_row` already sorted diagnoses by days_to_diagnosis,
    # so [0] is the earliest fallback.
    return diagnoses[0] if diagnoses else None


def _extract_last_contact_days(case: dict[str, Any]) -> int | float | None:
    """Liu's `last_contact_days_to` — max across primary dx and follow_ups."""
    candidates: list[int | float] = []
    primary = _primary_diagnosis(case)
    if primary is not None:
        v = primary.get("days_to_last_follow_up")
        if v is not None:
            candidates.append(v)
    for fu in case.get("follow_ups") or []:
        v = fu.get("days_to_follow_up")
        if v is not None:
            candidates.append(v)
    return max(candidates) if candidates else None


def _extract_tumor_status(case: dict[str, Any]) -> str | None:
    """Latest tumor-status signal across follow_ups + diagnoses, mapped to Liu's binary.

    Returns "WITH TUMOR" / "TUMOR FREE" / None. Three rules:

    1. Group `disease_response` annotations by `days_to_follow_up`. Among
       entries at the same day, prefer WT over TF (some cases like
       TCGA-OR-A5JU record both at the same timestamp).
    2. Track the latest day a new tumor event was recorded (recurrence /
       progression follow-up, or non-primary diagnosis classified as
       recurrence / metastasis / subsequent primary).
    3. Pick the *latest* of the two timestamps. If a new tumor event is
       at-or-after the latest disease_response, return WITH TUMOR; else
       return whatever disease_response says at its latest timestamp.

    This handles two competing patterns:
      - TCGA-3X-AAV9 / TCGA-OR-A5JU: cancer present at last contact, WT
        signal latest -> WITH TUMOR.
      - TCGA-AO-A03U: had recurrence early but a later TF follow-up means
        the patient was clinically tumor-free at death -> TUMOR FREE.
    """
    by_day: dict[int | float, set[str]] = {}
    for fu in case.get("follow_ups") or []:
        raw = fu.get("disease_response")
        day = fu.get("days_to_follow_up")
        if raw is None or day is None:
            continue
        s = str(raw)
        if "With Tumor" in s or s.upper().startswith("WT"):
            by_day.setdefault(day, set()).add("WITH TUMOR")
        elif "Tumor Free" in s or s.upper().startswith("TF"):
            by_day.setdefault(day, set()).add("TUMOR FREE")

    new_tumor_day = _extract_new_tumor_event_days(case)
    if not by_day and new_tumor_day is None:
        return None
    if not by_day:
        return "WITH TUMOR"
    latest_dr = max(by_day)
    if new_tumor_day is not None and new_tumor_day >= latest_dr:
        return "WITH TUMOR"
    statuses = by_day[latest_dr]
    if "WITH TUMOR" in statuses:
        return "WITH TUMOR"
    if "TUMOR FREE" in statuses:
        return "TUMOR FREE"
    return None


# Values of `diagnoses[].classification_of_tumor` that mark a *secondary*
# diagnosis as a new tumor event (recurrence, metastasis, or subsequent
# primary). Excluded values: "primary" (the original tumor) and "Prior
# primary"-style entries (cancer history that pre-dates the index date).
# Modern GDC casing varies, so we lower-case before matching.
_NEW_TUMOR_DIAGNOSIS_KINDS: set[str] = {
    "recurrence",
    "metastatic",
    "subsequent primary",
    "additional - new primary",
    "additional new primary",
}


def _follow_up_event_days(fu: dict[str, Any]) -> int | float | None:
    """Days-from-index of a new tumor event recorded on this follow-up, if any.

    Three sources, in order of priority (lowest day wins if multiple are
    populated): explicit `days_to_recurrence`, explicit `days_to_progression`,
    or — for cases where the GDC recorded `progression_or_recurrence == 'Yes'`
    without a dedicated days_to_* — the follow-up's own `days_to_follow_up`.
    The last branch covers patients like TCGA-UY-A78L where Liu's CDR
    correctly identified an event but the explicit recurrence timestamps
    are missing in modern GDC.
    """
    candidates: list[int | float] = []
    for key in ("days_to_recurrence", "days_to_progression"):
        v = fu.get(key)
        if v is not None:
            candidates.append(v)
    por = fu.get("progression_or_recurrence")
    if por and str(por).strip().lower() == "yes":
        d = fu.get("days_to_follow_up")
        if d is not None:
            candidates.append(d)
    return min(candidates) if candidates else None


def _diagnosis_event_days(dx: dict[str, Any]) -> int | float | None:
    """Days-from-index of a non-primary diagnosis recorded as a new tumor event.

    Modern GDC encodes some recurrences / metastases / subsequent primaries
    as additional diagnosis rows with `classification_of_tumor` != "primary"
    and a `days_to_diagnosis` of when the new tumor was found. This is a
    second source of `new_tumor_event_dx_days_to` beyond the follow_ups
    (the only source the user's mapping notes called out).
    """
    if dx.get("diagnosis_is_primary_disease"):
        return None
    kind = (dx.get("classification_of_tumor") or "").strip().lower()
    if kind not in _NEW_TUMOR_DIAGNOSIS_KINDS:
        return None
    days = dx.get("days_to_diagnosis")
    return days if days is not None else None


def _extract_new_tumor_event_days(case: dict[str, Any]) -> int | float | None:
    """Liu's `new_tumor_event_dx_days_to` — earliest event across all sources.

    Pools follow-up timestamps (days_to_recurrence/progression) and
    diagnosis timestamps (classification_of_tumor in the recurrence /
    metastatic / subsequent-primary set), and returns the minimum.
    """
    times: list[int | float] = []
    for fu in case.get("follow_ups") or []:
        t = _follow_up_event_days(fu)
        if t is not None:
            times.append(t)
    for dx in case.get("diagnoses") or []:
        t = _diagnosis_event_days(dx)
        if t is not None:
            times.append(t)
    return min(times) if times else None


def _extract_new_tumor_event_type(case: dict[str, Any]) -> str | None:
    """Type of the *earliest* new tumor event across follow_ups and diagnoses.

    Liu's DFI special-cases "New Primary in Other Organ" as a censoring
    event rather than a DFI event. The type comes from whichever source
    supplied the earliest timestamp — `progression_or_recurrence_type` for
    follow-up entries, `classification_of_tumor` for diagnosis entries
    (the latter normalized to title case for cross-source comparability).
    """
    best_type: str | None = None
    best_time: int | float | None = None
    for fu in case.get("follow_ups") or []:
        t = _follow_up_event_days(fu)
        if t is None:
            continue
        if best_time is None or t < best_time:
            best_time = t
            best_type = fu.get("progression_or_recurrence_type")
    primary = _primary_diagnosis(case)
    primary_tissue = (primary or {}).get("tissue_or_organ_of_origin") if primary else None
    for dx in case.get("diagnoses") or []:
        t = _diagnosis_event_days(dx)
        if t is None:
            continue
        if best_time is not None and t >= best_time:
            continue
        best_time = t
        kind = (dx.get("classification_of_tumor") or "").strip().lower()
        # A Subsequent Primary in a *different* organ from the primary
        # tumor is Liu's "new primary in other organ" -> censored. In the
        # *same* organ it's effectively a same-disease recurrence and CDR
        # treats it as a DFI event. Map only the different-organ case to
        # the Liu-style "New Primary" sentinel so the DFI censor below
        # fires correctly.
        is_subsequent = kind in (
            "subsequent primary",
            "additional - new primary",
            "additional new primary",
        )
        dx_tissue = dx.get("tissue_or_organ_of_origin")
        if is_subsequent and primary_tissue and dx_tissue and dx_tissue != primary_tissue:
            best_type = "New Primary"
        else:
            best_type = dx.get("classification_of_tumor")
    return best_type


def _extract_treatment_outcome_first_course(case: dict[str, Any]) -> str | None:
    """Liu's `treatment_outcome_first_course` — first populated value.

    Used by `_is_disease_free` for DFI plus by code paths that just want
    to surface "what's the recorded outcome". Returns the first populated
    value across (clinical supplement patient form, supplement follow-ups
    in order, harmonized treatments[]). For the disease-free *check* used
    by DFI, prefer `_has_disease_free_outcome` — it inspects every
    populated value rather than collapsing to one, so a patient with
    "Complete Remission/Response" on one follow-up isn't rejected because
    a later follow-up overwrote the slot with a non-CR value (TCGA-2G-AAGA
    is the canonical example — chemo CR on early treatments, "No
    Measureable Tumor or Tumor Markers" on a later follow-up).
    """
    supp = case.get("clinical_supplement")
    if supp:
        from ds2hf_pipeline.tcga.clinical_supplement import first_value

        val = first_value(supp, "treatment_outcome_first_course")
        if val:
            return val

    primary = _primary_diagnosis(case)
    if primary is None:
        return None
    for tx in primary.get("treatments") or []:
        outcome = tx.get("treatment_outcome")
        if outcome:
            return str(outcome)
    return None


def _has_disease_free_outcome(case: dict[str, Any]) -> bool:
    """True if any recorded `treatment_outcome_first_course` says disease-free.

    Liu's algorithm flags a patient as disease-free if `treatment_outcome_first_course`
    ever equals "Complete Remission/Response" (modernized GDC also uses
    "Complete Response"). The BCR biotab can carry multiple follow-up forms
    per patient, and we count *any* form recording a CR/CRR signal — same
    semantics Liu's CDR ends up with after consolidation.
    """
    supp = case.get("clinical_supplement")
    if supp:
        from ds2hf_pipeline.tcga.clinical_supplement import any_disease_free_signal

        if any_disease_free_signal(supp):
            return True

    primary = _primary_diagnosis(case)
    if primary is None:
        return False
    for tx in primary.get("treatments") or []:
        outcome = tx.get("treatment_outcome")
        if outcome in ("Complete Response", "Complete Remission/Response"):
            return True
    return False


def _extract_residual_disease(case: dict[str, Any]) -> str | None:
    """Liu's `residual_tumor` — renamed `residual_disease` in modern GDC."""
    primary = _primary_diagnosis(case)
    if primary is None:
        return None
    rd = primary.get("residual_disease")
    return str(rd) if rd else None


def _extract_margin_status(case: dict[str, Any]) -> str | None:
    """First populated `margin_status` across the primary diagnosis's treatments."""
    primary = _primary_diagnosis(case)
    if primary is None:
        return None
    for tx in primary.get("treatments") or []:
        ms = tx.get("margin_status")
        if ms:
            return str(ms)
    return None


def _is_stage_iv(case: dict[str, Any]) -> bool:
    """True if AJCC pathologic stage starts with 'Stage IV'.

    GDC values are 'Stage I', 'Stage IA', ..., 'Stage IV', 'Stage IVA', etc.
    Substring 'Stage IV' won't accidentally match 'Stage III' (different
    Roman numerals), so a lowercase prefix check is safe.
    """
    primary = _primary_diagnosis(case)
    if primary is None:
        return False
    stage = primary.get("ajcc_pathologic_stage")
    if not stage:
        return False
    return str(stage).strip().lower().startswith("stage iv")


def _is_disease_free(case: dict[str, Any], project_short: str | None) -> bool:
    """Did the patient achieve a disease-free state at end of first course?

    Liu's STAR Methods specifies an *OR* across three independent signals
    (modernized GDC vocabulary in parens):

      - `treatment_outcome_first_course == "Complete Remission/Response"`
        (modern: `treatments[].treatment_outcome == "Complete Response"`)
      - `residual_tumor == "R0"`
        (modern: `diagnoses[primary].residual_disease == "R0"`)
      - `margin_status == "negative"`
        (modern: `treatments[].margin_status == "Uninvolved"`)

    Any one positive signal qualifies. SARC is the documented exception —
    Liu chose `residual_tumor` only for SARC because both fields were
    populated and clinically that one reflects end-of-first-course state.
    Patients with no positive signal in any field get DFI=NA.
    """
    rd = _extract_residual_disease(case)
    if project_short in _DFI_RESIDUAL_ONLY:
        return rd == "R0"
    ms = _extract_margin_status(case)
    return rd == "R0" or ms == "Uninvolved" or _has_disease_free_outcome(case)


# ---------------------------------------------------------------------------
# Endpoint derivations — one function per endpoint, returns (event, time)
# ---------------------------------------------------------------------------


def derive_os(case: dict[str, Any]) -> tuple[int | None, float | None]:
    """Liu's OS: 1 if Dead, 0 if Alive; time = max(last_contact, death)."""
    demo = case.get("demographic") or {}
    vital = demo.get("vital_status")
    if vital not in ("Alive", "Dead"):
        return (None, None)
    event = 1 if vital == "Dead" else 0
    death_days = demo.get("days_to_death")
    last_contact = _extract_last_contact_days(case)
    candidates = [v for v in (death_days, last_contact) if v is not None]
    if not candidates:
        return (None, None)
    return (event, float(max(candidates)))


def derive_dss(case: dict[str, Any]) -> tuple[int | None, float | None]:
    """Liu's DSS: alive=0, dead-with-tumor or cancer-cause=1, dead-tumor-free=0, else NA."""
    demo = case.get("demographic") or {}
    vital = demo.get("vital_status")
    if vital not in ("Alive", "Dead"):
        return (None, None)

    death_days = demo.get("days_to_death")
    last_contact = _extract_last_contact_days(case)
    time_candidates = [v for v in (death_days, last_contact) if v is not None]
    if not time_candidates:
        return (None, None)
    time = float(max(time_candidates))

    if vital == "Alive":
        return (0, time)

    # Dead — distinguish events by tumor_status with cause_of_death as
    # tiebreaker, in the order Liu specified ("Dead AND TUMOR FREE -> 0"
    # comes before "cause_of_death indicates the cancer -> 1"). When
    # neither field disambiguates, default to event=1: Liu's CDR
    # populates DSS=1 for dead patients without a TUMOR FREE signal,
    # which is also the conservative survival-analysis default (assume
    # cancer death until proven otherwise).
    tumor_status = _extract_tumor_status(case)
    if tumor_status == "TUMOR FREE":
        return (0, time)
    if tumor_status == "WITH TUMOR":
        return (1, time)
    if demo.get("cause_of_death") == "Cancer Related":
        return (1, time)
    return (1, time)


def derive_pfi(case: dict[str, Any]) -> tuple[int | None, float | None]:
    """Liu's PFI: any new tumor event, or died-with-cancer-no-event."""
    demo = case.get("demographic") or {}
    vital = demo.get("vital_status")
    if vital not in ("Alive", "Dead"):
        return (None, None)

    death_days = demo.get("days_to_death")
    last_contact = _extract_last_contact_days(case)
    new_tumor_days = _extract_new_tumor_event_days(case)

    if new_tumor_days is not None:
        return (1, float(new_tumor_days))

    # No new tumor event recorded — for Dead patients, check if they died
    # with the cancer (Liu's "died with cancer without new tumor event"
    # branch counts as PFI event, anchored at death).
    if vital == "Dead":
        tumor_status = _extract_tumor_status(case)
        cause_of_death = demo.get("cause_of_death")
        died_with_cancer = tumor_status == "WITH TUMOR" or cause_of_death == "Cancer Related"
        if died_with_cancer and death_days is not None:
            return (1, float(death_days))
        # Dead but tumor-free / unrelated cause -> censored at death (or last contact).
        time_candidates = [v for v in (death_days, last_contact) if v is not None]
        if not time_candidates:
            return (None, None)
        return (0, float(max(time_candidates)))

    # Alive, no new tumor event -> censored at last contact.
    if last_contact is None:
        return (None, None)
    return (0, float(last_contact))


def derive_dfi(case: dict[str, Any], project_id: str | None) -> tuple[int | None, float | None]:
    """Liu's DFI: any new tumor event, restricted to patients who started disease-free.

    Returns (None, None) when DFI is undefined for this patient — three reasons:
    tumor type lacks any usable disease-free signal (SKCM/THYM/UVM); patient
    wasn't disease-free at end of first course; or patient is Stage IV /
    dead-with-tumor-no-event (both Liu exclusions).
    """
    project_short = _short_project(project_id)

    if project_short in _DFI_NO_FIELD:
        return (None, None)
    if _is_stage_iv(case):
        return (None, None)

    demo = case.get("demographic") or {}
    vital = demo.get("vital_status")
    if vital not in ("Alive", "Dead"):
        return (None, None)

    new_tumor_days = _extract_new_tumor_event_days(case)

    # Exclusion: dead with tumor and no new tumor event recorded.
    if vital == "Dead" and new_tumor_days is None:
        tumor_status = _extract_tumor_status(case)
        cause_of_death = demo.get("cause_of_death")
        died_with_cancer = tumor_status == "WITH TUMOR" or cause_of_death == "Cancer Related"
        if died_with_cancer:
            return (None, None)

    # Liu's "disease-free at end of first course" check (STAR Methods,
    # OR across the three signals; see `_is_disease_free`). Patients with
    # no positive signal in any of `residual_disease="R0"`,
    # `margin_status="Uninvolved"`, or `treatment_outcome="Complete
    # Response"` get DFI=NA — Liu calls them "never disease-free".
    if not _is_disease_free(case, project_short):
        return (None, None)

    # New primary in other organ -> censored, not an event. Modern GDC
    # encodes this only via follow-up `progression_or_recurrence_type`
    # (with values like "New Primary"); Liu's CDR doesn't censor diagnosis-
    # level "Subsequent Primary" rows even though the words are similar,
    # so we restrict the censor to the follow-up signal alone.
    if new_tumor_days is not None:
        nte_type = _extract_new_tumor_event_type(case)
        if nte_type and "New Primary" in nte_type and "+" not in nte_type:
            new_tumor_days = None

    if new_tumor_days is not None:
        return (1, float(new_tumor_days))

    # Censored — time is last_contact_days_to or death_days_to (whichever applies).
    death_days = demo.get("days_to_death")
    last_contact = _extract_last_contact_days(case)
    time_candidates = [v for v in (death_days, last_contact) if v is not None]
    if not time_candidates:
        return (None, None)
    return (0, float(max(time_candidates)))


# ---------------------------------------------------------------------------
# Pipeline integration
# ---------------------------------------------------------------------------


def attach_survival(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Mutate patient rows to populate the `survival_derived` struct column.

    Sets `row["survival_derived"]` to a dict with the 8 re-derived
    endpoint values. The struct shape matches `ds2hf.tcga.schema.SURVIVAL_DERIVED_FIELDS`,
    so pyarrow's `Table.from_pylist(rows, schema=PATIENTS)` lifts these
    into the struct column directly.
    """
    for row in rows:
        os_event, os_time = derive_os(row)
        dss_event, dss_time = derive_dss(row)
        pfi_event, pfi_time = derive_pfi(row)
        dfi_event, dfi_time = derive_dfi(row, row.get("project_id"))
        row["survival_derived"] = {
            "os_event": os_event,
            "os_time": os_time,
            "dss_event": dss_event,
            "dss_time": dss_time,
            "pfi_event": pfi_event,
            "pfi_time": pfi_time,
            "dfi_event": dfi_event,
            "dfi_time": dfi_time,
        }
    return rows
