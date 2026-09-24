## Section 3 — Liu's algorithm definitions and our implementation

All four endpoints anchor on **date of diagnosis = time zero** (Liu's choice; see *Choice of time zero* in the STAR Methods). For each event the time is days from diagnosis to either the event or the censoring contact.

### Verbatim from Liu's TCGA-CDR_Notes sheet

**OS** — *overall survival event, 1 for death from any cause, 0 for alive.* OS.time = `last_contact_days_to` or `death_days_to`, whichever is larger.

**DSS** — *disease-specific survival event, 1 for patient whose vital_status was Dead and tumor_status was WITH TUMOR. If a patient died from the disease shown in field of cause_of_death, the status of DSS would be 1 for the patient. 0 for patient whose vital_status was Alive or whose vital_status was Dead and tumor_status was TUMOR FREE. This is not a 100% accurate definition but is the best we could do with this dataset.* DSS.time = same as OS.time.

**PFI** — *progression-free interval event, 1 for patient having new tumor event whether it was a progression of disease, local recurrence, distant metastasis, new primary tumors all sites, or died with the cancer without new tumor event, including cases with a new tumor event whose type is N/A.* PFI.time = `new_tumor_event_dx_days_to` or `death_days_to` for events; `last_contact_days_to` or `death_days_to` for censored.

**DFI** — *disease-free interval event, 1 for patient having new tumor event whether it is a local recurrence, distant metastasis, new primary tumor of the cancer, including cases with a new tumor event whose type is N/A. Disease free was defined by:* `treatment_outcome_first_course == "Complete Remission/Response"` *OR* `residual_tumor == "R0"` *OR* `margin_status == "negative"`. *New primary tumor in other organ was censored; patients who were Dead with tumor without new tumor event are excluded; patients with stage IV are excluded too.*

### Field rename in modern GDC

Liu's old-TCGA field names → modern GDC paths:

| Liu | Modern GDC path |
|---|---|
| `vital_status` | `demographic.vital_status` |
| `death_days_to` | `demographic.days_to_death` |
| `cause_of_death` | `demographic.cause_of_death` |
| `last_contact_days_to` | `max(diagnoses[primary].days_to_last_follow_up, max(follow_ups[].days_to_follow_up))` |
| `tumor_status` | latest `follow_ups[].disease_response` (`WT-With Tumor` → `WITH TUMOR`, `TF-Tumor Free` → `TUMOR FREE`) |
| `new_tumor_event_dx_days_to` | min populated `follow_ups[].days_to_recurrence` / `.days_to_progression` / `progression_or_recurrence==Yes` follow-up day |
| `treatment_outcome_first_course` | **BCR biotab Clinical Supplement** `clinical_patient_<proj>.txt` + `clinical_follow_up_v*_<proj>.txt`, field `treatment_outcome_first_course`. Falls back to harmonized `diagnoses[primary].treatments[].treatment_outcome` when supplement unavailable. |
| `residual_tumor` (Liu) → `residual_disease` (modern GDC) | `diagnoses[primary].residual_disease` *(renamed)* |
| `margin_status` | `diagnoses[primary].treatments[].margin_status` |

### About the BCR biotab Clinical Supplements

The harmonized `/cases?expand=...` API drops or under-populates `treatment_outcome_first_course` (Liu's `primary_therapy_outcome_success`). We fetch the original BCR biotab files separately (one TSV per project per form: patient, follow_up, nte, drug, radiation, ablation, omf) and read the field directly. For the same patient TCGA-DK-A2I6: the harmonized API has `treatment_outcome=None` while the BCR biotab follow-up form has the canonical `"Complete Remission/Response"` Liu used. Adding the supplement integration moved DFI's overall match rate against Liu's CDR from 34.6% to 48.5%, dropping the under-population gap (`cdr_pop_der_na`) from 1,625 patients to 52 — a 97% reduction. The biotab data is also shipped as 7 new tables in the [`gabrielaltay/tcga-tabular-open`][tabular] HF dataset (per-project schemas).

The disease-free check uses **"any form ever recorded CR/CRR"** semantics, not "latest value wins". TCGA-2G-AAGA is the canonical example: CR on early chemo-treatment follow-ups, then "No Measureable Tumor or Tumor Markers" on a later follow-up. The any-form check correctly classifies as disease-free (matching Liu's CDR); the naive latest-value-wins approach would reject this patient.

[tabular]: https://huggingface.co/datasets/gabrielaltay/tcga-tabular-open

### Liu's special cases (STAR Methods, *Handling of special cases and problems*)

From the paper's *Handling of special cases and problems in clinical data files* (page e2). Liu reports having resolved "over 1000 cases with apparent or real problems":

1. **7 Stage 0 SKCM cases** — actually distant metastases, kept.
2. **483 Dead-with-Tumor cases without a defined NTE** — excluded from DFI; kept in OS, PFI, DSS.
3. **62 Dead-with-Tumor cases that *did* have a defined NTE** — censored for DSS at the NTE date (vital/tumor inconsistency).
4. **10 Dead-with-Tumor-Free cases with cancer-related cause-of-death** — DSS event resolved via cause_of_death.
5. **797 of 3,346 NTEs without a specified type** — kept as PFI events under the "relaxed" PFI definition.
6. **6 cases with negative `last_contact_days_to`, 6 with negative `new_tumor_event_dx_days_to`** — clamped to 0.
7. **46 patients aged 90** — HIPAA cap, ages are artificial.
8. **Multiple follow-up files** — file with the longer follow-up wins.
9. **One Alive-then-Dead patient with no death date** — used the enrollment-file Alive value.
10. **One OV with grade 4** — invalid grade for OV but kept as-is.

Our `ds2hf_pipeline.tcga.survival` follows Liu's algorithm verbatim where the modern GDC schema permits. The `cdr_*` columns are Liu's verbatim values (subject to Liu's special-case handling); the re-derived `_event`/`_time` columns use our implementation against current data plus the Clinical Supplement augmentation.

**Stage IV exclusion (DFI):** Liu excluded 1,095 stage IV patients from DFI. We do the same via `_is_stage_iv` checking `ajcc_pathologic_stage` starts with `"Stage IV"`.

### DFI tumor-type special cases

- **No DFI available**: SKCM, THYM, UVM, LAML — Liu's `_DFI_NO_FIELD` set. These tumor types lack any usable disease-free signal in the three fields (or in LAML's case, are liquid tumors where "disease-free interval" doesn't apply). Our `derive_dfi` returns `(None, None)` for these.
- **SARC**: residual_disease only (Liu chose it over margin_status when both populated; clinically the right field for end-of-first-course assessment).
- **All other tumor types**: disease-free is the OR of `treatment_outcome == "Complete Response"`, `residual_disease == "R0"`, or `margin_status == "Uninvolved"`. Liu's STAR Methods is explicit about this being an OR — three independent signals, any one suffices. Earlier drafts of our `survival.py` implemented a per-tumor-type chain instead, which over-excluded patients with a positive signal in a non-primary field.
