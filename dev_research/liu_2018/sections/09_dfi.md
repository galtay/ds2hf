## Section 9 — DFI deep-dive (the hardest endpoint)

DFI is structurally the hardest endpoint. Liu themselves flag it as the most difficult to derive and recommend it for only 27 of 33 cancer types. Three structural complications:

- **Three-way OR for disease-free**: `treatment_outcome_first_course == "Complete Remission/Response"` OR `residual_disease == "R0"` OR `margin_status == "Uninvolved"`. Patients with no positive signal in any of the three are excluded from DFI (Liu calls them "never disease-free").
- **Tumor-type exclusions**: SKCM, THYM, UVM, LAML have no usable disease-free signal in any TCGA sample. Liu excludes them entirely; we follow.
- **SARC override**: Liu used `residual_tumor` only for SARC (both fields populated; clinically the right one).
- **Stage IV excluded** (1,095 patients in Liu's count).
- **Dead-with-tumor-no-NTE excluded** (Liu's special case #2, 483 patients).
- **New primary in *other organ* censored**, not an event.

### How the BCR biotab fix moved the numbers

Adding the BCR biotab Clinical Supplement integration (see Section 3) was the largest single improvement to DFI:

| metric | before supplement | after supplement |
|---|---|---|
| **Agreement rate** (both NA, OR both populated and agreeing within 30 days) | 77.2% (8,620/11,160) | **90.1% (10,055/11,160)** |
| Under-population (Liu had it, we missed) | 1,625 | **52** (97% reduction) |
| Direction agreement when both populated | 99.7% | **99.7%** (unchanged — the gain is on coverage) |

The supplement fix recovered most of the under-population: Liu's 2018 algorithm read `treatment_outcome_first_course` from the BCR-original forms, but the modern GDC harmonized API drops or under-populates that field (e.g. `treatments[].treatment_outcome=None` for patients with `Complete Remission/Response` recorded on a follow-up form). We now read the BCR biotab directly.

### Five-bucket breakdown (full cohort)

Each of the 11,160 CDR-matched patients lands in one of five buckets when comparing Liu's `cdr_DFI`/`cdr_DFI_time` against our re-derived `dfi_event`/`dfi_time`:

| bucket | patients | % | meaning |
| --- | --- | --- | --- |
| match | 5412 | 48.5 | Both populated, event AND time agree (within 0.5 days) |
| both_na | 4640 | 41.6 | Both correctly say NA (Liu's documented exclusions) |
| der_pop_cdr_na | 989 | 8.9 | We populated, Liu didn't |
| mismatch | 65 | 0.6 | Both populated, disagree on event or time |
| cdr_pop_der_na | 54 | 0.5 | Liu populated, we don't |

### How to read the match rates

Several useful framings of "how often do we and Liu agree":

| framing | math | rate |
| --- | --- | --- |
| Agreement rate (both NA OR both populated within 30d) | 10054 / 11160 | 90.1% |
| Match where Liu populated | 5412 / 5531 | 97.8% |
| Match where both populated (exact event + time) | 5412 / 5477 | 98.8% |
| Event-direction agreement where both populated | 5461 / 5477 | 99.7% |

The headline metric is the **agreement rate** (90.1%): patients where we and Liu reach the same conclusion, defined as either (a) both correctly NA, or (b) both populated, event direction agrees, time within 30 days. This treats Liu's documented exclusions (SKCM/THYM/UVM/LAML cohorts; stage IV; never-disease-free) as correct agreement when we honor them — which we do — rather than penalizing the score for them.

When Liu had a value, we agree exactly **97.8%** of the time. When both have a value, the event direction (had-recurrence vs censored) agrees **99.7%** of the time — disagreements are almost always about *when* the event happened, not *whether* it happened.

### Per-project DFI populated counts

How many DFI patients each project has, comparing Liu's curated count vs ours:

| project | N | Liu_pop | ours_pop | delta |
| --- | --- | --- | --- | --- |
| BRCA | 1097 | 953 | 950 | -3 |
| UCEC | 548 | 426 | 474 | 48 |
| THCA | 507 | 358 | 390 | 32 |
| PRAD | 500 | 340 | 423 | 83 |
| LIHC | 377 | 323 | 322 | -1 |
| LUAD | 522 | 309 | 398 | 89 |
| LUSC | 504 | 304 | 435 | 131 |
| OV | 587 | 286 | 338 | 52 |
| STAD | 443 | 259 | 336 | 77 |
| COAD | 459 | 190 | 327 | 137 |
| BLCA | 412 | 189 | 190 | 1 |
| KIRP | 291 | 184 | 184 | 0 |
| CESC | 307 | 176 | 196 | 20 |
| PCPG | 179 | 160 | 160 | 0 |
| SARC | 261 | 153 | 153 | 0 |
| LGG | 515 | 134 | 139 | 5 |
| HNSC | 528 | 134 | 193 | 59 |
| KIRC | 537 | 117 | 117 | 0 |
| TGCT | 134 | 105 | 96 | -9 |
| ESCA | 185 | 89 | 143 | 54 |
| KICH | 113 | 71 | 71 | 0 |
| PAAD | 185 | 69 | 111 | 42 |
| ACC | 92 | 53 | 69 | 16 |
| READ | 170 | 48 | 127 | 79 |
| DLBC | 48 | 28 | 37 | 9 |
| CHOL | 45 | 28 | 28 | 0 |
| UCS | 57 | 27 | 37 | 10 |
| MESO | 87 | 15 | 17 | 2 |
| GBM | 596 | 3 | 5 | 2 |
| SKCM | 470 | 0 | 0 | 0 |
| LAML | 200 | 0 | 0 | 0 |
| THYM | 124 | 0 | 0 | 0 |
| UVM | 80 | 0 | 0 | 0 |

### What's causing the remaining differences

The five buckets break down structurally as follows:

#### `cdr_pop_der_na` (54 patients) — Liu had a value, we still don't

Concentrated in BRCA (25) and TGCT (15). Liu's 2018 BCR data had a disease-free signal for these patients that's been re-curated or removed in the modern data — both the harmonized API and the BCR biotab (the canonical source) currently show no signal. These are likely irrecoverable without a 2018 GDC snapshot.

#### `mismatch` (65 patients) — Both populated, disagree

Splits into two sub-categories:

- **16 event-direction mismatches** (12 Liu=0/ours=1, 4 Liu=1/ours=0). The Liu=0/ours=1 cases are the same data drift signature as OS: patients censored in Liu's 2018 freeze who have since had a recurrence/progression entered. The reverse (Liu=1/ours=0, 4 cases) is rarer — probably patients where Liu's hand-fix step classified as event but current data doesn't, OR a tumor event that's been re-classified/removed in the GDC since 2018.
- **46 time-only mismatches** (event direction agrees, time differs ≥0.5 days). Median diff is **-447 days** (ours earlier than Liu). The pattern: many patients have `ours_time = 0` while Liu had a real value. Example: BRCA TCGA-D8-A1JK had Liu DFI_time = 612 days; current biotab has `last_contact_days_to = 0`. The 612 isn't anywhere in modern biotab — the 2018 BCR data had a follow-up timestamp that the GDC's re-curation since then has lost or replaced with 0.

#### `der_pop_cdr_na` (989 patients) — We populate, Liu didn't

Concentrated in COAD (137), LUSC (132), LUAD (92), PRAD (83). 628 are event=0 (censored — we compute, Liu didn't); 361 are event=1 (we found a recurrence Liu couldn't classify because she lacked a baseline disease-free signal). The biotab gives us `treatment_outcome_first_course` data Liu didn't have in 2018 — for these patients, Liu's algorithm marked NA because no disease-free signal was found in the 2018 data; we now find one in the modernized BCR forms and can populate DFI. **This isn't a bug** — these are genuinely additional patients we can compute DFI for that Liu couldn't.

#### `both_na` (4640 patients) — Both correctly say NA

Liu's documented exclusions: SKCM/THYM/UVM/LAML cohorts entirely (no usable disease-free signal in any TCGA sample), plus stage IV patients, plus dead-with-tumor-no-NTE patients, plus "never disease-free" patients (no positive signal in any of `treatment_outcome / residual_disease / margin_status`). We honor all of these.

#### `match` (5412 patients) — Both populated, exact agreement

Event AND time agree (within 0.5 days).

### Implication

DFI is now usable as a re-derived endpoint for most projects. Where you need exact Liu reproduction, use `DFI` / `DFI.time` from `gabrielaltay/tcga-pancanatlas-cdr`. Where you want the broadest coverage including post-2018 patients, use `dfi_event` / `dfi_time` and accept the residual ~150 patient gap from Liu's own values.

Causal summary of the remaining ~1.1% disagreement when both populated:
- **About half is data drift** — modern GDC has re-curated, lost, or updated BCR timestamps that Liu had in 2018 (irrecoverable without a 2018 GDC snapshot).
- **About a quarter is algorithm-permissiveness** — our OR-of-three check classifies patients Liu's hand-curated step excluded.
- **About a quarter is post-2018 vital-status updates** — patients censored in Liu's freeze who've since had events recorded.
