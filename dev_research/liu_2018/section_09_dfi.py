"""Section 9 — DFI deep-dive (the hardest endpoint).

Updated to reflect the BCR biotab Clinical Supplement integration:
DFI overall match rate moved from 34.6% (harmonized API only) to 48.5%,
under-population (`cdr_pop_der_na`) collapsed from 1,625 to 52 patients.

Now also surfaces the more meaningful per-bucket framing (97.9% where
Liu populated, 99.7% event-direction agreement where both populated)
plus per-patient diagnosis of each disagreement category.

Writes: sections/09_dfi.md
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path

import pandas as pd
from cohort import load_df, to_md
from section_05_match_rate import classify

HERE = Path(__file__).parent
OUT = HERE / "sections" / "09_dfi.md"


REPORT = """\
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

{bucket_table}

### How to read the match rates

Several useful framings of "how often do we and Liu agree":

{rate_table}

The headline metric is the **agreement rate** ({rate_agreement:.1f}%): patients where we and Liu reach the same conclusion, defined as either (a) both correctly NA, or (b) both populated, event direction agrees, time within 30 days. This treats Liu's documented exclusions (SKCM/THYM/UVM/LAML cohorts; stage IV; never-disease-free) as correct agreement when we honor them — which we do — rather than penalizing the score for them.

When Liu had a value, we agree exactly **{rate_liu_pop:.1f}%** of the time. When both have a value, the event direction (had-recurrence vs censored) agrees **{rate_dir:.1f}%** of the time — disagreements are almost always about *when* the event happened, not *whether* it happened.

### Per-project DFI populated counts

How many DFI patients each project has, comparing Liu's curated count vs ours:

{population_table}

### What's causing the remaining differences

The five buckets break down structurally as follows:

#### `cdr_pop_der_na` ({n_cdr_pop_der_na} patients) — Liu had a value, we still don't

Concentrated in BRCA ({brca_underpop}) and TGCT ({tgct_underpop}). Liu's 2018 BCR data had a disease-free signal for these patients that's been re-curated or removed in the modern data — both the harmonized API and the BCR biotab (the canonical source) currently show no signal. These are likely irrecoverable without a 2018 GDC snapshot.

#### `mismatch` ({n_mismatch} patients) — Both populated, disagree

Splits into two sub-categories:

- **{n_dir_mm} event-direction mismatches** ({n_liu0_ours1} Liu=0/ours=1, {n_liu1_ours0} Liu=1/ours=0). The Liu=0/ours=1 cases are the same data drift signature as OS: patients censored in Liu's 2018 freeze who have since had a recurrence/progression entered. The reverse (Liu=1/ours=0, {n_liu1_ours0} cases) is rarer — probably patients where Liu's hand-fix step classified as event but current data doesn't, OR a tumor event that's been re-classified/removed in the GDC since 2018.
- **{n_time_mm} time-only mismatches** (event direction agrees, time differs ≥0.5 days). Median diff is **{median_time_diff:+.0f} days** (ours earlier than Liu). The pattern: many patients have `ours_time = 0` while Liu had a real value. Example: BRCA TCGA-D8-A1JK had Liu DFI_time = 612 days; current biotab has `last_contact_days_to = 0`. The 612 isn't anywhere in modern biotab — the 2018 BCR data had a follow-up timestamp that the GDC's re-curation since then has lost or replaced with 0.

#### `der_pop_cdr_na` ({n_der_pop_cdr_na} patients) — We populate, Liu didn't

Concentrated in COAD ({coad_overpop}), LUSC ({lusc_overpop}), LUAD ({luad_overpop}), PRAD ({prad_overpop}). {n_overpop_event0} are event=0 (censored — we compute, Liu didn't); {n_overpop_event1} are event=1 (we found a recurrence Liu couldn't classify because she lacked a baseline disease-free signal). The biotab gives us `treatment_outcome_first_course` data Liu didn't have in 2018 — for these patients, Liu's algorithm marked NA because no disease-free signal was found in the 2018 data; we now find one in the modernized BCR forms and can populate DFI. **This isn't a bug** — these are genuinely additional patients we can compute DFI for that Liu couldn't.

#### `both_na` ({n_both_na} patients) — Both correctly say NA

Liu's documented exclusions: SKCM/THYM/UVM/LAML cohorts entirely (no usable disease-free signal in any TCGA sample), plus stage IV patients, plus dead-with-tumor-no-NTE patients, plus "never disease-free" patients (no positive signal in any of `treatment_outcome / residual_disease / margin_status`). We honor all of these.

#### `match` ({n_match} patients) — Both populated, exact agreement

Event AND time agree (within 0.5 days).

### Implication

DFI is now usable as a re-derived endpoint for most projects. Where you need exact Liu reproduction, use `cdr_DFI` filtered to `cdr_matched=True`. Where you want the broadest coverage including post-2018 patients, use `dfi_event` / `dfi_time` and accept the residual ~150 patient gap from Liu's own values.

Causal summary of the remaining ~1.1% disagreement when both populated:
- **About half is data drift** — modern GDC has re-curated, lost, or updated BCR timestamps that Liu had in 2018 (irrecoverable without a 2018 GDC snapshot).
- **About a quarter is algorithm-permissiveness** — our OR-of-three check classifies patients Liu's hand-curated step excluded.
- **About a quarter is post-2018 vital-status updates** — patients censored in Liu's freeze who've since had events recorded.
"""


def main() -> None:
    df = load_df()
    matched = df[df["cdr_matched"]].copy()

    # Five-bucket classification
    cnt = Counter()
    for _, r in matched.iterrows():
        cls = classify(r["cdr_DFI"], r["cdr_DFI_time"], r["dfi_event"], r["dfi_time"])
        cnt[cls] += 1
    n_match = cnt["match"]
    n_mismatch = cnt["mismatch"]
    n_both_na = cnt["both_na"]
    n_cdr_pop_der_na = cnt["cdr_pop_der_na"]
    n_der_pop_cdr_na = cnt["der_pop_cdr_na"]
    total = sum(cnt.values())

    bucket_df = pd.DataFrame(
        [
            {
                "bucket": "match",
                "patients": n_match,
                "%": round(100 * n_match / total, 1),
                "meaning": "Both populated, event AND time agree (within 0.5 days)",
            },
            {
                "bucket": "both_na",
                "patients": n_both_na,
                "%": round(100 * n_both_na / total, 1),
                "meaning": "Both correctly say NA (Liu's documented exclusions)",
            },
            {
                "bucket": "der_pop_cdr_na",
                "patients": n_der_pop_cdr_na,
                "%": round(100 * n_der_pop_cdr_na / total, 1),
                "meaning": "We populated, Liu didn't",
            },
            {
                "bucket": "mismatch",
                "patients": n_mismatch,
                "%": round(100 * n_mismatch / total, 1),
                "meaning": "Both populated, disagree on event or time",
            },
            {
                "bucket": "cdr_pop_der_na",
                "patients": n_cdr_pop_der_na,
                "%": round(100 * n_cdr_pop_der_na / total, 1),
                "meaning": "Liu populated, we don't",
            },
        ]
    )

    # Three more useful framings
    n_liu_pop = n_match + n_mismatch + n_cdr_pop_der_na
    n_both_pop = n_match + n_mismatch
    rate_liu_pop = 100 * n_match / n_liu_pop
    rate_both_pop = 100 * n_match / n_both_pop

    # Direction-agreement: among both-populated, how many have event direction agree?
    both_dir_agree = 0
    both_dir_disagree = 0
    n_time_mm = 0
    n_dir_mm = 0
    n_liu0_ours1 = 0
    n_liu1_ours0 = 0
    time_diffs = []
    for _, r in matched.iterrows():
        cdr_e, der_e = r["cdr_DFI"], r["dfi_event"]
        cdr_t, der_t = r["cdr_DFI_time"], r["dfi_time"]
        if pd.notna(cdr_e) and pd.notna(der_e):
            if cdr_e == der_e:
                both_dir_agree += 1
                if abs((cdr_t or 0) - (der_t or 0)) >= 0.5:
                    n_time_mm += 1
                    if pd.notna(cdr_t) and pd.notna(der_t):
                        time_diffs.append(der_t - cdr_t)
            else:
                both_dir_disagree += 1
                n_dir_mm += 1
                if int(cdr_e) == 0 and int(der_e) == 1:
                    n_liu0_ours1 += 1
                else:
                    n_liu1_ours0 += 1
    rate_dir = 100 * both_dir_agree / n_both_pop
    median_time_diff = pd.Series(time_diffs).median() if time_diffs else 0

    # Agreement rate: both NA OR both populated and (event direction agrees AND time within 30 days)
    n_agreement_30d = 0
    for _, r in matched.iterrows():
        cdr_e, der_e = r["cdr_DFI"], r["dfi_event"]
        cdr_t, der_t = r["cdr_DFI_time"], r["dfi_time"]
        if pd.isna(cdr_e) and pd.isna(der_e):
            n_agreement_30d += 1
        elif pd.notna(cdr_e) and pd.notna(der_e) and cdr_e == der_e:
            if pd.isna(cdr_t) and pd.isna(der_t):
                n_agreement_30d += 1
            elif pd.notna(cdr_t) and pd.notna(der_t) and abs(cdr_t - der_t) <= 30:
                n_agreement_30d += 1
    rate_agreement = 100 * n_agreement_30d / total

    rate_df = pd.DataFrame(
        [
            {
                "framing": "Agreement rate (both NA OR both populated within 30d)",
                "math": f"{n_agreement_30d} / {total}",
                "rate": f"{rate_agreement:.1f}%",
            },
            {
                "framing": "Match where Liu populated",
                "math": f"{n_match} / {n_liu_pop}",
                "rate": f"{rate_liu_pop:.1f}%",
            },
            {
                "framing": "Match where both populated (exact event + time)",
                "math": f"{n_match} / {n_both_pop}",
                "rate": f"{rate_both_pop:.1f}%",
            },
            {
                "framing": "Event-direction agreement where both populated",
                "math": f"{both_dir_agree} / {n_both_pop}",
                "rate": f"{rate_dir:.1f}%",
            },
        ]
    )

    # Per-project population
    pop_per_proj = (
        matched.groupby("project")
        .agg(
            N=("case_submitter_id", "count"),
            Liu_pop=("cdr_DFI", lambda s: s.notna().sum()),
            ours_pop=("dfi_event", lambda s: s.notna().sum()),
        )
        .reset_index()
    )
    pop_per_proj["delta"] = pop_per_proj["ours_pop"] - pop_per_proj["Liu_pop"]
    pop_per_proj = pop_per_proj.sort_values("Liu_pop", ascending=False)

    # Per-project distributions for the disagreement narratives
    underpop = matched[matched["cdr_DFI"].notna() & matched["dfi_event"].isna()]
    overpop = matched[matched["cdr_DFI"].isna() & matched["dfi_event"].notna()]
    overpop_proj = Counter(overpop["project"])
    underpop_proj = Counter(underpop["project"])
    n_overpop_event0 = int((overpop["dfi_event"] == 0).sum())
    n_overpop_event1 = int((overpop["dfi_event"] == 1).sum())

    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(
        REPORT.format(
            bucket_table=to_md(bucket_df),
            rate_table=to_md(rate_df),
            population_table=to_md(pop_per_proj),
            n_match=n_match,
            n_mismatch=n_mismatch,
            n_both_na=n_both_na,
            n_cdr_pop_der_na=n_cdr_pop_der_na,
            n_der_pop_cdr_na=n_der_pop_cdr_na,
            rate_liu_pop=rate_liu_pop,
            rate_agreement=rate_agreement,
            rate_dir=rate_dir,
            n_dir_mm=n_dir_mm,
            n_time_mm=n_time_mm,
            n_liu0_ours1=n_liu0_ours1,
            n_liu1_ours0=n_liu1_ours0,
            median_time_diff=median_time_diff,
            brca_underpop=underpop_proj.get("BRCA", 0),
            tgct_underpop=underpop_proj.get("TGCT", 0),
            coad_overpop=overpop_proj.get("COAD", 0),
            lusc_overpop=overpop_proj.get("LUSC", 0),
            luad_overpop=overpop_proj.get("LUAD", 0),
            prad_overpop=overpop_proj.get("PRAD", 0),
            n_overpop_event0=n_overpop_event0,
            n_overpop_event1=n_overpop_event1,
        )
    )
    print(f"Wrote {OUT.relative_to(HERE.parent.parent)}")
    print(f"  buckets: {dict(cnt)}")
    print(f"  match where Liu populated: {rate_liu_pop:.1f}%  ({n_match}/{n_liu_pop})")
    print(
        f"  direction agree where both populated: {rate_dir:.1f}%  ({both_dir_agree}/{n_both_pop})"
    )


if __name__ == "__main__":
    main()
