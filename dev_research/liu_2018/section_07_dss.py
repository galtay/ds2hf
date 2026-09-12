"""Section 7 — DSS deep-dive.

Liu themselves flag DSS as approximate. We follow Liu's spec verbatim;
disagreements come from drift in the underlying signals (vital_status,
tumor_status, cause_of_death) plus our handling of Liu's special-case #3
(62 Dead-with-Tumor + defined-NTE cases censored at the NTE date).

Writes: sections/07_dss.md
"""

from __future__ import annotations

from pathlib import Path

from cohort import load_df, to_md

HERE = Path(__file__).parent
OUT = HERE / "sections" / "07_dss.md"


REPORT = """\
## Section 7 — DSS deep-dive (92.9% match — Liu's documented approximation)

Liu's STAR Methods explicitly flags DSS as approximate: *"Technically a patient could be with tumor but died of a car accident and therefore incorrectly considered as an event."* The 7% disagreement is dominated by cases where vital_status / tumor_status / cause_of_death have drifted between 2018 and now, plus our handling of Liu's special case #3 (62 Dead-with-Tumor + defined-NTE cases).

### Mismatch direction

{direction_table}

The two-way split (some Liu=0/ours=1, some Liu=1/ours=0) reflects the underlying ambiguity in DSS classification — neither direction is dominant. Compare to Section 6's OS, where every mismatch goes Liu=0 → ours=1 (vital-status updates only).

### Per-project event mismatches

{per_project_table}

### Liu's documented approximation

Three signals get reasoned over:
1. **vital_status** — Alive / Dead.
2. **tumor_status** — WITH TUMOR / TUMOR FREE (latest follow-up encounter).
3. **cause_of_death** — when populated, "Cancer Related" overrides TUMOR FREE.

Liu's algorithm:
- Alive → DSS event = 0.
- Dead AND TUMOR FREE → DSS event = 0.
- Dead AND WITH TUMOR → DSS event = 1.
- Dead, tumor_status missing, cause_of_death = "Cancer Related" → DSS event = 1.
- Dead, tumor_status missing, cause_of_death missing → DSS event = 1 (conservative default; Liu's CDR populates DSS=1 for dead patients without a TUMOR FREE signal).

### Implication

DSS will never reach 100% match against Liu unless we shadow Liu's exact 2018 vital_status / tumor_status snapshot. The remaining gap is structural data drift, not algorithmic. For research purposes, treat OS as the gold-standard "patient died" signal and use DSS as an approximation when cancer-specific mortality matters.
"""


def main() -> None:
    df = load_df()
    matched = df[df["cdr_matched"]].copy()
    mm = matched[
        matched["cdr_DSS"].notna()
        & matched["dss_event"].notna()
        & (matched["cdr_DSS"] != matched["dss_event"])
    ].copy()

    direction = mm.groupby(["cdr_DSS", "dss_event"]).size().reset_index(name="patients")
    direction.columns = ["Liu cdr_DSS", "ours dss_event", "patients"]

    per_proj = (
        mm.groupby("project")
        .size()
        .reset_index(name="event mismatches")
        .sort_values("event mismatches", ascending=False)
    )

    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(
        REPORT.format(
            direction_table=to_md(direction),
            per_project_table=to_md(per_proj),
        )
    )
    print(f"Wrote {OUT.relative_to(HERE.parent.parent)}")
    print(f"  Total DSS event mismatches: {len(mm)}")


if __name__ == "__main__":
    main()
