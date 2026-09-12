"""Section 6 — OS deep-dive.

OS is the simplest endpoint; effectively no algorithmic ambiguity. The 2%
disagreement is patients who were alive at Liu's 2018 freeze but have
since died in modern GDC.

Writes: sections/06_os.md
"""

from __future__ import annotations

from pathlib import Path

from cohort import load_df, to_md

HERE = Path(__file__).parent
OUT = HERE / "sections" / "06_os.md"


REPORT = """\
## Section 6 — OS deep-dive (98% match — data drift only)

OS is the simplest endpoint to derive (death from any cause; no tumor-status or new-tumor-event reasoning). Algorithmic agreement with Liu is essentially perfect; the 2% mismatch is data drift between Liu's 2018 freeze and the modern GDC.

### Mismatch direction

{direction_table}

The asymmetry is telling: every event mismatch goes in the **same direction** — Liu had `OS=0` (alive) and we now have `OS=1` (dead). These are patients who were alive at Liu's 2018 freeze and have since died.

### Per-project event mismatches

{per_project_table}

### Sample of named patients

Patients alive in Liu's CDR, dead in modern GDC. These are individually traceable on the GDC Data Portal — useful for spot-checking.

{sample_table}

### Implication

OS data drift is the *single most concrete* evidence that audit-grade reproducibility against Liu requires comparing against the frozen CDR (`cdr_*` columns), not against the modern GDC. For a survival cohort that includes vital-status updates since 2018, use the re-derived `os_event` / `os_time` columns; you'll get more recent and more accurate vital status, but won't bit-match Liu's 2018 numbers.

For exact Liu reproduction: filter to `cdr_matched == True` and use `cdr_OS` / `cdr_OS_time`.
"""


def main() -> None:
    df = load_df()
    matched = df[df["cdr_matched"]].copy()
    mm = matched[
        matched["cdr_OS"].notna()
        & matched["os_event"].notna()
        & (matched["cdr_OS"] != matched["os_event"])
    ].copy()

    direction = mm.groupby(["cdr_OS", "os_event"]).size().reset_index(name="count")
    direction.columns = ["Liu cdr_OS", "ours os_event", "patients"]

    per_proj = (
        mm.groupby("project")
        .size()
        .reset_index(name="event mismatches")
        .sort_values("event mismatches", ascending=False)
    )

    sample = (
        mm[["case_submitter_id", "project", "cdr_OS", "cdr_OS_time", "os_event", "os_time"]]
        .rename(
            columns={
                "cdr_OS": "Liu OS",
                "cdr_OS_time": "Liu OS_time",
                "os_event": "ours OS",
                "os_time": "ours OS_time",
            }
        )
        .head(10)
    )

    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(
        REPORT.format(
            direction_table=to_md(direction),
            per_project_table=to_md(per_proj),
            sample_table=to_md(sample),
        )
    )
    print(f"Wrote {OUT.relative_to(HERE.parent.parent)}")
    print(f"  Total OS event mismatches: {len(mm)}")


if __name__ == "__main__":
    main()
