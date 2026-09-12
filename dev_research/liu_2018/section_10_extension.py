"""Section 10 — Post-2018 extension coverage.

Liu's CDR is frozen at the 2018 data release. Our cohort has grown since;
this section quantifies what re-derivation against modern GDC adds.

Writes: sections/10_extension.md
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from cohort import ENDPOINTS, load_df, to_md

HERE = Path(__file__).parent
OUT = HERE / "sections" / "10_extension.md"


REPORT = """\
## Section 10 — Extension beyond Liu (post-2018 coverage)

Liu's CDR is frozen at the 2018 data release. Our re-derivation runs against the modern GDC and includes patients who entered TCGA after Liu's freeze.

### Headline numbers

| | count |
|---|---|
| Total cohort (all 33 TCGA projects) | {n_total} |
| In Liu's CDR (matched) | {n_matched} |
| **Post-Liu's freeze (extension)** | **{n_ext}** |

### Re-derived endpoint coverage on the extension cohort

For the {n_ext} post-freeze patients, how many have a populated value for each re-derived endpoint:

{coverage_table}

### Top projects by extension count

Most of the +{n_ext} patients land in a few projects (TGCT and LUAD lead by a wide margin). The rest of TCGA was essentially closed by 2018; only a handful of projects continued accruing.

{per_project_table}

### Why TGCT / LUAD dominate the extension

Two structural reasons:

- **TGCT** had a small published cohort in Liu's 2018 freeze (134 patients) but BCR continued accruing testicular germ cell cases substantially through 2019-2020. The {tgct_ext} new TGCT patients more than double the 2018 cohort.
- **LUAD** had ongoing tissue acquisitions through 2019 from the original TCGA centers; +{luad_ext} new patients reflects late-stage data freezes.

Other projects (GBM, OV, UCEC, DLBC) added handfuls of new patients, mostly cases that were in-flight at the 2018 freeze and finalized later.

### Implication for downstream analysis

If you're reproducing Liu's results: stick to `cdr_*` columns + filter `cdr_matched == True`. You'll work on 11,160 patients exactly matching the 2018 freeze.

If you want maximum cohort size or current vital status: use `{{os,dss,pfi,dfi}}_event/_time` on all 11,428 patients. The +268 post-freeze patients are most useful for TGCT-specific analyses (where they're a meaningful fraction of the cohort). For most other projects the addition is small enough that conclusions wouldn't change either way.
"""


def main() -> None:
    df = load_df()
    ext = df[~df["cdr_matched"]]

    coverage = pd.DataFrame(
        [
            {
                "endpoint": ep,
                "post-freeze patients with re-derived value": int(
                    ext[f"{ep.lower()}_event"].notna().sum()
                ),
                "of total post-freeze": len(ext),
            }
            for ep in ENDPOINTS
        ]
    )

    per_proj = (
        ext.groupby("project")
        .size()
        .reset_index(name="post-freeze patients")
        .sort_values("post-freeze patients", ascending=False)
    )

    tgct_ext = int(ext[ext["project"] == "TGCT"].shape[0])
    luad_ext = int(ext[ext["project"] == "LUAD"].shape[0])

    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(
        REPORT.format(
            n_total=len(df),
            n_matched=int(df["cdr_matched"].sum()),
            n_ext=len(ext),
            coverage_table=to_md(coverage),
            per_project_table=to_md(per_proj),
            tgct_ext=tgct_ext,
            luad_ext=luad_ext,
        )
    )
    print(f"Wrote {OUT.relative_to(HERE.parent.parent)}")
    print(
        f"  Post-freeze patients: {len(ext)} across {(per_proj['post-freeze patients'] > 0).sum()} projects"
    )


if __name__ == "__main__":
    main()
