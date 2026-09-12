"""Section 4 — Per-project survival landscape (numerical analog of Liu's Figure 1B-E).

Liu's Figure 1B-E showed Kaplan-Meier curves of OS / PFI / DFI / DSS for
each of the 33 tumor types overlaid on one axes per endpoint. We don't
ship figures in this report; instead we surface the same information as
event-counts plus survival rates at clinically meaningful timepoints
(1-year, 3-year, 5-year). Aggressive cancers (SKCM, OV, GBM, MESO) drop
fast; indolent cancers (TGCT, PRAD, THCA, KICH) stay flat — same pattern
the K-M panels showed visually.

Computed against our **re-derived** endpoints on the full current cohort
(11,428 patients). Tumor types Liu marks as having no usable signal for
an endpoint (e.g. SKCM/THYM/UVM/LAML for DFI) are omitted from that
column.

Writes: sections/04_survival_landscape.md
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from cohort import ENDPOINTS, load_df, to_md

HERE = Path(__file__).parent
OUT = HERE / "sections" / "04_survival_landscape.md"


def _km_rate_at(times: pd.Series, events: pd.Series, day: float) -> float | None:
    """Kaplan-Meier survival probability at `day`, with right-censoring.

    Standard product-limit estimator: at each event time t_i, survival is
    multiplied by (1 - d_i / n_i) where n_i is the at-risk count and d_i
    the events at t_i. We accept (time, event) pairs and return the K-M
    curve evaluated at the given day in days.
    """
    df = pd.DataFrame({"t": times, "e": events}).dropna()
    if df.empty:
        return None
    df = df.sort_values("t").reset_index(drop=True)
    surv = 1.0
    for t, sub in df[df["e"] == 1].groupby("t"):
        d = len(sub)
        at_risk = (df["t"] >= t).sum()
        if at_risk == 0:
            break
        surv *= 1 - d / at_risk
        if t >= day:
            return round(surv, 3)
    # If the curve never crosses `day`, return survival at the last event;
    # there's no late-time information beyond max(times).
    if df["t"].max() < day:
        # Censored before reaching `day` — return None to flag insufficient FU.
        return None
    return round(surv, 3)


def _format_pct(v: float | None) -> str:
    return "—" if v is None else f"{int(round(v * 100))}%"


def build_landscape(df: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """One DataFrame per endpoint: per-project (N, events, 1y/3y/5y survival rate)."""
    out: dict[str, pd.DataFrame] = {}
    for ep in ENDPOINTS:
        ev_col = f"{ep.lower()}_event"
        t_col = f"{ep.lower()}_time"
        rows = []
        for proj, sub in df.groupby("project"):
            populated = sub.dropna(subset=[ev_col, t_col])
            n = len(populated)
            if n == 0:
                continue
            events = int(populated[ev_col].sum())
            r1 = _km_rate_at(populated[t_col], populated[ev_col], 365.25)
            r3 = _km_rate_at(populated[t_col], populated[ev_col], 365.25 * 3)
            r5 = _km_rate_at(populated[t_col], populated[ev_col], 365.25 * 5)
            rows.append(
                {
                    "project": proj,
                    "N (populated)": n,
                    "events": events,
                    "1-yr": _format_pct(r1),
                    "3-yr": _format_pct(r3),
                    "5-yr": _format_pct(r5),
                }
            )
        out[ep] = pd.DataFrame(rows).sort_values("project").reset_index(drop=True)
    return out


REPORT = """\
## Section 4 — Per-project survival landscape

Liu's Figure 1B-E shows Kaplan-Meier curves per project for each endpoint. We surface the same information as event counts plus survival probabilities at 1-year, 3-year, and 5-year timepoints — the standard clinical reporting milestones. The same patterns Liu's panels show visually appear here numerically: aggressive cancers (SKCM, OV, GBM, MESO) drop fast; indolent cancers (TGCT, PRAD, THCA, KICH) stay flat.

Computed against our re-derived `{{os,dss,pfi,dfi}}_event/_time` columns on the full current cohort (11,428 patients). `—` means insufficient follow-up to estimate the rate at that timepoint (curve hasn't reached it). Projects with zero populated patients for an endpoint (e.g. LAML / SKCM / THYM / UVM for DFI — Liu's exclusions) are omitted from that table.

### OS — Overall Survival

{os_table}

### PFI — Progression-Free Interval

{pfi_table}

### DFI — Disease-Free Interval

{dfi_table}

### DSS — Disease-Specific Survival

{dss_table}

### Reading the tables

- **N (populated)** — patients with both event and time populated for that endpoint. For DFI especially this is much smaller than the project N (most patients lack a disease-free signal at end of first course; see Section 9 deep-dive).
- **events** — count of `event=1` patients in the populated subset.
- **1-yr / 3-yr / 5-yr** — Kaplan-Meier survival probability at that timepoint, accounting for right-censoring. `—` flags projects whose follow-up is shorter than the timepoint (e.g. mostly post-2018 cases).
- **OS vs DSS** — DSS censors cancer-unrelated deaths, so DSS rates are equal to or higher than OS rates at the same timepoint within a project.
- **OS vs PFI** — PFI events are generally earlier than OS events (Liu: *"PFI is generally considered a more informative endpoint"*); PFI rates at 1-year are typically lower than OS rates at 1-year.

### Cross-reference to Liu's reliability assessment

Liu's Table 3 marks each (project, endpoint) as recommended / acceptable / not recommended based on event count and assumption tests. The N-populated column above is the relevant input — projects with very low N for an endpoint don't pass Liu's reliability bar. See the original paper Table 3 for the per-cell recommendations; we don't reproduce that here because it's an editorial overlay rather than a reproducible computation.
"""


def main() -> None:
    df = load_df()
    landscape = build_landscape(df)
    OUT.parent.mkdir(exist_ok=True)
    body = REPORT.format(
        os_table=to_md(landscape["OS"]),
        dss_table=to_md(landscape["DSS"]),
        pfi_table=to_md(landscape["PFI"]),
        dfi_table=to_md(landscape["DFI"]),
    )
    OUT.write_text(body)
    print(f"Wrote {OUT.relative_to(HERE.parent.parent)}")
    for ep in ENDPOINTS:
        n_proj = len(landscape[ep])
        print(f"  {ep}: {n_proj} projects with populated cases")


if __name__ == "__main__":
    main()
