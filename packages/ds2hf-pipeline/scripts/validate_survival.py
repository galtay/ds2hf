"""Compare re-derived survival endpoints to Liu et al. 2018 CDR per project.

Reads `<data-dir>/raw/*/cases.json` plus the cached CDR workbook, runs
`survival.attach_survival` against each case, and reports per-project
match/mismatch/NA counts versus `cdr_*` for OS / DSS / PFI / DFI. Used
during Phase 2 development to iterate the algorithm; rerun after any
change to `ds2hf_pipeline.tcga.survival`.

Usage:
    uv run python -m ds2hf_pipeline.tcga.scripts.validate_survival
        [--data-dir PATH]
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

from ds2hf_pipeline.tcga import cdr, clinical, clinical_supplement, survival

ENDPOINTS = ["OS", "DSS", "PFI", "DFI"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path.home() / "data" / "tcga2hf",
        help="Root data dir (containing raw/<project>/cases.json and raw/cdr/).",
    )
    args = parser.parse_args()
    raw = args.data_dir / "raw"

    idx = cdr.load_cdr_index(raw)
    print(f"CDR rows indexed: {len(idx)}")

    projects = sorted(p.name for p in raw.iterdir() if p.is_dir() and p.name.startswith("TCGA-"))

    overall: dict[str, dict[str, int]] = defaultdict(
        lambda: {
            "match": 0,
            "mismatch": 0,
            "both_na": 0,
            "cdr_pop_der_na": 0,
            "der_pop_cdr_na": 0,
            "total_matched": 0,
        }
    )

    header_bits = "  ".join(f"{ep}={'match/total':>14}" for ep in ENDPOINTS)
    print(f"\n{'project':<12} {'n':>5}  {header_bits}")
    for proj in projects:
        cases_path = raw / proj / "cases.json"
        if not cases_path.exists():
            continue
        cases = json.loads(cases_path.read_text())
        rows = clinical.to_patient_rows(cases)
        # CDR values are no longer attached to rows by the build pipeline
        # (the published datasets are TCGA/GDC + derived only). For
        # validation against Liu's frozen 2018 values we attach them here
        # as transient `_cdr_*` keys, used only by this script.
        cdr.attach_cdr(rows, idx)
        supps = clinical_supplement.load_supplements_for_project(raw / proj / "clinical_supplement")
        if supps:
            clinical_supplement.attach_supplements(rows, supps)
        survival.attach_survival(rows)

        matched = [r for r in rows if r["cdr_matched"]]
        proj_stats = {
            ep: {"match": 0, "mismatch": 0, "both_na": 0, "cdr_pop_der_na": 0, "der_pop_cdr_na": 0}
            for ep in ENDPOINTS
        }
        for r in matched:
            sd = r.get("survival_derived") or {}
            for ep in ENDPOINTS:
                cdr_e, cdr_t = r[f"cdr_{ep}"], r[f"cdr_{ep}_time"]
                der_e, der_t = sd.get(f"{ep.lower()}_event"), sd.get(f"{ep.lower()}_time")
                overall[ep]["total_matched"] += 1
                if cdr_e is None and der_e is None:
                    proj_stats[ep]["both_na"] += 1
                elif cdr_e is None:
                    proj_stats[ep]["der_pop_cdr_na"] += 1
                elif der_e is None:
                    proj_stats[ep]["cdr_pop_der_na"] += 1
                elif cdr_e == der_e and (
                    cdr_t == der_t
                    or (cdr_t is not None and der_t is not None and abs(cdr_t - der_t) < 0.5)
                ):
                    proj_stats[ep]["match"] += 1
                else:
                    proj_stats[ep]["mismatch"] += 1
        for ep in ENDPOINTS:
            for k, v in proj_stats[ep].items():
                overall[ep][k] += v

        n = len(matched)
        bits = "  ".join(f"{ep}={proj_stats[ep]['match']:>5}/{n:<5}" for ep in ENDPOINTS)
        print(f"  {proj:<10} {n:>5}  {bits}")

    print("\n=== OVERALL ===")
    for ep in ENDPOINTS:
        s = overall[ep]
        total = s["total_matched"]
        pct = 100 * s["match"] / total if total else 0
        print(
            f"{ep}: {s['match']:>5}/{total} match ({pct:5.1f}%)  "
            f"mismatch={s['mismatch']:>4}  both_na={s['both_na']:>5}  "
            f"cdr_pop_der_na={s['cdr_pop_der_na']:>5}  der_pop_cdr_na={s['der_pop_cdr_na']:>5}"
        )


if __name__ == "__main__":
    main()
