"""List patients where any sample's `days_to_collection` exceeds the patient's
`days_to_death`. Prints one row per such patient with the relevant numbers
and a GDC Data Portal URL for spot-checking.

Usage:
    uv run python scripts/list_bcr_after_death.py
    uv run python scripts/list_bcr_after_death.py --data-dir /some/other/path
    uv run python scripts/list_bcr_after_death.py --project TCGA-CHOL

Reads `<data-dir>/processed/<PROJECT>/data.parquet` for each TCGA-* project
present, so run `ds2hf-pipeline tcga build` first.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pyarrow.parquet as pq

GDC_PORTAL_CASE = "https://portal.gdc.cancer.gov/cases"


def list_for_project(project_dir: Path) -> list[dict]:
    rows = pq.read_table(project_dir / "data.parquet").to_pylist()
    out = []
    for r in rows:
        demo = r.get("demographic") or {}
        dod = demo.get("days_to_death")
        if dod is None:
            continue
        dtcs = sorted(
            {
                s["days_to_collection"]
                for s in r["samples"]
                if s.get("days_to_collection") is not None
            }
        )
        if not dtcs or max(dtcs) <= dod:
            continue
        out.append(
            {
                "project_id": r["project_id"],
                "case_id": r["case_id"],
                "case_submitter_id": r["case_submitter_id"],
                "days_to_death": dod,
                "days_to_collection": dtcs,
                "delta": max(dtcs) - dod,
            }
        )
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--data-dir",
        type=Path,
        default=Path.home() / "data/tcga2hf",
        help="Root data dir (default: $HOME/data/tcga2hf)",
    )
    ap.add_argument(
        "--project",
        action="append",
        help="Filter to specific project ids (repeatable). Default: all under processed/",
    )
    args = ap.parse_args()

    processed = args.data_dir / "processed"
    project_dirs = sorted(processed.glob("TCGA-*"))
    if args.project:
        project_dirs = [p for p in project_dirs if p.name in args.project]
    if not project_dirs:
        raise SystemExit(f"no project parquets found under {processed}")

    rows: list[dict] = []
    for pd in project_dirs:
        rows.extend(list_for_project(pd))

    if not rows:
        print("(no patients found where days_to_collection > days_to_death)")
        return

    rows.sort(key=lambda r: -r["delta"])  # largest gap first

    header = (
        f"{'project':<10} {'case':<14} {'d_death':>8} {'d_collection':>14} {'delta':>7}  gdc_portal"
    )
    print(header)
    print("-" * len(header))
    for r in rows:
        url = f"{GDC_PORTAL_CASE}/{r['case_id']}"
        dtcs = r["days_to_collection"]
        coll = dtcs[0] if len(dtcs) == 1 else dtcs
        print(
            f"{r['project_id']:<10} {r['case_submitter_id']:<14} "
            f"{r['days_to_death']:>8} {str(coll):>14} {r['delta']:>+7}  {url}"
        )

    print(f"\n{len(rows)} patients total")


if __name__ == "__main__":
    main()
