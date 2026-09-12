"""Build and cache the Liu 2018 validation cohort.

Reads `<RAW>/TCGA-*/cases.json` for every project, attaches Liu's curated
CDR values (from the workbook — cdr_* are *not* in the published HF
datasets; we attach them here just for validation) plus our re-derived
survival values from the `survival_derived` struct, and writes three
caches:

    _cache/df.parquet      — slim per-patient DataFrame (scalar cols only)
    _cache/demo.parquet    — Table 1 demographics (incl. resolved stage/grade)
    _cache/all_rows.jsonl  — full nested patient rows (one JSON per line)

The slim DataFrame keeps a flat schema (cdr_OS, os_event, ...) so the
section scripts don't need to know about the nested struct layout used
by the published patients dataset. The struct gets unpacked here.

Section scripts read from `_cache/` via `cohort.load_*()` instead of
re-walking raw — the slow part is the 33-project JSON load + survival
attachment, which only needs to run when raw data or pipeline code
changes.

Usage:
    uv run python dev_research/liu_2018/load_cohort.py [--raw PATH]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
from tcga2hf_pipeline import cdr, clinical, clinical_supplement, survival

HERE = Path(__file__).parent
CACHE = HERE / "_cache"

# Liu's Table 1 stage-field fallback chain. AJCC pathologic for most cancers;
# clinical fallbacks for CESC/DLBC/OV/UCEC/UCS; ENSAT for ACC; Masaoka for THYM.
_STAGE_FIELDS = (
    "ajcc_pathologic_stage",
    "ajcc_clinical_stage",
    "ensat_pathologic_stage",
    "ann_arbor_pathologic_stage",
    "ann_arbor_clinical_stage",
    "masaoka_stage",
    "figo_stage",
    "igcccg_stage",
)

# Top-level scalar fields we surface to the slim per-patient DataFrame.
# The 8 derived survival fields live in the `survival_derived` struct on
# each row; we unpack them into flat columns at slim-DataFrame build time
# so section scripts can pandas-filter on `os_event` etc. directly.
_TOP_LEVEL_COLS = [
    "case_submitter_id",
    "project_id",
    "primary_site",
    "disease_type",
    "cdr_matched",
    "cdr_redaction",
    "cdr_OS",
    "cdr_OS_time",
    "cdr_DSS",
    "cdr_DSS_time",
    "cdr_PFI",
    "cdr_PFI_time",
    "cdr_DFI",
    "cdr_DFI_time",
    "cdr_survival_complete",
]
_DERIVED_SURVIVAL_COLS = [
    "os_event",
    "os_time",
    "dss_event",
    "dss_time",
    "pfi_event",
    "pfi_time",
    "dfi_event",
    "dfi_time",
]


def _resolve_stage(row: dict) -> str | None:
    primary = survival._primary_diagnosis(row) or {}
    for f in _STAGE_FIELDS:
        v = primary.get(f)
        if v:
            return v
    return None


def build_cohort(raw: Path) -> tuple[pd.DataFrame, pd.DataFrame, list[dict]]:
    cdr_index = cdr.load_cdr_index(raw)
    print(f"CDR rows indexed: {len(cdr_index)}")

    all_rows: list[dict] = []
    for proj_dir in sorted(raw.glob("TCGA-*")):
        cases_path = proj_dir / "cases.json"
        if not cases_path.exists():
            continue
        cases = json.loads(cases_path.read_text())
        rows = clinical.to_patient_rows(cases)
        cdr.attach_cdr(rows, cdr_index)
        supp_dir = proj_dir / "clinical_supplement"
        supps = clinical_supplement.load_supplements_for_project(supp_dir)
        if supps:
            clinical_supplement.attach_supplements(rows, supps)
        survival.attach_survival(rows)
        all_rows.extend(rows)
        print(f"  {proj_dir.name}: {len(rows)} patients")

    def _slim_row(r: dict) -> dict:
        out = {c: r.get(c) for c in _TOP_LEVEL_COLS}
        sd = r.get("survival_derived") or {}
        for c in _DERIVED_SURVIVAL_COLS:
            out[c] = sd.get(c)
        return out

    df = pd.DataFrame([_slim_row(r) for r in all_rows])
    df["project"] = df["project_id"].str.replace("TCGA-", "", regex=False)

    demo = pd.DataFrame(
        [
            {
                "case_submitter_id": r["case_submitter_id"],
                "project": r["project_id"].replace("TCGA-", ""),
                "vital_status": (r.get("demographic") or {}).get("vital_status"),
                "days_to_birth": (r.get("demographic") or {}).get("days_to_birth"),
                "sex_at_birth": (r.get("demographic") or {}).get("sex_at_birth"),
                "race": (r.get("demographic") or {}).get("race"),
                "ajcc_stage": (survival._primary_diagnosis(r) or {}).get("ajcc_pathologic_stage"),
                "tumor_grade": (survival._primary_diagnosis(r) or {}).get("tumor_grade"),
                "stage_raw": _resolve_stage(r),
            }
            for r in all_rows
        ]
    )

    return df, demo, all_rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--raw",
        type=Path,
        default=Path.home() / "data" / "tcga2hf" / "raw",
        help="Path to <data-dir>/raw/ holding TCGA-*/cases.json + cdr/.",
    )
    args = parser.parse_args()

    CACHE.mkdir(exist_ok=True)
    df, demo, all_rows = build_cohort(args.raw)

    df_path = CACHE / "df.parquet"
    demo_path = CACHE / "demo.parquet"
    rows_path = CACHE / "all_rows.jsonl"

    df.to_parquet(df_path)
    demo.to_parquet(demo_path)
    with rows_path.open("w") as f:
        for r in all_rows:
            f.write(json.dumps(r, default=str))
            f.write("\n")

    print()
    print(f"Wrote {len(df):>6} rows  -> {df_path.relative_to(HERE)}")
    print(f"Wrote {len(demo):>6} rows  -> {demo_path.relative_to(HERE)}")
    print(f"Wrote {len(all_rows):>6} rows  -> {rows_path.relative_to(HERE)}")
    print(f"\nCDR-matched: {df['cdr_matched'].sum()}  /  post-freeze: {(~df['cdr_matched']).sum()}")


if __name__ == "__main__":
    main()
