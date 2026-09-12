"""Reader interface for the Liu 2018 validation cache.

Section scripts call `load_df()` / `load_demo()` / `load_all_rows()` to
get the cached cohort without re-walking raw data. The cache is built by
`load_cohort.py`; if the parquet files don't exist, instruct the user to
run that first.
"""

from __future__ import annotations

import json
from collections.abc import Iterator
from pathlib import Path

import pandas as pd

HERE = Path(__file__).parent
CACHE = HERE / "_cache"

DAYS_PER_MONTH = 30.44
ENDPOINTS = ["OS", "DSS", "PFI", "DFI"]


def _check_cache() -> None:
    if not (CACHE / "df.parquet").exists():
        raise FileNotFoundError(
            f"Cache not found at {CACHE}. Build it first with:\n"
            f"    uv run python dev_research/liu_2018/load_cohort.py"
        )


def load_df() -> pd.DataFrame:
    """Slim per-patient DataFrame (scalar cols + cdr_*/<ep>_event/<ep>_time)."""
    _check_cache()
    return pd.read_parquet(CACHE / "df.parquet")


def load_demo() -> pd.DataFrame:
    """Table 1 demographics (vital status, age proxy, race, stage_raw, grade)."""
    _check_cache()
    return pd.read_parquet(CACHE / "demo.parquet")


def iter_all_rows() -> Iterator[dict]:
    """Stream the full nested patient rows (one per JSON line)."""
    _check_cache()
    path = CACHE / "all_rows.jsonl"
    with path.open() as f:
        for line in f:
            yield json.loads(line)


def load_all_rows() -> list[dict]:
    """Load every nested patient row into memory (~420MB JSONL → list)."""
    return list(iter_all_rows())


def to_md(df: pd.DataFrame, index: bool = False) -> str:
    """Render a DataFrame as a GitHub-flavored markdown table.

    Avoids the `tabulate` dependency that pandas' `to_markdown` requires.
    """
    cols = ([df.index.name or ""] + list(df.columns)) if index else list(df.columns)
    header = "| " + " | ".join(str(c) for c in cols) + " |"
    sep = "| " + " | ".join("---" for _ in cols) + " |"
    lines = [header, sep]
    for idx, row in df.iterrows():
        vals = ([str(idx)] if index else []) + ["" if pd.isna(v) else str(v) for v in row.tolist()]
        lines.append("| " + " | ".join(vals) + " |")
    return "\n".join(lines)
