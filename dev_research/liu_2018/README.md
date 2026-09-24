# Liu et al. 2018 — TCGA-CDR reproduction

Notes, transcribed reference data, and reproduction scripts for [Liu J, Lichtenberg T, Hoadley KA, et al. *An Integrated TCGA Pan-Cancer Clinical Data Resource to Drive High-Quality Survival Outcome Analytics.* Cell 173(2):400–416 (2018)](https://doi.org/10.1016/j.cell.2018.02.052).

## Files

- [`notes.md`](notes.md) — paper summary + how this paper enters the pipeline.
- [`PIIS0092867418302290.pdf`](PIIS0092867418302290.pdf) — local copy of the paper.
- [`liu_table1.csv`](liu_table1.csv) — Liu's published Table 1 transcribed verbatim (used as ground-truth for Section 1's bucketing-sanity check).
- [`liu_table2.csv`](liu_table2.csv) — Liu's published Table 2 transcribed verbatim (used as ground-truth for Section 2's bucketing-sanity check).
- [`load_cohort.py`](load_cohort.py) — walks raw GDC data + CDR workbook + BCR biotab Clinical Supplements once, writes `_cache/{df.parquet, demo.parquet, all_rows.jsonl}`. Run after raw data or pipeline code changes; sections re-run from the cache in <1s each.
- [`cohort.py`](cohort.py) — reader interface for the cache (`load_df`, `load_demo`, `iter_all_rows`, `to_md` markdown helper).
- [`section_NN_*.py`](.) — one script per report section. Each reads the cache and writes its slice to `sections/NN_*.md` (markdown is the editable, git-friendly source).
- [`sections/`](sections/) — per-section markdown outputs.
- [`assemble_html.py`](assemble_html.py) — converts `sections/*.md` into a single self-contained [`report.html`](report.html) with sticky TOC, sortable tables, embedded CSS/JS (no external assets).
- [`report.html`](report.html) — assembled report. Open in any browser; works offline.

## Run

```sh
cd dev_research/liu_2018

# One-time (rebuild after raw data or pipeline code changes):
uv run python load_cohort.py

# Iterate per section (each takes <1s once the cache exists):
uv run python section_01_cohort.py
uv run python section_02_followup.py
# ... section_03 through 11 ...

# Assemble report.html (the `report` extra adds python-markdown):
uv run --extra report python assemble_html.py
```

The cache (`_cache/*.parquet`, `_cache/*.jsonl`) is gitignored. Section markdown outputs (`sections/*.md`) and the assembled `report.html` are committed.

## Section structure

Each section script follows the same pattern:

1. Imports `cohort.load_*` to read the cache.
2. Runs an analysis (table reproduction, drift comparison, deep-dive).
3. Writes `sections/NN_<name>.md` with narrative + tables.
4. Prints a one-line summary stat to stdout.

Adding a new section: copy an existing `section_NN_*.py`, change the analysis logic. `assemble_html.py` picks up `sections/[0-9][0-9]_*.md` in lexicographic order automatically — no registration needed.

## Why scripts + HTML instead of a notebook

The original `cdr_validation.ipynb` (now removed) reloaded all 33 projects' raw data on every kernel restart (~30 seconds). The script-based version caches the cohort once via `load_cohort.py` and lets each section re-run in <1s, making iteration much faster. The HTML output gives nicer rendering than markdown (sticky TOC, sortable tables, syntax-highlighted code) while staying portable as a single file.

The trade-off is no in-line plots — section 4 surfaces survival rates as numerical tables (1-yr / 3-yr / 5-yr survival per project per endpoint) instead of K-M curves, which is enough signal for the validation goal. If we ever want figures back, `assemble_html.py` would be a natural place to add embedded SVG output from matplotlib (write SVG to a string, inline into the HTML).
