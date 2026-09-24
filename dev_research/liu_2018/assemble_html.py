"""Render sections/*.md into a single self-contained report.html.

Wraps each section's converted markdown in semantic HTML, with:
  - sticky TOC sidebar driven by section H2s
  - sortable tables (click any column header to sort; click again to reverse)
  - clean responsive typography (system font stack, no external fonts)
  - inline CSS + JS so report.html can be emailed / opened offline

Usage:
    uv run python dev_research/liu_2018/assemble_html.py
"""

from __future__ import annotations

import re
from datetime import UTC, datetime
from pathlib import Path

import markdown

HERE = Path(__file__).parent
SECTIONS = HERE / "sections"
OUT = HERE / "report.html"

INTRO_MD = """\
Liu et al. curated four survival endpoints — Overall Survival (OS), Disease-Specific Survival (DSS), Disease-Free Interval (DFI), and Progression-Free Interval (PFI) — for 11,160 TCGA patients across 33 cancer types. Their result is the canonical *TCGA-CDR* table, frozen at a 2018 data freeze.

This report compares two sources of the same four endpoints:

1. **Liu's values (curated, frozen)** — verbatim from `TCGA-CDR-SupplementalTableS1.xlsx`. Published as-is in [`gabrielaltay/tcga-pancanatlas-cdr`](https://huggingface.co/datasets/gabrielaltay/tcga-pancanatlas-cdr). Throughout this report they appear as `cdr_OS`, `cdr_OS_time`, etc.; that prefix exists only in the report's working cache, not in any published dataset. ~268 of our 11,428 patients post-date Liu's freeze and have no CDR row.
2. **Our re-derived values (live)** — `{os,dss,pfi,dfi}_event` / `_time`, recomputed from current GDC data using Liu's documented algorithm (`ds2hf_pipeline.tcga.survival`), augmented with `treatment_outcome_first_course` from the BCR biotab Clinical Supplements. Published in the `survival_derived` struct of `gabrielaltay/tcga-patients-open` and the `survival_derived` table of `gabrielaltay/tcga-tabular-open`.

This report walks through three things:

- **Reproductions** — Tables 1, 2 and a per-project survival landscape (numerical analog of Liu's Figure 1).
- **Validations** — per-endpoint per-project agreement between curated and re-derived values, plus deep-dives into each endpoint's failure modes.
- **Extensions** — what re-deriving from current GDC data buys us that Liu's frozen 2018 table doesn't.

## TL;DR — agreement against Liu's curated values

On the 11,160 CDR-matched patients, defining **agreement** as either (a) both correctly NA, or (b) both populated and event direction agrees within 30 days:

| Endpoint | Agreement rate | Notes |
|---|---|---|
| **OS**  | **98.2%** | Remaining 2% is data drift — patients alive at Liu's 2018 freeze who have since died. |
| **DSS** | **93.2%** | Liu flagged DSS as approximate; we follow the same definition. |
| **PFI** | **96.3%** | Past 95% — meets Liu's reliability bar. |
| **DFI** | **90.1%** | Up from 77.2% pre-Clinical-Supplement integration. Most of the remaining 10% is patients we have *extra coverage* on, not active contradictions. |

The DFI gain came from adding the BCR biotab Clinical Supplement fetcher: under-population (patients Liu populated but we didn't) collapsed from 1,625 to 52, a 97% reduction. The harmonized GDC `/cases?expand=...` API drops Liu's `treatment_outcome_first_course` field; the BCR biotab files preserve it. See Section 3 for the field-rename map and Section 9 for the DFI deep-dive.
"""

CSS = """\
:root {
    --bg: #fafbfc;
    --fg: #1a1f24;
    --muted: #6a737d;
    --accent: #0366d6;
    --border: #e1e4e8;
    --code-bg: #f6f8fa;
    --table-stripe: #fafbfc;
    --table-hover: #eef4fb;
    --callout-bg: #f0f6ff;
    --callout-border: #cce4ff;
    --max-content: 880px;
    --sidebar-width: 280px;
}

@media (prefers-color-scheme: dark) {
    :root {
        --bg: #0d1117;
        --fg: #e6edf3;
        --muted: #8b949e;
        --accent: #58a6ff;
        --border: #30363d;
        --code-bg: #161b22;
        --table-stripe: #0f141a;
        --table-hover: #1c2733;
        --callout-bg: #0f1a2c;
        --callout-border: #1f3a5f;
    }
}

* { box-sizing: border-box; }

html { scroll-behavior: smooth; scroll-padding-top: 1rem; }

body {
    margin: 0;
    background: var(--bg);
    color: var(--fg);
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", system-ui,
                 "Helvetica Neue", Arial, sans-serif;
    line-height: 1.6;
    font-size: 16px;
}

.layout {
    display: grid;
    grid-template-columns: var(--sidebar-width) minmax(0, 1fr);
    max-width: calc(var(--sidebar-width) + var(--max-content) + 4rem);
    margin: 0 auto;
    gap: 2rem;
    padding: 2rem 1.5rem;
}

@media (max-width: 800px) {
    .layout { grid-template-columns: 1fr; padding: 1rem; }
    nav.toc { position: static !important; max-height: none !important; }
}

nav.toc {
    position: sticky;
    top: 1.5rem;
    align-self: start;
    max-height: calc(100vh - 3rem);
    overflow-y: auto;
    font-size: 0.92rem;
    padding-right: 0.5rem;
    border-right: 1px solid var(--border);
}

nav.toc h2 {
    font-size: 0.78rem;
    font-weight: 600;
    text-transform: uppercase;
    letter-spacing: 0.05em;
    color: var(--muted);
    margin: 0 0 0.6rem 0;
}

nav.toc ol {
    list-style: none;
    padding: 0;
    margin: 0;
}

nav.toc li { margin: 0.25rem 0; }

nav.toc a {
    color: var(--fg);
    text-decoration: none;
    display: block;
    padding: 0.25rem 0.5rem;
    border-radius: 4px;
    border-left: 2px solid transparent;
}

nav.toc a:hover { background: var(--table-hover); border-left-color: var(--accent); }

main {
    min-width: 0;  /* enables tables to shrink-to-fit instead of overflowing */
}

main > header {
    margin-bottom: 2rem;
    padding-bottom: 1.5rem;
    border-bottom: 1px solid var(--border);
}

main > header h1 {
    margin: 0 0 0.5rem 0;
    font-size: 2rem;
    font-weight: 700;
    line-height: 1.2;
}

main > header .subtitle {
    color: var(--muted);
    font-size: 0.95rem;
    margin: 0;
}

main > header .meta {
    color: var(--muted);
    font-size: 0.82rem;
    margin-top: 0.75rem;
}

h2, h3, h4 {
    font-weight: 600;
    line-height: 1.3;
    margin-top: 2rem;
    margin-bottom: 0.75rem;
}

h2 {
    font-size: 1.5rem;
    padding-bottom: 0.4rem;
    border-bottom: 1px solid var(--border);
}

h3 { font-size: 1.2rem; }
h4 { font-size: 1.0rem; color: var(--muted); }

p, li { line-height: 1.6; }

a { color: var(--accent); text-decoration: none; }
a:hover { text-decoration: underline; }

code {
    font-family: ui-monospace, SFMono-Regular, "SF Mono", Menlo, Monaco,
                 Consolas, "Liberation Mono", "Courier New", monospace;
    font-size: 0.88em;
    background: var(--code-bg);
    padding: 0.15em 0.35em;
    border-radius: 3px;
}

pre {
    background: var(--code-bg);
    padding: 1rem;
    border-radius: 6px;
    overflow-x: auto;
    border: 1px solid var(--border);
}

pre code {
    background: none;
    padding: 0;
    font-size: 0.85rem;
}

/* Tables: sortable on header click; alternating rows; horizontal-scroll on overflow */
.table-wrap {
    overflow-x: auto;
    margin: 1rem 0;
    border: 1px solid var(--border);
    border-radius: 6px;
}

table {
    border-collapse: collapse;
    width: 100%;
    font-size: 0.92rem;
}

table th, table td {
    text-align: left;
    padding: 0.5rem 0.75rem;
    border-bottom: 1px solid var(--border);
    vertical-align: top;
}

table th {
    background: var(--code-bg);
    font-weight: 600;
    cursor: pointer;
    user-select: none;
    position: relative;
    white-space: nowrap;
}

table th:hover { background: var(--callout-bg); }

table th::after {
    content: "↕";
    color: var(--muted);
    font-size: 0.7em;
    margin-left: 0.4em;
    opacity: 0.5;
}

table th.sort-asc::after { content: "▲"; opacity: 1; }
table th.sort-desc::after { content: "▼"; opacity: 1; }

table tbody tr:nth-child(even) { background: var(--table-stripe); }
table tbody tr:hover { background: var(--table-hover); }

/* Make the first column visually distinct in wide tables (project / endpoint) */
table tbody td:first-child { font-weight: 500; }

blockquote {
    margin: 1rem 0;
    padding: 0.75rem 1rem;
    background: var(--callout-bg);
    border-left: 3px solid var(--accent);
    border-radius: 0 4px 4px 0;
    color: var(--fg);
}

blockquote p { margin: 0; }

hr {
    border: none;
    border-top: 1px solid var(--border);
    margin: 2.5rem 0;
}

ul, ol { padding-left: 1.5rem; }

footer {
    margin-top: 3rem;
    padding-top: 1.5rem;
    border-top: 1px solid var(--border);
    color: var(--muted);
    font-size: 0.85rem;
    text-align: center;
}
"""

JS_SORTABLE = """\
// Click-to-sort for every table on the page. Numeric columns auto-detected
// (parseFloat first cell); fallback to locale string compare. Click again
// on the same column to reverse direction. Stable sort using Array.sort,
// which is stable in Node 12+ / modern browsers (V8 / SpiderMonkey).
document.querySelectorAll("table").forEach(table => {
    const tbody = table.querySelector("tbody");
    if (!tbody) return;
    table.querySelectorAll("th").forEach((th, colIdx) => {
        th.addEventListener("click", () => {
            const wasAsc = th.classList.contains("sort-asc");
            // Clear all sort indicators on this table
            table.querySelectorAll("th").forEach(h => {
                h.classList.remove("sort-asc", "sort-desc");
            });
            const dir = wasAsc ? "desc" : "asc";
            th.classList.add(dir === "asc" ? "sort-asc" : "sort-desc");

            const rows = Array.from(tbody.querySelectorAll("tr"));
            // Detect numeric column: any sample cell parses as number?
            const samples = rows.slice(0, Math.min(rows.length, 5))
                .map(r => r.children[colIdx]?.textContent.trim() ?? "");
            const isNumeric = samples.every(s => s === "" || s === "—" || s === "NA"
                || !isNaN(parseFloat(s)));
            const cmp = (a, b) => {
                const av = a.children[colIdx]?.textContent.trim() ?? "";
                const bv = b.children[colIdx]?.textContent.trim() ?? "";
                if (isNumeric) {
                    const an = parseFloat(av);
                    const bn = parseFloat(bv);
                    if (isNaN(an) && isNaN(bn)) return 0;
                    if (isNaN(an)) return 1;
                    if (isNaN(bn)) return -1;
                    return an - bn;
                }
                return av.localeCompare(bv);
            };
            rows.sort((a, b) => dir === "asc" ? cmp(a, b) : cmp(b, a));
            rows.forEach(r => tbody.appendChild(r));
        });
    });
});
"""


def _render_markdown(md_text: str) -> str:
    """Convert markdown to HTML using python-markdown with table + fenced-code support."""
    return markdown.markdown(
        md_text,
        extensions=["tables", "fenced_code", "sane_lists", "smarty"],
    )


def _wrap_tables(html: str) -> str:
    """Wrap each <table> in a <div class='table-wrap'> for horizontal scroll on overflow."""
    return re.sub(
        r"(<table[^>]*>.*?</table>)",
        r'<div class="table-wrap">\1</div>',
        html,
        flags=re.DOTALL,
    )


def _slugify(s: str) -> str:
    s = s.lower()
    s = re.sub(r"[^\w\s-]", "", s)
    return re.sub(r"[\s_-]+", "-", s).strip("-")


def _add_section_anchors(html: str) -> tuple[str, list[tuple[str, str]]]:
    """Add id="..." to every <h2> and return [(slug, title), ...] for the TOC."""
    toc: list[tuple[str, str]] = []

    def repl(m: re.Match[str]) -> str:
        title = re.sub(r"<[^>]+>", "", m.group(1))
        slug = _slugify(title)
        toc.append((slug, title))
        return f'<h2 id="{slug}">{m.group(1)}</h2>'

    return re.sub(r"<h2>(.*?)</h2>", repl, html), toc


def _build_toc(toc_entries: list[tuple[str, str]]) -> str:
    items = "\n".join(
        f'        <li><a href="#{slug}">{title}</a></li>' for slug, title in toc_entries
    )
    return f"""\
<nav class="toc" aria-label="Table of contents">
    <h2>Contents</h2>
    <ol>
{items}
    </ol>
</nav>
"""


def main() -> None:
    files = sorted(SECTIONS.glob("[0-9][0-9]_*.md"))
    if not files:
        raise SystemExit(f"No section files found in {SECTIONS}")

    parts: list[str] = [_render_markdown(INTRO_MD)]
    for f in files:
        parts.append("<hr>")
        parts.append(_render_markdown(f.read_text()))
    body_html = "\n".join(parts)
    body_html = _wrap_tables(body_html)
    body_html, toc_entries = _add_section_anchors(body_html)
    toc_html = _build_toc(toc_entries)

    timestamp = datetime.now(UTC).strftime("%Y-%m-%d %H:%M UTC")
    title = "Reproducing Liu et al. 2018 — TCGA Pan-Cancer Clinical Data Resource"
    subtitle = (
        "Liu J, Lichtenberg T, Hoadley KA, et al. "
        "<em>An Integrated TCGA Pan-Cancer Clinical Data Resource to Drive "
        "High-Quality Survival Outcome Analytics.</em> "
        "<strong>Cell</strong> 173(2):400–416 (2018). "
        '<a href="https://doi.org/10.1016/j.cell.2018.02.052">'
        "DOI 10.1016/j.cell.2018.02.052</a>"
    )

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>{title}</title>
    <style>
{CSS}    </style>
</head>
<body>
    <div class="layout">
{toc_html}
        <main>
            <header>
                <h1>{title}</h1>
                <p class="subtitle">{subtitle}</p>
                <p class="meta">Rendered {timestamp} from
                    <code>dev_research/liu_2018/sections/*.md</code> via
                    <code>assemble_html.py</code>.</p>
            </header>
{body_html}
            <footer>Built from the
                <a href="https://github.com/galtay/ds2hf">ds2hf</a> pipeline.
                Source markdown lives in
                <code>dev_research/liu_2018/sections/</code>.
            </footer>
        </main>
    </div>
    <script>
{JS_SORTABLE}    </script>
</body>
</html>
"""
    OUT.write_text(html)
    print(
        f"Wrote {OUT.relative_to(HERE.parent.parent)}  ({len(files)} sections, {len(toc_entries)} TOC entries)"
    )
    print(f"  open: file://{OUT}")


if __name__ == "__main__":
    main()
