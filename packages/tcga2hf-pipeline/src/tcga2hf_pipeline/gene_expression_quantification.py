"""Builder for the expression-only HF dataset (`tcga-gene-expression-quantification-open`).

The per-project datasets serve `gene_expression_quantification` as one row
per (aliquot, gene) — the shape that keeps each project standalone and
joinable against the rest of its cohort. Transcriptomics wants the other
shape: the whole cohort as one matrix, sample-major, so that a batch of
samples is a batch of rows.

That is what this module writes. One row per aliquot, the 60,660 gene
values as a single list column, and the gene axis pinned by position
rather than by a join key:

    genes[i]                <-> values[i]      for every value config
    aliquots[j]             <-> row j          within a value config

which is exactly AnnData's `var` / `obs` / `X` split. Reading the three
configs back into an `AnnData` takes about a second for the full cohort;
see the dataset card for the four lines that do it.

The row is an aliquot, not a sample: GDC publishes one STAR file per
aliquot, and a sample sequenced twice (TCGA-GBM has many) has two.

One config per GDC quantification, each named for the column it carries in
GDC's own STAR-counts TSV (`unstranded`, `tpm_unstranded`, ...) rather than
a name of our invention. A consumer downloads the one they model on and
none of the others.

Two narrowings, both verified lossless against what the GDC publishes:
counts are `int32` (max observed 4.5M against int32's 2.1B) and the
normalized values are `float32` (round-trip error 6e-8 against the <=4
decimal places GDC prints). Nothing here re-derives a value — every number
is the one GDC serves, and the source TSV row is recoverable by pairing
`genes` with a row's `values`.

Gene *selection* is deliberately absent. All 60,660 GENCODE v36 features
ship, with `gene_type` on the `genes` config, so restricting by biotype
stays the consumer's choice rather than ours.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

from tcga2hf_pipeline import expression as _expression_mod

# One config per column of the GDC STAR-counts TSV, GDC's spelling kept.
# The dtype is the narrowest that holds every published value exactly.
QUANTIFICATIONS: dict[str, pa.DataType] = {
    "unstranded": pa.int32(),
    "stranded_first": pa.int32(),
    "stranded_second": pa.int32(),
    "tpm_unstranded": pa.float32(),
    "fpkm_unstranded": pa.float32(),
    "fpkm_uq_unstranded": pa.float32(),
}

# Sample-major rows are ~200 KB each, so a row group of 64 is ~13 MB — small
# enough to stream a minibatch without decoding the file, large enough that
# the per-row-group footer stays negligible. Measured flat on both size and
# throughput from 32 to 512, so this is a comfort choice, not a tuned one.
#
# Note this is a *maximum*, not a stride: each project is written by its own
# `write_table` call, which starts fresh row groups, so a 79-sample project
# contributes groups of 64 and 15. Row-group boundaries therefore do not sit
# on a global multiple of 64, and locating sample `i` by `i // ROW_GROUP_SIZE`
# is wrong — read the `aliquot_index` column instead.
ROW_GROUP_SIZE = 64

# Label columns duplicated onto every value config. They cost ~40 bytes per
# row and remove the obs join from the training loop entirely; the full
# record still lives once in `aliquots`.
INLINE_ALIQUOT_FIELDS: list[pa.Field] = [
    pa.field("aliquot_index", pa.int32()),
    pa.field("aliquot_id", pa.string()),
    pa.field("case_submitter_id", pa.string()),
    pa.field("project_id", pa.string()),
    pa.field("sample_type", pa.string()),
]

ALIQUOTS_FIELDS: list[pa.Field] = [
    pa.field("aliquot_index", pa.int32()),
    pa.field("aliquot_id", pa.string()),
    pa.field("aliquot_submitter_id", pa.string()),
    # The aliquot's parents in GDC's biospecimen graph, innermost first:
    # aliquot -> analyte (the extracted RNA) -> portion -> sample -> case.
    # Read from the graph, not parsed out of the barcode: for a few dozen
    # aliquots the recorded parent's barcode is not a prefix of the
    # aliquot's own, so barcode slicing would name the wrong parent.
    pa.field("analyte_id", pa.string()),
    pa.field("analyte_submitter_id", pa.string()),
    pa.field("analyte_type", pa.string()),
    pa.field("portion_id", pa.string()),
    pa.field("portion_submitter_id", pa.string()),
    pa.field("sample_id", pa.string()),
    pa.field("sample_submitter_id", pa.string()),
    pa.field("sample_type", pa.string()),
    pa.field("case_id", pa.string()),
    pa.field("case_submitter_id", pa.string()),
    pa.field("project_id", pa.string()),
    # The GDC file the row was read from. GDC serves only its current
    # release, so the release we downloaded at does not identify the bytes;
    # `file_id` + `md5sum` do. `version` and `first_release` come from
    # `/files/versions`, and the URL is templated from the id: it serves
    # the exact TSV, so any row can be re-derived and checked against it.
    pa.field("source_file_id", pa.string()),
    pa.field("source_file_md5sum", pa.string()),
    pa.field("source_file_version", pa.string()),
    pa.field("source_file_first_release", pa.string()),
    pa.field("source_file_url", pa.string()),
    # Derived here, not from GDC: stranded_first / (stranded_first +
    # stranded_second) over the whole library. ~0.5 means the library is not
    # strand-specific and `unstranded` is the column to use; values near 0
    # mean a reverse-stranded (dUTP) protocol, where `stranded_second`
    # carries the signal. Most of TCGA sits at 0.5, but a real minority does
    # not — chiefly TCGA-GBM — so this is published per library rather than
    # asserted in prose. Null when the library has no stranded reads.
    pa.field("strand_balance", pa.float32()),
    # STAR's four unassigned-read tallies, copied from the pseudo-gene rows
    # at the top of the source TSV. Together with the gene counts they
    # account for every read in the library, so `assigned / total` is a
    # per-library QC measure — and it is not a constant: it tracks the
    # project (TCGA-LAML's median is lowest), which makes it a confounder
    # worth being able to condition on rather than discover. The tallies
    # differ by strand column in the TSV; these are the `unstranded` ones.
    # The range the card quotes is measured at build (`assigned_summary`).
    pa.field("n_unmapped", pa.int64()),
    pa.field("n_multimapping", pa.int64()),
    pa.field("n_nofeature", pa.int64()),
    pa.field("n_ambiguous", pa.int64()),
]

# Source TSV name -> our column. Snake-cased the way every other GDC column
# in this project is; `N_noFeature` loses its interior capital.
QC_COLUMNS: dict[str, str] = {
    "N_unmapped": "n_unmapped",
    "N_multimapping": "n_multimapping",
    "N_noFeature": "n_nofeature",
    "N_ambiguous": "n_ambiguous",
}

# The tallies are the one thing here read from `raw/` rather than from the
# per-project tables, which carry exactly the 60,660 gene rows. That costs
# nothing: the N_* rows sit at lines 3-6 of each TSV, ahead of the gene
# rows, so only the head of each file is read (~0.7 s for 500 files) rather
# than the ~49 GB the modality occupies. They stay off the matrix and on the
# aliquot, which is also why every `values` list is exactly 60,660 long and
# needs no masking.

GENES_FIELDS: list[pa.Field] = [
    pa.field("gene_index", pa.int32()),
    pa.field("gene_id", pa.string()),
    pa.field("gene_name", pa.string()),
    pa.field("gene_type", pa.string()),
    pa.field("chromosome", pa.string()),
    pa.field("start", pa.int64()),
    pa.field("end", pa.int64()),
]


def _values_column(matrix: np.ndarray, value_type: pa.DataType) -> pa.Array:
    """A list<value_type> column, one list per row of `matrix`.

    Built from the flat buffer with an explicit offsets array rather than
    from a list of rows: `pa.array(list(matrix))` boxes every element and is
    orders of magnitude slower at 60,660 values a row.
    """
    n_rows, width = matrix.shape
    offsets = pa.array(np.arange(n_rows + 1, dtype=np.int32) * width, type=pa.int32())
    flat = pa.array(matrix.reshape(-1), type=value_type)
    return pa.ListArray.from_arrays(offsets, flat)


def gene_axis(project_dir: Path) -> pa.Table:
    """The `genes` config: GDC's gene model with its array position.

    Read from a project's `gene_model` table, which is the GENCODE v36
    model GDC repeats identically in every per-gene file. `gene_index` is
    just the row number, made explicit so the positional contract with
    `values` is visible in the data rather than only in prose.
    """
    # Kept in Arrow rather than round-tripped through pandas: `gene_type`
    # and `chromosome` carry real nulls (the 37 chrM genes have no
    # chromosome), and pandas would turn those into NaN floats.
    model = pq.read_table(project_dir / "gene_model" / "data.parquet")
    columns = {"gene_index": pa.array(np.arange(model.num_rows, dtype=np.int32))}
    for field in GENES_FIELDS[1:]:
        columns[field.name] = model.column(field.name).combine_chunks().cast(field.type)
    return pa.Table.from_pydict(columns, schema=pa.schema(GENES_FIELDS))


def qc_tallies(project_raw_dir: Path) -> dict[str, dict[str, int | None]]:
    """aliquot_id -> STAR's four unassigned-read tallies, from the raw TSVs.

    Reads only the head of each file. The N_* rows are lines 3-6 — after the
    `# gene-model` comment and the header, before the first ENSG row — so
    this stops at the first gene row and never touches the other 60,660
    lines. Parsing is positional (`gene_id`, `gene_name`, `gene_type`,
    `unstranded`); the two name columns are empty on these rows.

    Missing files or missing rows yield nulls rather than zeros: a library
    with no unmapped reads and a library we failed to read are different
    facts, and only one of them is a measurement.
    """
    expr_dir = project_raw_dir / "expression"
    manifest_path = expr_dir / "manifest.json"
    if not manifest_path.exists():
        return {}

    out: dict[str, dict[str, int | None]] = {}
    for entry in json.loads(manifest_path.read_text()):
        path = expr_dir / entry.get("file_name", "")
        if not path.exists():
            continue
        _, aliquot_id = _expression_mod._file_aliquot_and_case(entry)
        if not aliquot_id:
            continue
        tallies: dict[str, int | None] = dict.fromkeys(QC_COLUMNS.values())
        with path.open() as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                name = parts[0]
                if name in QC_COLUMNS:
                    value = parts[3] if len(parts) > 3 else ""
                    tallies[QC_COLUMNS[name]] = int(value) if value else None
                elif name != "gene_id":
                    break  # first gene row: the QC block is behind us
        out[aliquot_id] = tallies
    return out


def _biospecimen_index(cases_path: Path) -> dict[str, dict[str, Any]]:
    """aliquot_id -> every parent's id and submitter id, for one project.

    Expression rows carry case and aliquot FKs only, so the levels between
    them come back through the biospecimen tree in `cases.json`, taking
    each parent as the node the aliquot is actually nested under.
    """
    out: dict[str, dict[str, Any]] = {}
    for case in json.loads(cases_path.read_text()):
        for sample in case.get("samples") or []:
            for portion in sample.get("portions") or []:
                for analyte in portion.get("analytes") or []:
                    for aliquot in analyte.get("aliquots") or []:
                        aid = aliquot.get("aliquot_id")
                        if not aid:
                            continue
                        out[aid] = {
                            "aliquot_submitter_id": aliquot.get("submitter_id"),
                            "analyte_id": analyte.get("analyte_id"),
                            "analyte_submitter_id": analyte.get("submitter_id"),
                            "analyte_type": analyte.get("analyte_type"),
                            "portion_id": portion.get("portion_id"),
                            "portion_submitter_id": portion.get("submitter_id"),
                            "sample_id": sample.get("sample_id"),
                            "sample_submitter_id": sample.get("submitter_id"),
                            "sample_type": sample.get("sample_type"),
                        }
    return out


GDC_DATA_URL = "https://api.gdc.cancer.gov/data/{file_id}"


def source_files(project_raw_dir: Path) -> dict[str, dict[str, Any]]:
    """source file_id -> its GDC provenance, from the fetch manifest.

    The manifest recorded what `/files` and `/files/versions` said about
    each file when it was downloaded; nothing here re-queries the GDC.
    """
    manifest_path = project_raw_dir / "expression" / "manifest.json"
    if not manifest_path.exists():
        return {}
    return {
        entry["file_id"]: {
            "source_file_md5sum": entry.get("md5sum"),
            "source_file_version": entry.get("gdc_version"),
            "source_file_first_release": entry.get("gdc_first_release"),
            "source_file_url": GDC_DATA_URL.format(file_id=entry["file_id"]),
        }
        for entry in json.loads(manifest_path.read_text())
    }


def download_release(raw_dir: Path, projects: list[str]) -> dict[str, Any]:
    """The GDC release the expression files were downloaded at.

    Read from the status the expression fetch wrote alongside its files,
    not the project-level one, which belongs to the clinical fetch and can
    predate or postdate it. The card states one release, so projects
    fetched at different releases are an error rather than a footnote.
    """
    releases: dict[str, dict[str, Any]] = {}
    for project in projects:
        path = raw_dir / project / "expression" / "gdc_status.json"
        if not path.exists():
            raise FileNotFoundError(f"{project}: no expression fetch status at {path}")
        releases[project] = json.loads(path.read_text())
    distinct = {r["data_release"] for r in releases.values()}
    if len(distinct) != 1:
        by_release: dict[str, list[str]] = {}
        for project, status in releases.items():
            by_release.setdefault(status["data_release"], []).append(project)
        raise ValueError(f"expression was fetched at more than one GDC release: {by_release}")
    return next(iter(releases.values()))


def _project_block(
    expr_path: Path,
    cases_path: Path,
    project_id: str,
    n_genes: int,
) -> tuple[list[dict[str, Any]], dict[str, np.ndarray]]:
    """One project's aliquot rows plus its {quantification: (n, n_genes)} block.

    Read once and reshaped, rather than once per quantification: the source
    table is the largest thing this build touches and six passes over it
    would dominate the runtime.
    """
    columns = [
        "aliquot_id",
        "aliquot_submitter_id",
        "case_id",
        "case_submitter_id",
        "source_file_id",
        "gene_id",
        *QUANTIFICATIONS,
    ]
    df = pq.read_table(expr_path, columns=columns).to_pandas()

    genes = df
    del df

    aliquots = genes["aliquot_id"].drop_duplicates().to_numpy()
    n = len(aliquots)
    if n and len(genes) != n * n_genes:
        raise ValueError(
            f"{project_id}: expected {n} aliquots x {n_genes} genes = "
            f"{n * n_genes:,} rows, found {len(genes):,}"
        )

    blocks = {
        name: genes[name].to_numpy().reshape(n, n_genes).astype(dtype.to_pandas_dtype())
        for name, dtype in QUANTIFICATIONS.items()
    }

    biospecimen = _biospecimen_index(cases_path)
    tallies = qc_tallies(cases_path.parent)
    files = source_files(cases_path.parent)
    first = genes.groupby("aliquot_id", sort=False).first()
    rows = []
    for aliquot in aliquots:
        rec = first.loc[aliquot]
        meta = biospecimen.get(aliquot, {})
        # Every row must name the bytes it came from; a row whose file the
        # manifest does not know would publish an unverifiable value.
        provenance = files.get(rec["source_file_id"])
        if provenance is None:
            raise ValueError(
                f"{project_id}: {aliquot} was read from {rec['source_file_id']}, "
                "which is not in the expression manifest"
            )
        rows.append(
            {
                "aliquot_id": aliquot,
                "aliquot_submitter_id": meta.get("aliquot_submitter_id")
                or rec["aliquot_submitter_id"],
                "analyte_id": meta.get("analyte_id"),
                "analyte_submitter_id": meta.get("analyte_submitter_id"),
                "analyte_type": meta.get("analyte_type"),
                "portion_id": meta.get("portion_id"),
                "portion_submitter_id": meta.get("portion_submitter_id"),
                "sample_id": meta.get("sample_id"),
                "sample_submitter_id": meta.get("sample_submitter_id"),
                "sample_type": meta.get("sample_type"),
                "case_id": rec["case_id"],
                "case_submitter_id": rec["case_submitter_id"],
                "project_id": project_id,
                "source_file_id": rec["source_file_id"],
                **provenance,
                **tallies.get(aliquot, dict.fromkeys(QC_COLUMNS.values())),
            }
        )
    return rows, blocks


def strand_summary(first: np.ndarray, second: np.ndarray) -> dict[str, float]:
    """Whether the cohort's libraries are strand-specific, measured.

    STAR emits three count columns because it cannot know the library
    protocol: `unstranded`, and one for each strand assignment. Picking the
    right one is the analyst's job, and picking a stranded column for an
    unstranded library yields noise rather than signal.

    The tell is `first / (first + second)` per sample. A strand-specific
    protocol drives it towards 0 or 1; an unstranded library leaves it at
    0.5, because the reads are just being split by whichever strand they
    happened to align to. Reported so the card can state what the shipped
    data actually is instead of asserting it.
    """
    both = first + second
    balance = np.divide(first, both, out=np.full_like(first, np.nan), where=both > 0)
    balance = balance[np.isfinite(balance)]
    return {
        "n_samples": int(balance.size),
        "balance_median": float(np.median(balance)) if balance.size else float("nan"),
        "balance_p01": float(np.percentile(balance, 1)) if balance.size else float("nan"),
        "balance_p99": float(np.percentile(balance, 99)) if balance.size else float("nan"),
        # Samples whose balance leaves the band an unstranded library sits in.
        "n_strand_specific": int(((balance < 0.4) | (balance > 0.6)).sum()),
    }


def assigned_summary(assigned: np.ndarray, rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Range of the fraction of each library's reads STAR assigned to a gene.

    `assigned` is the row sum of the `unstranded` counts. The four `N_*`
    tallies stored on the row are the ones from that same column, so
    `assigned + sum(N_*)` is every read in the library. Rows without
    tallies are left out rather than treated as fully assigned.
    """
    unassigned = np.array(
        [
            np.nan
            if any(r.get(c) is None for c in QC_COLUMNS.values())
            else sum(r[c] for c in QC_COLUMNS.values())
            for r in rows
        ],
        dtype=np.float64,
    )
    total = assigned + unassigned
    fraction = np.divide(assigned, total, out=np.full_like(total, np.nan), where=total > 0)
    fraction = fraction[np.isfinite(fraction)]
    if not fraction.size:
        return {"assigned_n": 0, "assigned_min": float("nan"), "assigned_max": float("nan")}
    return {
        "assigned_n": int(fraction.size),
        "assigned_min": float(fraction.min()),
        "assigned_max": float(fraction.max()),
    }


def build(
    processed_project_dir: Path,
    raw_dir: Path,
    out_dir: Path,
    projects: list[str] | None = None,
    echo: Any = None,
) -> tuple[dict[str, int], dict[str, Any]]:
    """Write the expression dataset tree, one config per directory.

    Streams project by project into six open `ParquetWriter`s rather than
    accumulating the cohort: aliquots are disjoint across projects, so each
    project's block is simply the next set of row groups. Peak memory is
    one project's six matrices (~1.8 GB for TCGA-BRCA, the largest), not
    the ~17 GB the full cohort would need.

    Returns `({config_name: rows written}, stats)`. The stats are measured
    here, while the matrices are already in memory, so the card can
    describe the data it actually ships rather than quote a remembered
    constant: strandedness (see `strand_summary`) and the range of the
    fraction of reads assigned to genes.
    """
    say = echo or (lambda _m: None)
    dirs = sorted(
        d
        for d in processed_project_dir.glob("TCGA-*")
        if (d / "gene_expression_quantification" / "data.parquet").exists()
        and (projects is None or d.name in projects)
    )
    if not dirs:
        raise ValueError(f"no built projects with expression under {processed_project_dir}")

    genes_table = gene_axis(dirs[0])
    n_genes = genes_table.num_rows
    gene_ids = genes_table.column("gene_id").to_pylist()

    out_dir.mkdir(parents=True, exist_ok=True)
    writers: dict[str, pq.ParquetWriter] = {}
    schemas = {
        name: pa.schema([*INLINE_ALIQUOT_FIELDS, pa.field("values", pa.list_(dtype))])
        for name, dtype in QUANTIFICATIONS.items()
    }
    for name, schema in schemas.items():
        path = out_dir / name / "data.parquet"
        path.parent.mkdir(parents=True, exist_ok=True)
        writers[name] = pq.ParquetWriter(path, schema, compression="zstd", write_page_index=True)

    aliquot_rows: list[dict[str, Any]] = []
    strand_totals: dict[str, list[np.ndarray]] = {"first": [], "second": [], "assigned": []}
    try:
        for project_dir in dirs:
            project_id = project_dir.name
            cases_path = raw_dir / project_id / "cases.json"
            if not cases_path.exists():
                raise FileNotFoundError(f"{project_id}: no cases.json at {cases_path}")

            # Every project must repeat the same gene model in the same
            # order, or a positional gene axis would silently misalign.
            other = (
                pq.read_table(project_dir / "gene_model" / "data.parquet", columns=["gene_id"])
                .column("gene_id")
                .to_pylist()
            )
            if other != gene_ids:
                raise ValueError(
                    f"{project_id}: gene_model differs from {dirs[0].name}; "
                    "a positional gene axis requires an identical model"
                )

            rows, blocks = _project_block(
                project_dir / "gene_expression_quantification" / "data.parquet",
                cases_path,
                project_id,
                n_genes,
            )
            base = len(aliquot_rows)
            for offset, row in enumerate(rows):
                row["aliquot_index"] = base + offset
            aliquot_rows.extend(rows)

            inline = {
                "aliquot_index": np.arange(base, base + len(rows), dtype=np.int32),
                "aliquot_id": [r["aliquot_id"] for r in rows],
                "case_submitter_id": [r["case_submitter_id"] for r in rows],
                "project_id": [r["project_id"] for r in rows],
                "sample_type": [r["sample_type"] for r in rows],
            }
            for name, dtype in QUANTIFICATIONS.items():
                table = pa.Table.from_pydict(
                    {**inline, "values": _values_column(blocks[name], dtype)},
                    schema=schemas[name],
                )
                writers[name].write_table(table, row_group_size=ROW_GROUP_SIZE)
            # Per-sample library totals, kept while the matrices are here.
            # Three float64 scalars a sample; the arrays are freed below.
            strand_totals["first"].append(blocks["stranded_first"].sum(axis=1, dtype=np.float64))
            strand_totals["second"].append(blocks["stranded_second"].sum(axis=1, dtype=np.float64))
            strand_totals["assigned"].append(blocks["unstranded"].sum(axis=1, dtype=np.float64))
            del blocks
            say(f"  {project_id:<12}{len(rows):>6} aliquots  (total {len(aliquot_rows):,})")
    finally:
        for writer in writers.values():
            writer.close()

    all_first = np.concatenate(strand_totals["first"])
    all_second = np.concatenate(strand_totals["second"])
    both = all_first + all_second
    balances = np.divide(all_first, both, out=np.full_like(all_first, np.nan), where=both > 0)
    for row, balance in zip(aliquot_rows, balances, strict=True):
        row["strand_balance"] = None if not np.isfinite(balance) else float(balance)

    counts = {name: len(aliquot_rows) for name in QUANTIFICATIONS}
    stats: dict[str, Any] = strand_summary(all_first, all_second)
    stats.update(assigned_summary(np.concatenate(strand_totals["assigned"]), aliquot_rows))

    (out_dir / "genes").mkdir(parents=True, exist_ok=True)
    pq.write_table(
        genes_table,
        out_dir / "genes" / "data.parquet",
        compression="zstd",
        write_page_index=True,
    )
    counts["genes"] = n_genes

    (out_dir / "aliquots").mkdir(parents=True, exist_ok=True)
    pq.write_table(
        pa.Table.from_pylist(aliquot_rows, schema=pa.schema(ALIQUOTS_FIELDS)),
        out_dir / "aliquots" / "data.parquet",
        compression="zstd",
        write_page_index=True,
    )
    counts["aliquots"] = len(aliquot_rows)
    return counts, stats
