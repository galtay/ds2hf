"""Check a built project dataset against the GDC, not against our own beliefs.

Every other test in this repo asserts that the pipeline does what the
pipeline's author intended. These checks ask a different question: does the
published tree agree with what the GDC actually serves *right now*? They
therefore talk to the live API and re-hash local bytes rather than trusting
any manifest we wrote.

Five checks, cheapest first:

  1. `file_coverage`   — every open-access GDC file for the project is either
                         in our `files` table or on the documented exclusion
                         list. Catches a modality silently not fetched.
  2. `manifest_vs_gdc` — our recorded `file_size` / `md5sum` match what the
                         API reports today. Catches a file GDC has revised
                         under the same id.
  3. `local_md5`       — raw bytes on disk hash to GDC's `md5sum`. Catches a
                         truncated or corrupted download. Sampled.
  4. `case_coverage`   — the `cases` table's case_id set equals the GDC's
                         for the project. Catches a paging or merge bug.
  5. `fk_integrity`    — every `aliquot_id` in every molecular table resolves
                         inside that patient's own biospecimen tree. Purely
                         local, but it is the join every consumer relies on.

Each returns a `Check`; `verify_project` runs them and the CLI renders them.
Nothing here imports from the build path, so a bug shared with the builder
cannot hide itself from these.
"""

from __future__ import annotations

import hashlib
import json
import random
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from ds2hf_pipeline.tcga.gdc import GDCClient, and_, eq

# Open-access data types we deliberately do not fetch. Listed here rather
# than inferred so that a modality vanishing from the pipeline shows up as a
# failure instead of quietly joining this list.
EXPECTED_EXCLUSIONS: dict[str, str] = {
    "Slide Image": "whole-slide .svs images; ~89 GiB for one project, not tabular",
    "Masked Intensities": "raw .idat arrays; the SeSAMe betas are the analysis-ready form",
}


@dataclass
class Check:
    name: str
    passed: bool
    summary: str
    details: list[str] = field(default_factory=list)


def _gdc_open_file_counts(client: GDCClient, project: str) -> dict[str, int]:
    """{data_type: count} for every open-access file GDC holds for the project."""
    payload = {
        "filters": and_(
            eq("cases.project.project_id", project),
            eq("access", "open"),
        ),
        "facets": "data_type",
        "size": 0,
        "format": "JSON",
    }
    data = client._post("/files", payload)["data"]
    return {b["key"]: b["doc_count"] for b in data["aggregations"]["data_type"]["buckets"]}


def _our_files(processed_project_dir: Path) -> list[dict[str, Any]]:
    import pyarrow.parquet as pq

    path = processed_project_dir / "files" / "data.parquet"
    if not path.exists():
        return []
    table = pq.read_table(
        path,
        columns=[
            "file_id",
            "file_name",
            "file_size",
            "md5sum",
            "data_type",
            "modality",
            "in_dataset",
            "gdc_download_url",
        ],
    )
    seen: set[str] = set()
    out = []
    for row in table.to_pylist():
        if row["file_id"] in seen:
            continue
        seen.add(row["file_id"])
        out.append(row)
    return out


def check_file_coverage(client: GDCClient, project: str, processed: Path) -> Check:
    """Every open GDC file has a row, and the `in_dataset` flag is honest.

    Since the `files` table indexes the project's whole open footprint, row
    parity is now exact — a missing row means the index is stale, not that a
    modality was skipped. What varies is `in_dataset`, so that is checked
    separately: nothing on the exclusion list may claim to be published.
    """
    gdc = _gdc_open_file_counts(client, project)
    rows = _our_files(processed)
    ours: dict[str, int] = {}
    published: dict[str, int] = {}
    for row in rows:
        dt = row.get("data_type")
        ours[dt] = ours.get(dt, 0) + 1
        if row.get("in_dataset"):
            published[dt] = published.get(dt, 0) + 1

    details, bad = [], False
    for data_type, expected in sorted(gdc.items(), key=lambda kv: -kv[1]):
        got = ours.get(data_type, 0)
        n_pub = published.get(data_type, 0)
        if got != expected:
            bad = True
        excluded = data_type in EXPECTED_EXCLUSIONS
        if excluded and n_pub:
            bad = True
            note = f"MISMATCH {n_pub} rows claim in_dataset but this type is excluded"
        elif got != expected:
            note = "MISMATCH"
        else:
            note = f"{n_pub} in dataset" + (" (excluded by design)" if excluded else "")
        details.append(f"  {data_type:<38}{expected:>6} gdc {got:>6} rows  {note}")

    missing_url = sum(1 for r in rows if not r.get("gdc_download_url"))
    if missing_url:
        bad = True
        details.append(f"  {missing_url} row(s) lack a gdc_download_url")
    nulls = ours.get(None, 0)
    if nulls:
        bad = True
        details.append(f"  {nulls} files carry a null data_type")

    total = sum(gdc.values())
    n_pub_total = sum(published.values())
    return Check(
        "file_coverage",
        not bad,
        f"{len(rows)}/{total} open GDC files indexed, {n_pub_total} carried in tables",
        details,
    )


def check_manifest_vs_gdc(client: GDCClient, project: str, processed: Path) -> Check:
    """Our recorded size/md5 must equal what the API reports for the same id."""
    # Only rows whose bytes we hold: an indexed-but-absent file has no
    # local copy whose size/md5 could have drifted.
    ours = {r["file_id"]: r for r in _our_files(processed) if r.get("modality")}
    if not ours:
        return Check("manifest_vs_gdc", False, "no downloaded files to check")

    remote: dict[str, dict[str, Any]] = {}
    ids = list(ours)
    for i in range(0, len(ids), 500):
        chunk = ids[i : i + 500]
        payload = {
            "filters": {"op": "in", "content": {"field": "file_id", "value": chunk}},
            "fields": "file_id,file_size,md5sum",
            "format": "JSON",
            "size": len(chunk),
        }
        for hit in client._post("/files", payload)["data"]["hits"]:
            remote[hit["file_id"]] = hit

    missing = [f for f in ours if f not in remote]
    size_bad, md5_bad = [], []
    for fid, row in ours.items():
        r = remote.get(fid)
        if not r:
            continue
        if row.get("file_size") is not None and r.get("file_size") != row["file_size"]:
            size_bad.append(fid)
        if row.get("md5sum") and r.get("md5sum") != row["md5sum"]:
            md5_bad.append(fid)

    details = []
    if missing:
        details.append(f"  {len(missing)} file_id(s) no longer resolve at GDC: {missing[:3]}")
    if size_bad:
        details.append(f"  {len(size_bad)} file_size mismatch(es): {size_bad[:3]}")
    if md5_bad:
        details.append(f"  {len(md5_bad)} md5sum mismatch(es): {md5_bad[:3]}")
    ok = not (missing or size_bad or md5_bad)
    return Check(
        "manifest_vs_gdc",
        ok,
        f"{len(ours)} downloaded files checked against the live API"
        + ("" if ok else f"; {len(missing)} missing, {len(size_bad)} size, {len(md5_bad)} md5"),
        details,
    )


def check_local_md5(project_raw_dir: Path, sample: int, seed: int = 0) -> Check:
    """Re-hash raw bytes on disk against the md5 GDC recorded in the manifest.

    Sampled per modality so one cheap run still touches every modality; a
    corrupt download in a rarely-read modality is exactly what this is for.
    """
    rng = random.Random(seed)
    checked = failed = 0
    details = []
    for manifest_path in sorted(project_raw_dir.glob("*/manifest.json")):
        entries = [
            e
            for e in json.loads(manifest_path.read_text())
            if e.get("md5sum") and (manifest_path.parent / e["file_name"]).exists()
        ]
        if not entries:
            continue
        picks = entries if len(entries) <= sample else rng.sample(entries, sample)
        for entry in picks:
            path = manifest_path.parent / entry["file_name"]
            digest = hashlib.md5()  # noqa: S324 - matching GDC's own checksum
            with path.open("rb") as fh:
                for block in iter(lambda: fh.read(1 << 20), b""):
                    digest.update(block)
            checked += 1
            if digest.hexdigest() != entry["md5sum"]:
                failed += 1
                details.append(f"  {manifest_path.parent.name}/{entry['file_name']}: md5 differs")
    return Check(
        "local_md5",
        failed == 0,
        f"{checked} raw files re-hashed, {failed} mismatched",
        details,
    )


def check_case_coverage(client: GDCClient, project: str, processed: Path) -> Check:
    import pyarrow.parquet as pq

    path = processed / "cases" / "data.parquet"
    if not path.exists():
        return Check("case_coverage", False, "no cases table")
    ours = set(pq.read_table(path, columns=["case_id"]).column("case_id").to_pylist())

    payload = {
        "filters": eq("project.project_id", project),
        "fields": "case_id",
        "format": "JSON",
        "size": 10000,
    }
    remote = {h["case_id"] for h in client._post("/cases", payload)["data"]["hits"]}

    missing, extra = remote - ours, ours - remote
    details = []
    if missing:
        details.append(f"  {len(missing)} GDC case(s) absent from our table: {sorted(missing)[:3]}")
    if extra:
        details.append(f"  {len(extra)} case(s) in our table but not at GDC: {sorted(extra)[:3]}")
    return Check(
        "case_coverage",
        not (missing or extra),
        f"{len(ours)} cases ours / {len(remote)} at GDC",
        details,
    )


# Molecular tables keyed by aliquot. Two are deliberately absent:
# `protein_expression_quantification` attaches to a *portion*, and
# `masked_somatic_mutation` carries `tumor_sample_id` rather than an
# aliquot. The check re-confirms the column exists before reading, so this
# tuple staying in sync is a convenience, not a correctness requirement.
_ALIQUOT_TABLES = (
    "gene_expression_quantification",
    "mirna_expression_quantification",
    "isoform_expression_quantification",
    "methylation_beta_value",
    "allele_specific_copy_number_segment",
    "masked_copy_number_segment",
    "copy_number_segment",
    "gene_level_copy_number",
)


def check_fk_integrity(processed: Path) -> Check:
    """Every molecular `aliquot_id` must exist in that case's biospecimen tree."""
    import pyarrow.parquet as pq

    cases_path = processed / "cases" / "data.parquet"
    if not cases_path.exists():
        return Check("fk_integrity", False, "no cases table")
    known: set[str] = set()
    for row in pq.read_table(cases_path, columns=["samples"]).to_pylist():
        for sample in row["samples"] or []:
            for portion in sample.get("portions") or []:
                for analyte in portion.get("analytes") or []:
                    for aliquot in analyte.get("aliquots") or []:
                        if aliquot.get("aliquot_id"):
                            known.add(aliquot["aliquot_id"])

    details, bad = [], 0
    checked = 0
    for table in _ALIQUOT_TABLES:
        path = processed / table / "data.parquet"
        if not path.exists():
            continue
        if "aliquot_id" not in pq.ParquetFile(path).schema_arrow.names:
            details.append(f"  {table}: no aliquot_id column, skipped")
            continue
        col = pq.read_table(path, columns=["aliquot_id"]).column("aliquot_id")
        ids = {v for v in col.to_pylist() if v}
        checked += 1
        orphans = ids - known
        if orphans:
            bad += len(orphans)
            details.append(
                f"  {table}: {len(orphans)} aliquot(s) not in any case tree, "
                f"e.g. {sorted(orphans)[:2]}"
            )
    return Check(
        "fk_integrity",
        bad == 0,
        f"{checked} aliquot-keyed tables checked against {len(known)} known aliquots, "
        f"{bad} orphan(s)",
        details,
    )


def verify_project(
    project: str,
    raw_dir: Path,
    processed_dir: Path,
    sample: int = 3,
) -> list[Check]:
    """Run every check for one project. Network checks share one client."""
    project_raw = raw_dir / project
    processed = processed_dir / project
    checks = [check_fk_integrity(processed), check_local_md5(project_raw, sample)]
    with GDCClient() as client:
        checks.insert(0, check_file_coverage(client, project, processed))
        checks.insert(1, check_manifest_vs_gdc(client, project, processed))
        checks.insert(2, check_case_coverage(client, project, processed))
    return checks


# ===========================================================================
# Expression dataset
#
# Checked against the GDC, not against the per-project tables it was
# reshaped from, so that the dataset's claims stand on its own sources: every
# row names a GDC file, and the checks below re-read those files -- from the
# md5-verified bytes on disk, from `source_file_url`, and from the live API --
# and compare them with what ships. The positional contract (genes[i] <->
# values[i], aliquots[j] <-> row j) is checked first, because every other
# check addresses cells by position and would be meaningless without it.
# ===========================================================================

_EXPR_ID_COLUMNS = (
    "aliquot_id",
    "aliquot_submitter_id",
    "analyte_id",
    "analyte_submitter_id",
    "analyte_type",
    "portion_id",
    "portion_submitter_id",
    "sample_id",
    "sample_submitter_id",
    "sample_type",
    "case_id",
    "case_submitter_id",
    "project_id",
    "source_file_id",
    "source_file_md5sum",
    "source_file_version",
    "source_file_first_release",
    "source_file_url",
)
_EXPR_TALLIES = {
    "N_unmapped": "n_unmapped",
    "N_multimapping": "n_multimapping",
    "N_noFeature": "n_nofeature",
    "N_ambiguous": "n_ambiguous",
}
_GDC_DATA_URL = "https://api.gdc.cancer.gov/data/"


def _row_group_for(parquet_file: Any, row_index: int) -> tuple[int, int]:
    """Locate a global row index as (row_group, offset within that group).

    Row groups are per-project write batches capped at ROW_GROUP_SIZE rather
    than a fixed stride, so `row_index // ROW_GROUP_SIZE` lands on the wrong
    row as soon as any project's sample count is not a multiple of it. Walk
    the metadata instead.
    """
    start = 0
    for group in range(parquet_file.metadata.num_row_groups):
        rows = parquet_file.metadata.row_group(group).num_rows
        if start + rows > row_index:
            return group, row_index - start
        start += rows
    raise IndexError(f"row {row_index} beyond {start} rows")


def _row_values(processed: Path, row_index: int) -> dict[str, Any]:
    """{quantification: numpy array} for one row, read from each value config."""
    import numpy as np
    import pyarrow.parquet as pq

    from ds2hf_pipeline.tcga.gene_expression_quantification import QUANTIFICATIONS

    out = {}
    for name in QUANTIFICATIONS:
        parquet_file = pq.ParquetFile(processed / name / "data.parquet")
        group, offset = _row_group_for(parquet_file, row_index)
        column = parquet_file.read_row_group(group, columns=["values"]).column("values")
        out[name] = np.asarray(column[offset].values)
    return out


def _read_star_tsv(raw: bytes) -> tuple[Any, Any]:
    """(tally rows, gene rows) of a GDC STAR-Counts TSV, as DataFrames."""
    import io

    import pandas as pd

    frame = pd.read_csv(io.BytesIO(raw), sep="\t", comment="#")
    is_tally = frame["gene_id"].str.startswith("N_")
    return frame[is_tally].set_index("gene_id"), frame[~is_tally].reset_index(drop=True)


def _compare_row_to_source(
    raw: bytes,
    row: dict[str, Any],
    values: dict[str, Any],
    gene_ids: list[str],
) -> list[str]:
    """Every way one published row can disagree with the file it names.

    Exact comparisons throughout. The counts are integers, and the float32
    values are compared with the TSV's text parsed and narrowed the same
    way, so any difference at all is a real one.
    """
    import numpy as np

    label = f"{row['project_id']} {row['aliquot_id']}"
    digest = hashlib.md5(raw).hexdigest()
    if digest != row["source_file_md5sum"]:
        return [f"  {label}: md5 {digest} != {row['source_file_md5sum']}"]

    tallies, genes = _read_star_tsv(raw)
    problems = []
    if genes["gene_id"].tolist() != gene_ids:
        problems.append(f"  {label}: gene order differs from `genes`")
        return problems

    for name, got in values.items():
        want = genes[name].to_numpy()
        want = want.astype(got.dtype) if got.dtype.kind == "f" else want.astype(np.int64)
        if not np.array_equal(got.astype(want.dtype), want):
            n = int((got.astype(want.dtype) != want).sum())
            problems.append(f"  {label}: {name} differs in {n:,} of {len(want):,} genes")

    for tsv_name, column in _EXPR_TALLIES.items():
        want = int(tallies.loc[tsv_name, "unstranded"])
        if row[column] != want:
            problems.append(f"  {label}: {column} {row[column]} != {want}")

    first = float(genes["stranded_first"].sum())
    second = float(genes["stranded_second"].sum())
    if first + second > 0:
        balance = first / (first + second)
        if row["strand_balance"] is None or abs(row["strand_balance"] - balance) > 1e-6:
            problems.append(f"  {label}: strand_balance {row['strand_balance']} != {balance}")
    return problems


def check_axis_alignment(processed: Path) -> Check:
    """Every config must agree on the aliquot axis, in the same order.

    The dataset joins by position, not by key, so a config whose rows are
    ordered differently would silently attribute one patient's expression to
    another. Nothing downstream would notice.
    """
    import pyarrow.parquet as pq

    from ds2hf_pipeline.tcga.gene_expression_quantification import QUANTIFICATIONS

    aliquots = pq.read_table(processed / "aliquots" / "data.parquet")
    expected = aliquots.column("aliquot_id").to_pylist()
    n_genes = pq.ParquetFile(processed / "genes" / "data.parquet").metadata.num_rows
    details, bad = [], False

    if aliquots.column("aliquot_index").to_pylist() != list(range(len(expected))):
        details.append("  aliquots.aliquot_index is not 0..n-1")
        bad = True

    for name in QUANTIFICATIONS:
        path = processed / name / "data.parquet"
        if not path.exists():
            details.append(f"  {name}: missing")
            bad = True
            continue
        table = pq.read_table(path, columns=["aliquot_index", "aliquot_id"])
        if table.column("aliquot_id").to_pylist() != expected:
            details.append(f"  {name}: aliquot order differs from `aliquots`")
            bad = True
        if table.column("aliquot_index").to_pylist() != list(range(len(expected))):
            details.append(f"  {name}: aliquot_index is not 0..n-1")
            bad = True
        # One row's list length stands for the config: a short list would
        # shift every gene after it.
        parquet_file = pq.ParquetFile(path)
        first = parquet_file.read_row_group(0, columns=["values"]).column("values")
        width = len(first[0].as_py())
        if width != n_genes:
            details.append(f"  {name}: values length {width:,} != {n_genes:,} genes")
            bad = True

    return Check(
        "axis_alignment",
        not bad,
        f"{len(QUANTIFICATIONS)} value configs over {len(expected):,} aliquots x {n_genes:,} genes",
        details,
    )


def check_gene_axis(processed: Path, raw_dir: Path) -> Check:
    """`genes` must equal the model in the GDC files it was read from.

    `gene_id` / `gene_name` / `gene_type` against a raw STAR TSV, and the
    coordinates against a raw gene-level copy number TSV, the two sources
    the card names.
    """
    import pandas as pd
    import pyarrow.parquet as pq

    ours = pq.read_table(processed / "genes" / "data.parquet").to_pandas()
    star = next(iter(sorted(raw_dir.glob("TCGA-*/expression/*.tsv"))), None)
    cnv = next(iter(sorted(raw_dir.glob("TCGA-*/gene_level_copy_number/*.tsv"))), None)
    if star is None or cnv is None:
        return Check("gene_axis", False, f"no raw STAR or gene-level CN TSV under {raw_dir}")

    details = []
    _, genes = _read_star_tsv(star.read_bytes())
    if genes["gene_id"].tolist() != ours["gene_id"].tolist():
        details.append(f"  gene_id order differs from {star.name}")
    else:
        for column in ("gene_name", "gene_type"):
            want = [None if pd.isna(v) else v for v in genes[column]]
            got = [None if pd.isna(v) else v for v in ours[column]]
            n = sum(a != b for a, b in zip(got, want, strict=True))
            if n:
                details.append(f"  {column}: {n:,} genes differ from {star.name}")

    coords = pd.read_csv(cnv, sep="\t", usecols=["gene_id", "chromosome", "start", "end"])
    merged = ours.merge(coords, on="gene_id", how="left", suffixes=("", "_cnv"))
    for column in ("chromosome", "start", "end"):
        present = merged[f"{column}_cnv"].notna()
        n = int((merged.loc[present, column] != merged.loc[present, f"{column}_cnv"]).sum())
        if n:
            details.append(f"  {column}: {n:,} genes differ from {cnv.name}")
        stray = int((merged[column].notna() & ~present).sum())
        if stray:
            details.append(f"  {column}: {stray:,} genes have a value the CN file does not")

    return Check(
        "gene_axis",
        not details,
        f"{len(ours):,} genes against {star.parent.parent.name} STAR + gene-level CN TSVs",
        details,
    )


def check_aliquot_metadata(processed: Path) -> Check:
    """Every row must be identified, unique, and name a well-formed source."""
    import re

    import pyarrow.parquet as pq

    frame = pq.read_table(processed / "aliquots" / "data.parquet").to_pandas()
    details = []

    for column in _EXPR_ID_COLUMNS:
        nulls = int(frame[column].isna().sum())
        if nulls:
            details.append(f"  {column}: {nulls:,} null")
    for column in ("aliquot_id", "source_file_id"):
        dupes = int(frame[column].duplicated().sum())
        if dupes:
            details.append(f"  {column}: {dupes:,} duplicated")

    md5 = re.compile(r"^[0-9a-f]{32}$")
    bad_md5 = int((~frame["source_file_md5sum"].fillna("").str.match(md5)).sum())
    if bad_md5:
        details.append(f"  source_file_md5sum: {bad_md5:,} not a 32-hex md5")
    bad_url = int((frame["source_file_url"] != _GDC_DATA_URL + frame["source_file_id"]).sum())
    if bad_url:
        details.append(f"  source_file_url: {bad_url:,} not the data URL of source_file_id")

    balance = frame["strand_balance"].dropna()
    if balance.empty:
        details.append("  strand_balance: entirely null")
    elif balance.min() < 0.0 or balance.max() > 1.0:
        details.append(f"  strand_balance outside [0,1]: [{balance.min()}, {balance.max()}]")
    for column in _EXPR_TALLIES.values():
        if (frame[column].dropna() < 0).any():
            details.append(f"  {column}: negative value")

    return Check(
        "aliquot_metadata",
        not details,
        f"{len(frame):,} aliquots, {len(_EXPR_ID_COLUMNS)} identifier columns checked",
        details,
    )


def _pick_rows(frame: Any, per_project: int, seed: int) -> Any:
    """`per_project` random rows from every project, so no project goes unread."""
    index = [
        i
        for _, group in frame.groupby("project_id")
        for i in group.sample(min(per_project, len(group)), random_state=seed).index
    ]
    return frame.loc[index]


def check_source_bytes(
    processed: Path, raw_dir: Path, per_project: int = 1, seed: int = 0
) -> Check:
    """Re-read whole rows from the md5-verified GDC files on disk.

    For each sampled row: the raw file it names must hash to its
    `source_file_md5sum`, and every one of its 60,660 x 6 values, its four
    read tallies and its `strand_balance` must equal what that file says.
    """
    import pyarrow.parquet as pq

    frame = pq.read_table(processed / "aliquots" / "data.parquet").to_pandas()
    gene_ids = pq.read_table(processed / "genes" / "data.parquet", columns=["gene_id"])
    gene_ids = gene_ids.column("gene_id").to_pylist()

    names: dict[str, Path] = {}
    for manifest in raw_dir.glob("TCGA-*/expression/manifest.json"):
        for entry in json.loads(manifest.read_text()):
            names[entry["file_id"]] = manifest.parent / entry["file_name"]

    picks = _pick_rows(frame, per_project, seed)
    problems: list[str] = []
    for row in picks.to_dict("records"):
        path = names.get(row["source_file_id"])
        if path is None or not path.exists():
            problems.append(f"  {row['project_id']} {row['aliquot_id']}: raw file missing")
            continue
        problems += _compare_row_to_source(
            path.read_bytes(), row, _row_values(processed, int(row["aliquot_index"])), gene_ids
        )

    return Check(
        "source_bytes",
        not problems,
        f"{len(picks):,} whole rows re-read from raw GDC files across "
        f"{picks['project_id'].nunique()} projects: {len(problems)} problem(s)",
        problems[:10],
    )


def check_source_urls(processed: Path, n: int = 2, seed: int = 0) -> Check:
    """Fetch rows' `source_file_url` and compare, as a consumer would.

    The same comparison as `check_source_bytes`, but through the URL the
    dataset publishes rather than a local copy, so it also proves the URL
    serves the file the row was read from.
    """
    import httpx
    import pyarrow.parquet as pq

    frame = pq.read_table(processed / "aliquots" / "data.parquet").to_pandas()
    gene_ids = pq.read_table(processed / "genes" / "data.parquet", columns=["gene_id"])
    gene_ids = gene_ids.column("gene_id").to_pylist()

    picks = frame.sample(min(n, len(frame)), random_state=seed)
    problems: list[str] = []
    with httpx.Client(timeout=120.0, follow_redirects=True) as http:
        for row in picks.to_dict("records"):
            response = http.get(row["source_file_url"])
            if response.status_code != 200:
                problems.append(f"  {row['source_file_url']}: HTTP {response.status_code}")
                continue
            problems += _compare_row_to_source(
                response.content, row, _row_values(processed, int(row["aliquot_index"])), gene_ids
            )

    return Check(
        "source_urls",
        not problems,
        f"{len(picks)} row(s) fetched from source_file_url and re-read: {len(problems)} problem(s)",
        problems[:10],
    )


def check_gdc_current(client: GDCClient, processed: Path) -> Check:
    """Is the source file set still exactly what the GDC serves?

    Compares by identity, not count: the set of file ids, and each file's
    md5 and version, against a live `/files` query for the dataset's scope.
    Then recomputes the source digest from both sides. Equal digests mean a
    rebuild today would read identical bytes; this is the check to run after
    a GDC release before deciding whether to rebuild.
    """
    import pyarrow.parquet as pq

    ours = pq.read_table(
        processed / "aliquots" / "data.parquet",
        columns=["source_file_id", "source_file_md5sum", "source_file_version"],
    ).to_pandas()
    hits = client.files(
        filters=and_(
            eq("cases.project.program.name", "TCGA"),
            eq("data_type", "Gene Expression Quantification"),
            eq("analysis.workflow_type", "STAR - Counts"),
            eq("access", "open"),
        ),
        fields=["file_id", "md5sum", "version"],
        page_size=500,
    )
    theirs = {h["file_id"]: h for h in hits}
    mine = {r.source_file_id: r for r in ours.itertuples(index=False)}

    added = sorted(set(theirs) - set(mine))
    withdrawn = sorted(set(mine) - set(theirs))
    changed = [
        f
        for f in sorted(set(mine) & set(theirs))
        if theirs[f]["md5sum"] != mine[f].source_file_md5sum
        or str(theirs[f].get("version")) != mine[f].source_file_version
    ]

    def digest(pairs: list[tuple[str, str]]) -> str:
        lines = sorted(f"{f}\t{m}\n" for f, m in pairs)
        return "sha256:" + hashlib.sha256("".join(lines).encode()).hexdigest()

    ours_digest = digest([(f, r.source_file_md5sum) for f, r in mine.items()])
    gdc_digest = digest([(f, h["md5sum"]) for f, h in theirs.items()])

    details = [f"  ours {ours_digest}", f"  gdc  {gdc_digest}"]
    details += [f"  added at GDC: {f}" for f in added[:5]]
    details += [f"  withdrawn from GDC: {f}" for f in withdrawn[:5]]
    details += [f"  md5 or version changed: {f}" for f in changed[:5]]
    return Check(
        "gdc_current",
        ours_digest == gdc_digest and not changed,
        f"{len(mine):,} ours / {len(theirs):,} at GDC: {len(added)} added, "
        f"{len(withdrawn)} withdrawn, {len(changed)} changed",
        details,
    )


def verify_expression(
    processed_dir: Path,
    raw_dir: Path,
    per_project: int = 1,
    remote: int = 2,
) -> list[Check]:
    """Run every expression-dataset check. Local ones first, then the network."""
    checks = [
        check_axis_alignment(processed_dir),
        check_gene_axis(processed_dir, raw_dir),
        check_aliquot_metadata(processed_dir),
        check_source_bytes(processed_dir, raw_dir, per_project=per_project),
    ]
    if remote:
        checks.append(check_source_urls(processed_dir, n=remote))
    with GDCClient() as client:
        checks.append(check_gdc_current(client, processed_dir))
    return checks


# ---------------------------------------------------------------------------
# Liu et al. 2018 CDR dataset.
#
# The source is one frozen workbook, so verification is about identity:
# the bytes shipped are the bytes GDC serves, and every cell of every
# config is the cell in that workbook. The workbook is re-read here with
# pandas rather than the builder's openpyxl walk, so a bug in one reader
# does not vouch for itself.
# ---------------------------------------------------------------------------

# GDC's open-access manifest for the PanCanAtlas publication page. It is the
# only place GDC states this file's md5: the file is served by `/data` but
# not indexed by `/files`.
PANCAN_OPEN_MANIFEST_URL = (
    "https://gdc.cancer.gov/system/files/public/file/PanCan-General_Open_GDC-Manifest_2.txt"
)

# A TCGA patient barcode: project site code, then participant code.
_PATIENT_BARCODE = r"^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}$"


def check_cdr_source_pinned(processed: Path) -> Check:
    """The shipped workbook hashes to the md5 the builder pinned."""
    from ds2hf_pipeline.tcga.cdr import CDR_FILE_MD5, CDR_FILE_NAME

    path = processed / "source" / CDR_FILE_NAME
    if not path.exists():
        return Check("source_pinned", False, f"{path} missing")
    md5 = hashlib.md5(path.read_bytes()).hexdigest()
    return Check("source_pinned", md5 == CDR_FILE_MD5, f"shipped {md5}, pinned {CDR_FILE_MD5}")


def check_cdr_source_at_gdc(processed: Path) -> Check:
    """GDC still lists and serves these exact bytes.

    Two independent statements from GDC: the md5 in its open-access
    manifest, and the md5 of what `/data/<uuid>` returns today without a
    token. Both must equal the shipped workbook's.
    """
    import httpx

    from ds2hf_pipeline.tcga.cdr import CDR_FILE_NAME, CDR_FILE_UUID, CDR_SOURCE_URL

    shipped = hashlib.md5((processed / "source" / CDR_FILE_NAME).read_bytes()).hexdigest()
    manifest = httpx.get(PANCAN_OPEN_MANIFEST_URL, timeout=60.0, follow_redirects=True)
    manifest.raise_for_status()
    listed = {
        fields[0]: fields
        for fields in (line.split("\t") for line in manifest.text.splitlines()[1:])
        if fields and fields[0]
    }
    entry = listed.get(CDR_FILE_UUID)
    manifest_md5 = entry[2] if entry and len(entry) > 2 else None

    served = httpx.get(CDR_SOURCE_URL, timeout=120.0, follow_redirects=True)
    served_md5 = hashlib.md5(served.content).hexdigest() if served.status_code == 200 else None

    details = [
        f"  shipped   {shipped}",
        f"  manifest  {manifest_md5} ({PANCAN_OPEN_MANIFEST_URL})",
        f"  served    {served_md5} (HTTP {served.status_code}, {CDR_SOURCE_URL})",
    ]
    return Check(
        "source_at_gdc",
        shipped == manifest_md5 == served_md5,
        "manifest and /data both match" if shipped == manifest_md5 == served_md5 else "mismatch",
        details,
    )


def check_cdr_cells(processed: Path) -> Check:
    """Every config cell equals its workbook cell; nothing added or dropped.

    The one allowed difference is Excel's `#N/A` error cell as null; pandas
    reads error cells as NaN, and with default NA parsing off nothing else
    becomes NaN. Column names and order must equal the sheet's header after
    its unlabelled row-number column, and row order must equal the sheet's.
    """
    import math

    import pandas as pd
    import pyarrow.parquet as pq

    from ds2hf_pipeline.tcga.cdr import CDR_CONFIGS, CDR_FILE_NAME

    workbook = processed / "source" / CDR_FILE_NAME
    details, compared = [], 0
    for config, (sheet, _) in CDR_CONFIGS.items():
        # dtype=object and no NA parsing: every cell arrives as Excel stored it.
        frame = pd.read_excel(workbook, sheet_name=sheet, dtype=object, keep_default_na=False)
        expected = frame.iloc[:, 1:]
        ours = pq.read_table(processed / config / "data.parquet")
        if list(expected.columns) != ours.column_names:
            details.append(f"  {config}: columns differ from sheet {sheet!r}")
            continue
        if len(expected) != ours.num_rows:
            details.append(f"  {config}: {ours.num_rows} rows, sheet has {len(expected)}")
            continue
        for name in ours.column_names:
            want = [
                None if isinstance(v, float) and math.isnan(v) else v
                for v in expected[name].tolist()
            ]
            got = ours.column(name).to_pylist()
            bad = [i for i, (w, g) in enumerate(zip(want, got, strict=True)) if w != g]
            compared += len(got)
            if bad:
                i = bad[0]
                details.append(
                    f"  {config}.{name}: {len(bad)} cell(s) differ, e.g. row {i}: "
                    f"sheet {want[i]!r} vs ours {got[i]!r}"
                )
    return Check(
        "cells_match_workbook",
        not details,
        f"{compared:,} cells compared across {len(CDR_CONFIGS)} configs",
        details,
    )


def check_cdr_keys(processed: Path) -> Check:
    """One row per patient, a well-formed barcode, the same patients in every config.

    Row *i* of every config is the same patient: the configs share row
    order as Liu published them, which is what makes a positional or a
    keyed join equally safe.
    """
    import re

    import pyarrow.parquet as pq

    from ds2hf_pipeline.tcga.cdr import CDR_CONFIGS

    keys = {
        config: pq.read_table(processed / config / "data.parquet", columns=["bcr_patient_barcode"])
        .column("bcr_patient_barcode")
        .to_pylist()
        for config in CDR_CONFIGS
    }
    details = []
    for config, barcodes in keys.items():
        malformed = [b for b in barcodes if b is None or not re.match(_PATIENT_BARCODE, b)]
        if malformed:
            details.append(f"  {config}: {len(malformed)} malformed, e.g. {malformed[:3]}")
        if len(set(barcodes)) != len(barcodes):
            details.append(f"  {config}: {len(barcodes) - len(set(barcodes))} duplicate barcode(s)")
    first, *rest = keys
    for other in rest:
        if keys[other] != keys[first]:
            details.append(f"  {other}: patients or their order differ from {first}")
    return Check(
        "patient_keys",
        not details,
        f"{len(keys[first]):,} patients, unique and aligned across {len(keys)} configs",
        details,
    )


def check_cdr_cases_at_gdc(client: GDCClient, processed: Path) -> Check:
    """Every patient is a TCGA case GDC serves today, by `submitter_id`.

    This is the join the card promises. A barcode GDC no longer serves
    would make that promise false for its row, so the check fails rather
    than reports: the card, not just the data, would need a look.
    """
    import pyarrow.parquet as pq

    barcodes = set(
        pq.read_table(processed / "cdr" / "data.parquet", columns=["bcr_patient_barcode"])
        .column("bcr_patient_barcode")
        .to_pylist()
    )
    hits = client.cases(
        filters=eq("project.program.name", "TCGA"), fields=["submitter_id"], page_size=5000
    )
    at_gdc = {h["submitter_id"] for h in hits}
    missing = sorted(barcodes - at_gdc)
    return Check(
        "cases_at_gdc",
        not missing,
        f"{len(barcodes) - len(missing):,} of {len(barcodes):,} patients are GDC TCGA cases",
        [f"  not at GDC: {b}" for b in missing[:10]],
    )


def verify_cdr(processed_dir: Path) -> list[Check]:
    """Run every CDR-dataset check. Local ones first, then the network."""
    checks = [
        check_cdr_source_pinned(processed_dir),
        check_cdr_cells(processed_dir),
        check_cdr_keys(processed_dir),
        check_cdr_source_at_gdc(processed_dir),
    ]
    with GDCClient() as client:
        checks.append(check_cdr_cases_at_gdc(client, processed_dir))
    return checks
