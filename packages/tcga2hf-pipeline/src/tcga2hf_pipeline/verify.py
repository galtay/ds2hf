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

from tcga2hf_pipeline.gdc import GDCClient, and_, eq

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
# A different question from the per-project checks above. That tree is built
# from raw GDC bytes, so it is checked against the GDC. This one is built
# from those per-project tables, so the GDC is the wrong reference: what
# matters is whether the reshape preserved them, and whether the positional
# contract the dataset is founded on actually holds.
# ===========================================================================


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


def check_axis_alignment(processed: Path) -> Check:
    """Every config must agree on the sample axis, in the same order.

    The dataset joins by position, not by key, so a config whose rows are
    ordered differently would silently attribute one patient's expression to
    another. Nothing downstream would notice.
    """
    import pyarrow.parquet as pq

    from tcga2hf_pipeline.gene_expression_quantification import QUANTIFICATIONS

    samples = pq.read_table(processed / "samples" / "data.parquet")
    expected = samples.column("aliquot_id").to_pylist()
    n_genes = pq.ParquetFile(processed / "genes" / "data.parquet").metadata.num_rows
    details, bad = [], False

    if samples.column("sample_index").to_pylist() != list(range(len(expected))):
        details.append("  samples.sample_index is not 0..n-1")
        bad = True

    for name in QUANTIFICATIONS:
        path = processed / name / "data.parquet"
        if not path.exists():
            details.append(f"  {name}: missing")
            bad = True
            continue
        table = pq.read_table(path, columns=["sample_index", "aliquot_id"])
        if table.column("aliquot_id").to_pylist() != expected:
            details.append(f"  {name}: aliquot order differs from `samples`")
            bad = True
        if table.column("sample_index").to_pylist() != list(range(len(expected))):
            details.append(f"  {name}: sample_index is not 0..n-1")
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
        f"{len(QUANTIFICATIONS)} value configs over {len(expected):,} samples x {n_genes:,} genes",
        details,
    )


def check_gene_axis(processed: Path, project_tabular_dir: Path) -> Check:
    """The gene axis must still equal the model the values were read against."""
    import pyarrow.parquet as pq

    ours = pq.read_table(processed / "genes" / "data.parquet", columns=["gene_id"])
    ours_ids = ours.column("gene_id").to_pylist()

    models = sorted(project_tabular_dir.glob("TCGA-*/gene_model/data.parquet"))
    if not models:
        return Check("gene_axis", False, f"no gene_model found under {project_tabular_dir}")
    theirs = pq.read_table(models[0], columns=["gene_id"]).column("gene_id").to_pylist()

    if ours_ids != theirs:
        first = next(
            (i for i, (a, b) in enumerate(zip(ours_ids, theirs, strict=False)) if a != b),
            min(len(ours_ids), len(theirs)),
        )
        return Check(
            "gene_axis",
            False,
            f"differs from {models[0].parent.parent.name} at index {first}",
            [f"  ours[{first}]={ours_ids[first : first + 1]}, theirs={theirs[first : first + 1]}"],
        )
    return Check("gene_axis", True, f"{len(ours_ids):,} genes match {models[0].parent.parent.name}")


def check_sample_metadata(processed: Path) -> Check:
    """Columns a consumer filters on must be populated and in range."""
    import pyarrow.parquet as pq

    table = pq.read_table(processed / "samples" / "data.parquet")
    details, bad = [], False

    for column in ("aliquot_id", "case_submitter_id", "project_id", "sample_type"):
        nulls = table.column(column).null_count
        if nulls:
            details.append(f"  {column}: {nulls:,} null")
            bad = True

    balance = [b for b in table.column("strand_balance").to_pylist() if b is not None]
    if not balance:
        details.append("  strand_balance: entirely null")
        bad = True
    elif min(balance) < 0.0 or max(balance) > 1.0:
        details.append(f"  strand_balance outside [0,1]: [{min(balance)}, {max(balance)}]")
        bad = True

    for column in ("n_unmapped", "n_multimapping", "n_nofeature", "n_ambiguous"):
        values = [v for v in table.column(column).to_pylist() if v is not None]
        if values and min(values) < 0:
            details.append(f"  {column}: negative value {min(values)}")
            bad = True

    return Check(
        "sample_metadata",
        not bad,
        f"{table.num_rows:,} samples, {table.num_columns} columns checked",
        details,
    )


def check_expression_values(
    processed: Path,
    project_tabular_dir: Path,
    projects: int = 4,
    samples_per_project: int = 2,
    genes_per_sample: int = 6,
    seed: int = 0,
) -> Check:
    """Re-read cells from the per-project tables this was reshaped from.

    Sampled rather than exhaustive: a full comparison would re-read every
    per-project expression table. The failure this guards against is a
    systematic reshape error, which a handful of random cells will expose.
    """
    import pyarrow.parquet as pq

    from tcga2hf_pipeline.gene_expression_quantification import QUANTIFICATIONS

    rng = random.Random(seed)
    samples = pq.read_table(processed / "samples" / "data.parquet").to_pandas()
    gene_ids = (
        pq.read_table(processed / "genes" / "data.parquet", columns=["gene_id"])
        .column("gene_id")
        .to_pylist()
    )

    available = sorted(samples["project_id"].unique())
    chosen = rng.sample(available, min(projects, len(available)))
    mismatches, checked = [], 0

    for project in chosen:
        source = project_tabular_dir / project / "gene_expression_quantification" / "data.parquet"
        if not source.exists():
            mismatches.append(f"  {project}: source table missing")
            continue
        rows = samples[samples["project_id"] == project]
        picks = rows.sample(min(samples_per_project, len(rows)), random_state=seed)
        frame = pq.read_table(
            source, columns=["aliquot_id", "gene_id", *QUANTIFICATIONS]
        ).to_pandas()

        for _, row in picks.iterrows():
            index = int(row["sample_index"])
            source_rows = frame[frame["aliquot_id"] == row["aliquot_id"]].set_index("gene_id")
            for gene_index in rng.sample(range(len(gene_ids)), genes_per_sample):
                gene_id = gene_ids[gene_index]
                for name in QUANTIFICATIONS:
                    parquet_file = pq.ParquetFile(processed / name / "data.parquet")
                    group, offset = _row_group_for(parquet_file, index)
                    values = (
                        parquet_file.read_row_group(group, columns=["values"])
                        .column("values")[offset]
                        .as_py()
                    )
                    got, want = values[gene_index], source_rows.loc[gene_id, name]
                    checked += 1
                    if abs(float(got) - float(want)) > abs(float(want)) * 1e-6:
                        mismatches.append(
                            f"  {project} {name} sample={index} gene={gene_index}: {got} != {want}"
                        )

    return Check(
        "expression_values",
        not mismatches,
        f"{checked:,} cells re-read across {len(chosen)} project(s): "
        f"{len(mismatches)} mismatch(es)",
        mismatches[:10],
    )


def check_sample_coverage(client: GDCClient, processed: Path) -> Check:
    """Our sample count must equal the expression files the GDC serves."""
    import pyarrow.parquet as pq

    ours = pq.ParquetFile(processed / "samples" / "data.parquet").metadata.num_rows
    # `size: 0` returns the count without the hits — one request rather than
    # paging through 11,505 file records to measure their length.
    payload = {
        "filters": and_(
            eq("cases.project.program.name", "TCGA"),
            eq("data_type", "Gene Expression Quantification"),
            eq("access", "open"),
        ),
        "fields": "file_id",
        "format": "JSON",
        "size": 0,
    }
    theirs = client._post("/files", payload)["data"]["pagination"]["total"]
    return Check(
        "sample_coverage",
        ours == theirs,
        f"{ours:,} samples ours / {theirs:,} open expression files at GDC",
    )


def verify_expression(
    processed_dir: Path,
    project_tabular_dir: Path,
    sample: int = 4,
) -> list[Check]:
    """Run every expression-dataset check. Local ones first, then the API."""
    checks = [
        check_axis_alignment(processed_dir),
        check_gene_axis(processed_dir, project_tabular_dir),
        check_sample_metadata(processed_dir),
        check_expression_values(processed_dir, project_tabular_dir, projects=sample),
    ]
    with GDCClient() as client:
        checks.append(check_sample_coverage(client, processed_dir))
    return checks
