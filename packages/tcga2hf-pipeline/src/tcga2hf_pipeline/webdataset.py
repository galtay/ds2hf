"""Row/member emitters for the WebDataset HF dataset (`tcga-wds-open`).

Third view of the same GDC source data, alongside the consolidated patient
Parquets (`tcga2hf_pipeline.clinical`) and the flat tabular tables
(`tcga2hf_pipeline.tabular`). Where those two ship *parsed* data, this one
ships the **files GDC serves**, grouped one sample per patient:

    data/<project_id>/<project>-<NNNNNN>.tar
        TCGA-3X-AAV9.case.json          this patient's GDC /cases record
        TCGA-3X-AAV9.files.jsonl        one line per member below
        TCGA-3X-AAV9.gene_expression_quantification.<file_id>.tsv.gz
        TCGA-3X-AAV9.masked_somatic_mutation.<file_id>.maf.gz
        TCGA-3X-AAV9.pathology_report.<file_id>.pdf
        ...
        TCGA-3X-AAVC.case.json
        ...

Naming rules, all of them GDC's rather than ours:

  - The sample **key** is `case.submitter_id` (`TCGA-3X-AAV9`). WebDataset
    splits a member name at the *first* dot, so the key must contain none —
    TCGA barcodes are hyphenated, so it never does.
  - The member **stem** is the snake_cased GDC `data_type`
    (`gene_expression_quantification`, not our `expression` raw-dir
    shorthand), matching the table names already published by
    `tabular.TABULAR_TABLES`.
  - Multiple files of one data_type per patient (copy number runs 2-7) are
    disambiguated by GDC `file_id`, the stable handle GDC versions against —
    never by a positional index of our own invention.
  - The suffix is the GDC file's own, with `.gz` appended where we compress.

Members are the GDC bytes verbatim, gzipped (`gzip_members=True`, the
default) except where GDC already ships them compressed — MAFs arrive as
`.maf.gz` and PDFs are internally compressed, so neither is re-gzipped.
`wds.WebDataset(...).decode()` and the HF `datasets` loader both gunzip
transparently. Because that applies to GDC's own gzip too, `files.jsonl`
carries `md5sum` (GDC's, over the file GDC serves) *and* `md5sum_member`
(over the exact bytes in the tar), with `gzipped_by_pipeline` naming the
transformation between them.

The one exception to raw-bytes is the pair of BCR biotab supplements: GDC
ships those per *project*, not per patient, so a per-patient member is a row
subset. `.clinical_supplement.` / `.biospecimen_supplement.` members keep the
biotab's 3-row header and only the patient's data rows, and `files.jsonl`
marks them `"subset_of_gdc_file": true` so the derivation is never implicit.
"""

from __future__ import annotations

import gzip
import hashlib
import io
import json
import tarfile
from collections import defaultdict
from pathlib import Path
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq

from tcga2hf_pipeline.biospecimen_supplement import _case_submitter_id as _bio_case_id

# Our raw-directory shorthand -> snake_cased GDC `data_type`. The values are
# the vocabulary already published by the tabular dataset; the keys exist only
# on local disk and never reach a shard.
RAW_DIR_TO_DATA_TYPE: dict[str, str] = {
    "expression": "gene_expression_quantification",
    "mutations": "masked_somatic_mutation",
    "mirna": "mirna_expression_quantification",
    "protein_expression": "protein_expression_quantification",
    "copy_number_masked": "masked_copy_number_segment",
    "copy_number_allele_specific": "allele_specific_copy_number_segment",
    "pathology_reports": "pathology_report",
}

# Project-scoped BCR biotabs. GDC `data_type` is "Clinical Supplement" /
# "Biospecimen Supplement"; GDC's own `type` enum matches these dir names.
SUPPLEMENT_DIRS: tuple[str, ...] = ("clinical_supplement", "biospecimen_supplement")

# GDC's `data_type` for the biotabs. Absent from our supplement manifests
# (the fetch doesn't request the field); confirmed against `GET /files`.
SUPPLEMENT_DATA_TYPE: dict[str, str] = {
    "clinical_supplement": "Clinical Supplement",
    "biospecimen_supplement": "Biospecimen Supplement",
}

# Manifest fields copied onto `files.jsonl`. `gdc_version` /
# `gdc_first_release` / `gdc_superseded` are the byte-provenance triple from
# `POST /files/versions` (see the pipeline README).
_FILE_FIELDS: tuple[str, ...] = (
    "file_id",
    "file_name",
    "file_size",
    "md5sum",
    "data_category",
    "data_type",
    "data_format",
    "experimental_strategy",
    "workflow_type",
    "access",
    "gdc_version",
    "gdc_first_release",
    "gdc_superseded",
)

# Compression suffixes GDC applies itself, so `_gdc_suffix` keeps them.
_COMPRESSION_SUFFIXES: frozenset[str] = frozenset({"gz", "bz2", "zip"})

# Formats we never re-compress: already compressed by GDC, or (PDF) an
# internally-compressed container that gzip cannot meaningfully shrink.
_ALREADY_COMPRESSED: frozenset[str] = frozenset({"gz", "bz2", "zip", "pdf"})

SHARD_TARGET_BYTES = 1_000_000_000


# ---------------------------------------------------------------------------
# Indexing raw/ by patient
# ---------------------------------------------------------------------------


def _gdc_suffix(file_name: str) -> str:
    """The GDC file's actual extension, lowercased.

    The last dot-component, plus the one before it when the last is a
    compression suffix GDC applied itself:

        `9a5e....rna_seq.augmented_star_gene_counts.tsv` -> `tsv`
        `00c2....wxs.aliquot_ensemble_masked.maf.gz`     -> `maf.gz`
        `BASIC_p_....nocnv_grch38.seg.v2.txt`            -> `txt`
        `TCGA-3X-AAV9.6E1F4753-....PDF`                  -> `pdf`

    Interior components (`rna_seq`, `seg`, `v2`) are GDC filename structure
    rather than extension, and dropping them keeps one field name per
    data_type. Nothing is lost: `files.jsonl` carries GDC's `file_name`
    verbatim for every member.
    """
    parts = file_name.lower().split(".")
    if len(parts) < 2:
        return "dat"
    if parts[-1] in _COMPRESSION_SUFFIXES and len(parts) > 2:
        return f"{parts[-2]}.{parts[-1]}"
    return parts[-1]


def index_project(project_raw_dir: Path) -> dict[str, dict[str, Any]]:
    """Map `case_submitter_id` -> per-patient file entries under one project.

    Mirrors `tabular._files_rows`: reads every modality `manifest.json` and
    explodes each entry's `cases` list. Entries whose bytes aren't on disk
    (manifest-only, i.e. discovered but never downloaded) are skipped.
    """
    by_case: dict[str, dict[str, Any]] = defaultdict(lambda: {"files": []})
    for raw_dir, data_type in RAW_DIR_TO_DATA_TYPE.items():
        modality_dir = project_raw_dir / raw_dir
        manifest_path = modality_dir / "manifest.json"
        if not manifest_path.exists():
            continue
        for entry in json.loads(manifest_path.read_text()):
            path = modality_dir / entry["file_name"]
            if not path.exists():
                continue
            for case in entry.get("cases") or []:
                submitter_id = case.get("submitter_id")
                if not submitter_id:
                    continue
                by_case[submitter_id]["files"].append(
                    {"path": path, "data_type": data_type, "entry": entry}
                )
    return by_case


def slice_supplements(project_raw_dir: Path) -> dict[str, list[dict[str, Any]]]:
    """Per-patient row subsets of the project-scoped BCR biotabs.

    GDC ships one set of forms per project, so this is the one derived
    member in a shard. Each slice keeps the biotab's 3-row header verbatim
    and appends only the rows belonging to that patient, so the result is
    still a well-formed biotab that parses with `parse_biotab`.
    """
    out: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for supp_dir_name in SUPPLEMENT_DIRS:
        supp_dir = project_raw_dir / supp_dir_name
        manifest_path = supp_dir / "manifest.json"
        if not manifest_path.exists():
            continue
        entries = {e["file_name"]: e for e in json.loads(manifest_path.read_text())}
        for path in sorted(supp_dir.glob("*.txt")):
            entry = entries.get(path.name)
            if entry is None:
                continue
            with path.open(newline="", encoding="latin-1") as handle:
                lines = handle.read().split("\n")
            if len(lines) < 4:
                continue
            header, body = lines[:3], lines[3:]
            columns = header[0].split("\t")
            rows_by_case: dict[str, list[str]] = defaultdict(list)
            for line in body:
                if not line.strip():
                    continue
                record = dict(zip(columns, line.split("\t"), strict=False))
                submitter_id = _bio_case_id(record)
                if submitter_id:
                    rows_by_case[submitter_id].append(line)
            for submitter_id, rows in rows_by_case.items():
                text = "\n".join(header + rows) + "\n"
                out[submitter_id].append(
                    {
                        "data_type": supp_dir_name,
                        "entry": entry,
                        "bytes": text.encode("latin-1"),
                        "n_rows": len(rows),
                    }
                )
    return out


# ---------------------------------------------------------------------------
# Shard writing
# ---------------------------------------------------------------------------


def _add(tar: tarfile.TarFile, name: str, data: bytes) -> None:
    """Append one member. `mtime=0` / fixed mode keep shards reproducible."""
    info = tarfile.TarInfo(name)
    info.size = len(data)
    info.mtime = 0
    info.mode = 0o444
    info.uid = info.gid = 0
    info.uname = info.gname = ""
    tar.addfile(info, io.BytesIO(data))


def _file_record(
    item: dict[str, Any],
    member: str,
    member_bytes: bytes,
    *,
    gzipped: bool,
    subset: bool = False,
) -> dict:
    """One `files.jsonl` record.

    Two checksums, because they answer different questions and conflating
    them is a trap:

      - `md5sum` is **GDC's**, over the file GDC serves. For a `.maf.gz`
        that file is already gzipped, so this is the md5 of compressed
        bytes; for a `.tsv` it is the md5 of plain text.
      - `md5sum_member` is over the exact bytes stored in the tar member,
        whatever we did to them. It verifies with no decode step at all.

    `gzipped_by_pipeline` says which transformation sits between the two.
    Note that `wds.decode()` gunzips *any* `.gz` member, GDC's own included,
    so a decoded MAF will match neither checksum until it is re-gzipped —
    read members undecoded when checksumming.
    """
    entry = item["entry"]
    record = {field: entry.get(field) for field in _FILE_FIELDS}
    record["member"] = member
    record["data_type_snake"] = item["data_type"]
    record["gzipped_by_pipeline"] = gzipped
    record["md5sum_member"] = hashlib.md5(member_bytes).hexdigest()
    if subset:
        # The member is a row subset of the named GDC file, not its bytes,
        # so GDC's checksum and size describe a file we did not ship.
        record["subset_of_gdc_file"] = True
        record["n_rows"] = item["n_rows"]
        record["md5sum"] = None
        record["file_size"] = None
        record["data_type"] = SUPPLEMENT_DATA_TYPE[item["data_type"]]
    else:
        record["subset_of_gdc_file"] = False
        record["n_rows"] = None
    return record


def write_sample(
    tar: tarfile.TarFile,
    key: str,
    case: dict[str, Any],
    items: list[dict[str, Any]],
    supplements: list[dict[str, Any]],
    *,
    gzip_members: bool = True,
) -> tuple[int, list[dict[str, Any]]]:
    """Write one patient's members, contiguously.

    Returns `(bytes_written, records)` where `records` is exactly what went
    into this sample's `files.jsonl` — the file-grain index is built from the
    same list, so the two cannot drift apart.

    WebDataset groups successive members sharing a prefix into one sample, so
    every member of a patient must be written before the next patient starts.
    """
    records: list[dict[str, Any]] = []
    payloads: list[tuple[str, bytes]] = []

    for item in sorted(items, key=lambda i: (i["data_type"], i["entry"]["file_id"])):
        entry = item["entry"]
        suffix = _gdc_suffix(entry["file_name"])
        data = item["path"].read_bytes()
        member = f"{key}.{item['data_type']}.{entry['file_id']}.{suffix}"
        gzipped = gzip_members and suffix.rsplit(".", 1)[-1] not in _ALREADY_COMPRESSED
        if gzipped:
            data = gzip.compress(data, 6, mtime=0)
            member += ".gz"
        records.append(_file_record(item, member, data, gzipped=gzipped))
        payloads.append((member, data))

    for item in sorted(supplements, key=lambda i: (i["data_type"], i["entry"]["file_id"])):
        entry = item["entry"]
        data = item["bytes"]
        member = f"{key}.{item['data_type']}.{entry['file_id']}.txt"
        if gzip_members:
            data = gzip.compress(data, 6, mtime=0)
            member += ".gz"
        records.append(_file_record(item, member, data, gzipped=gzip_members, subset=True))
        payloads.append((member, data))

    # Named for what each member actually holds. `case.json` is a single
    # object — this patient's `/cases` record — so it is singular and plain
    # JSON. `files.jsonl` is one object per line, one line per member, so it
    # is plural and JSON Lines. (The repo-root Parquet configs are `cases`
    # and `files`, plural in both cases: those are tables of many rows.)
    _add(tar, f"{key}.case.json", json.dumps(case, sort_keys=True).encode() + b"\n")
    _add(
        tar,
        f"{key}.files.jsonl",
        b"".join(json.dumps(r, sort_keys=True).encode() + b"\n" for r in records),
    )
    written = 0
    for member, data in payloads:
        _add(tar, member, data)
        written += len(data)
    return written, records


def build_project(
    cases: list[dict[str, Any]],
    project_raw_dir: Path,
    out_dir: Path,
    project_id: str,
    *,
    gzip_members: bool = True,
    shard_target_bytes: int = SHARD_TARGET_BYTES,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Write one project's shards. Returns `(index_rows, file_rows)`.

    `index_rows` is patient-grain (one per case), `file_rows` member-grain
    (one per tar member that holds a file). Both feed a Parquet config in the
    published repo.

    Patients are emitted in `case_submitter_id` order and a shard is closed
    once it passes `shard_target_bytes`, so a patient is never split across
    shards and a project never shares one with its neighbours. Every case in
    `cases.json` gets a sample, including the handful with no open-access
    files at all — those carry `case.json` and an empty `files.jsonl`.
    """
    by_case = index_project(project_raw_dir)
    supplements = slice_supplements(project_raw_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    index_rows: list[dict[str, Any]] = []
    file_rows: list[dict[str, Any]] = []
    shard_idx, written = 0, 0
    tar: tarfile.TarFile | None = None
    shard_name = ""

    def _open(idx: int) -> tuple[tarfile.TarFile, str]:
        name = f"{project_id.lower()}-{idx:06d}.tar"
        return tarfile.open(out_dir / name, "w", format=tarfile.GNU_FORMAT), name

    try:
        for case in sorted(cases, key=lambda c: c.get("submitter_id") or ""):
            key = case.get("submitter_id")
            if not key:
                continue
            items = by_case.get(key, {}).get("files", [])
            supps = supplements.get(key, [])
            if tar is None:
                tar, shard_name = _open(shard_idx)
            elif written >= shard_target_bytes:
                tar.close()
                shard_idx, written = shard_idx + 1, 0
                tar, shard_name = _open(shard_idx)
            sample_bytes, records = write_sample(
                tar, key, case, items, supps, gzip_members=gzip_members
            )
            written += sample_bytes

            index_rows.append(
                {
                    "case_id": case.get("case_id"),
                    "case_submitter_id": key,
                    "project_id": project_id,
                    "gdc_portal_url": (
                        f"https://portal.gdc.cancer.gov/cases/{case['case_id']}"
                        if case.get("case_id")
                        else None
                    ),
                    "shard": f"data/{project_id}/{shard_name}",
                    "n_files": len(items),
                    "n_bytes": sum(i["entry"].get("file_size") or 0 for i in items),
                }
            )
            patient_context = {
                "case_id": case.get("case_id"),
                "case_submitter_id": key,
                "project_id": project_id,
                "shard": f"data/{project_id}/{shard_name}",
            }
            # A row here is about a *file*, so its portal link points at the
            # file page. Join back to the case page through `case_id`.
            file_rows.extend(
                {
                    **patient_context,
                    **record,
                    "gdc_portal_url": (
                        f"https://portal.gdc.cancer.gov/files/{record['file_id']}"
                        if record.get("file_id")
                        else None
                    ),
                }
                for record in records
            )
    finally:
        if tar is not None:
            tar.close()
    return index_rows, file_rows


def cases_schema() -> pa.Schema:
    """Arrow schema for `cases.parquet` — all scalars, nothing to truncate.

    Deliberately narrow: identifiers, a portal link, the shard to fetch, and
    how much is in it. Clinical attributes are not copied up here — the
    patient's `case.json` member is the authority on those, and the
    `files` config covers anything at file grain.
    """
    return pa.schema(
        [
            pa.field("case_id", pa.string()),
            pa.field("case_submitter_id", pa.string()),
            pa.field("project_id", pa.string()),
            pa.field("gdc_portal_url", pa.string()),
            pa.field("shard", pa.string()),
            # Counts of GDC files only; supplement slices are derived and not
            # counted here. Which data_types a patient actually has lives in
            # that patient's `files.jsonl`, not denormalised into the index.
            pa.field("n_files", pa.int64()),
            pa.field("n_bytes", pa.int64()),
        ]
    )


def write_cases_index(rows: list[dict[str, Any]], processed_dir: Path) -> Path:
    """Write the repo-root `cases.parquet` (patient-grain config).

    Keeping the shards *out* of the declared configs is deliberate: HF
    `datasets` resolves one builder module per repo, so a card declaring a
    Parquet config alongside a tar config makes it read the shards as Parquet
    and fail. Two Parquet configs coexist fine; the shards stay undeclared and
    `wds.WebDataset` reads them straight off their resolve URLs.
    """
    schema = cases_schema()
    table = pa.Table.from_pylist(
        sorted(rows, key=lambda r: (r["project_id"], r["case_submitter_id"] or "")),
        schema=schema,
    )
    out = processed_dir / "cases.parquet"
    pq.write_table(table, out, compression="zstd")
    return out


def files_schema() -> pa.Schema:
    """Arrow schema for `files.parquet` — one row per shard member.

    A complete table of contents for the tars: every member is findable here,
    with `subset_of_gdc_file` separating the derived supplement slices from
    the files GDC serves. Column order follows `files.jsonl`, prefixed with
    the patient context and the shard the member lives in.
    """
    return pa.schema(
        [
            pa.field("case_id", pa.string()),
            pa.field("case_submitter_id", pa.string()),
            pa.field("project_id", pa.string()),
            # The file's own portal page. For a supplement slice this links
            # the whole-project biotab the slice was cut from.
            pa.field("gdc_portal_url", pa.string()),
            pa.field("shard", pa.string()),
            pa.field("member", pa.string()),
            pa.field("data_type_snake", pa.string()),
            pa.field("file_id", pa.string()),
            pa.field("file_name", pa.string()),
            # Null for a supplement slice: GDC's size and checksum describe
            # the whole-project form, which is not what the member holds.
            pa.field("file_size", pa.int64()),
            pa.field("md5sum", pa.string()),
            pa.field("md5sum_member", pa.string()),
            pa.field("gzipped_by_pipeline", pa.bool_()),
            pa.field("subset_of_gdc_file", pa.bool_()),
            pa.field("n_rows", pa.int64()),
            pa.field("data_category", pa.string()),
            pa.field("data_type", pa.string()),
            pa.field("data_format", pa.string()),
            pa.field("experimental_strategy", pa.string()),
            pa.field("workflow_type", pa.string()),
            pa.field("access", pa.string()),
            pa.field("gdc_version", pa.string()),
            pa.field("gdc_first_release", pa.string()),
            pa.field("gdc_superseded", pa.bool_()),
        ]
    )


def write_files_index(rows: list[dict[str, Any]], processed_dir: Path) -> Path:
    """Write the repo-root `files.parquet` (member-grain config)."""
    schema = files_schema()
    # Project onto the schema: `rows` also carry keys the index doesn't ship.
    table = pa.Table.from_pylist(
        [
            {name: row.get(name) for name in schema.names}
            for row in sorted(
                rows, key=lambda r: (r["project_id"], r["case_submitter_id"] or "", r["member"])
            )
        ],
        schema=schema,
    )
    out = processed_dir / "files.parquet"
    pq.write_table(table, out, compression="zstd")
    return out
