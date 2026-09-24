"""Load open-access GDC Pathology Report PDFs into per-case report records.

Unlike `mutations` / `expression`, there is nothing to parse here: the
payload is the source document. Each report is one scanned PDF of the
surgical pathology report for a case's tumor sample, shipped by GDC with
patient identifiers redacted out of the page image. We carry the bytes
verbatim in `pdf_bytes`.

## Why bytes and not text

These are page scans. Most carry an OCR text layer added upstream of GDC,
so `pypdf` / `pdfminer` / poppler all return several hundred to a few
thousand characters per report — but that layer transcribes the barcode
strip and handwritten margin notes as noise, and its fidelity varies by
submitting institution. Any extraction is therefore a tool- and
version-specific derivation, which is exactly the kind of thing this
project keeps re-derivable rather than baked in (cf. `survival.py`
re-deriving Liu's endpoints on every build instead of shipping his frozen
values). Shipping the bytes means a consumer can run a better parser next
year without re-downloading from GDC, and a canonical parse — if we add one
— becomes an additional clearly-labeled column rather than a replacement
for the source.

## Linkage

GDC attaches each report to a `sample` (not an aliquot), reported in the
file's `associated_entities`. The report's own UUID is the middle
component of the file name:

    TCGA-W5-AA2X.D2B18607-E16D-4570-9E96-5A7CBAFD79FC.PDF
    ^case_submitter_id ^pathology_report_uuid

which is the same value GDC reports on `sample.pathology_report_uuid` —
already present in this dataset's `SAMPLE_FIELDS`. So each report joins to
its sample two independent ways, and `attach` prefers the explicit
`associated_entities` pointer, falling back to the filename UUID.

Layout on disk:

    <data-dir>/raw/<project_id>/pathology_reports/
        TCGA-W5-AA2X.D2B18607-....PDF
        manifest.json
"""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path
from typing import Any

# GDC file names are `<case_submitter_id>.<REPORT_UUID>.PDF`. Splitting on
# "." is enough (submitter ids and UUIDs never contain one), but we still
# require the shape rather than assuming it, so a future rename surfaces as
# a null FK instead of a wrong one.
_UUID_HYPHEN_GROUPS = (8, 4, 4, 4, 12)


def _pathology_report_uuid(file_name: str) -> str | None:
    """Extract the report UUID from a GDC pathology report file name."""
    parts = file_name.split(".")
    if len(parts) < 3:
        return None
    candidate = parts[1]
    groups = candidate.split("-")
    if len(groups) != len(_UUID_HYPHEN_GROUPS):
        return None
    if any(len(g) != n for g, n in zip(groups, _UUID_HYPHEN_GROUPS, strict=True)):
        return None
    if not all(c in "0123456789abcdefABCDEF" for g in groups for c in g):
        return None
    return candidate


def _sample_from_entities(entry: dict[str, Any]) -> tuple[str | None, str | None]:
    """Return (sample_id, sample_submitter_id) from the manifest's associated_entities.

    Returns (None, None) when the file names no sample-level entity or names
    more than one — the caller then falls back to the filename UUID.
    """
    samples = [
        e
        for e in (entry.get("associated_entities") or [])
        if e.get("entity_type") == "sample" and e.get("entity_id")
    ]
    if len(samples) != 1:
        return None, None
    return samples[0]["entity_id"], samples[0].get("entity_submitter_id")


def _case_id(entry: dict[str, Any]) -> str | None:
    """The single case this report belongs to; None if absent or ambiguous."""
    case_ids = {c["case_id"] for c in (entry.get("cases") or []) if c.get("case_id")}
    if len(case_ids) != 1:
        return None
    return next(iter(case_ids))


def _build_record(
    entry: dict[str, Any],
    pdf_bytes: bytes,
) -> dict[str, Any]:
    sample_id, sample_submitter_id = _sample_from_entities(entry)
    return {
        "sample_id": sample_id,
        "sample_submitter_id": sample_submitter_id,
        "pathology_report_uuid": _pathology_report_uuid(entry["file_name"]),
        "source_file_id": entry["file_id"],
        "file_name": entry["file_name"],
        "file_size": entry.get("file_size"),
        "md5sum": entry.get("md5sum"),
        "pdf_bytes": pdf_bytes,
    }


def load_for_project(project_raw_dir: Path) -> dict[str, list[dict[str, Any]]]:
    """Return {case_id: [report_record, ...]} for raw/<PROJECT>/pathology_reports/.

    Empty dict if the directory or its manifest is missing — letting `build`
    proceed without pathology reports when they haven't been fetched.
    Manifest entries whose PDF isn't on disk (`_status="manifest_only"`) are
    skipped, matching how the other modalities treat a capped fetch.
    """
    reports_dir = project_raw_dir / "pathology_reports"
    manifest_path = reports_dir / "manifest.json"
    if not manifest_path.exists():
        return {}

    manifest = json.loads(manifest_path.read_text())
    by_case: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for entry in manifest:
        file_path = reports_dir / entry["file_name"]
        if not file_path.exists():
            continue
        case_id = _case_id(entry)
        if not case_id:
            continue
        by_case[case_id].append(_build_record(entry, file_path.read_bytes()))
    return dict(by_case)


def attach(
    rows: list[dict[str, Any]], by_case: dict[str, list[dict[str, Any]]]
) -> list[dict[str, Any]]:
    """Mutate `rows` to populate `samples_pathology_report`.

    Reports whose `associated_entities` didn't name a sample fall back to
    matching `pathology_report_uuid` against the patient's own
    `samples[].pathology_report_uuid`. Sorted by file_name for deterministic
    output; rows with no reports get [].
    """
    for row in rows:
        records = by_case.get(row["case_id"], [])
        uuid_to_sample: dict[str, dict[str, str | None]] = {}
        for s in row.get("samples") or []:
            report_uuid = s.get("pathology_report_uuid")
            if report_uuid:
                uuid_to_sample[report_uuid] = {
                    "sample_id": s.get("sample_id"),
                    "sample_submitter_id": s.get("submitter_id"),
                }
        for r in records:
            if r["sample_id"] is None:
                fallback = uuid_to_sample.get(r["pathology_report_uuid"] or "")
                if fallback:
                    r["sample_id"] = fallback["sample_id"]
                    r["sample_submitter_id"] = fallback["sample_submitter_id"]
        records.sort(key=lambda r: r.get("file_name") or "")
        row["samples_pathology_report"] = records
    return rows
