"""Clinical Supplement biotab fetcher + parser.

The harmonized GDC `/cases?expand=...` API drops or under-populates a
handful of clinical fields that Liu et al. 2018 used (notably
`treatment_outcome_first_course`, aliased in BCR XML as
`primary_therapy_outcome_success`). Those fields *are* preserved in the
GDC's Clinical Supplement files — the original BCR forms shipped in BCR
XML or BCR biotab (TSV) format.

We use **biotab** as the source: per-project per-form TSVs, ~6-7 files
per project (~150 files total for TCGA), versus ~11,167 per-patient
BCR XMLs. Coverage check (2026-05): all 33 TCGA projects have a
`clinical_patient_<proj>.txt`; 32/33 have `clinical_follow_up_v*` and
`clinical_nte` (LAML missing both — consistent with Liu's exclusion of
LAML from DSS / PFI / DFI).

GDC portal note: the larger `data_format=tsv` Clinical Supplement bucket
(~19,566 cases) is entirely non-TCGA — 18,004 cases from Foundation
Medicine (FM) and 1,562 from MP2PRT, neither of which ships in `bcr
biotab`. For TCGA the format choice is bcr biotab vs bcr xml; we picked
biotab because per-project tabular is 50× fewer files than per-patient
XML and carries the same Liu-relevant fields (verified by cross-check on
TCGA-DK-A2I6's `treatment_outcome_first_course`).

For DFI's `treatment_outcome_first_course` specifically, the patient
form may carry "[Not Available]" while the follow-up form carries
"Complete Remission/Response" for the same patient — Liu read both, and
so do we (`combine_supplements_per_case` keeps the latest non-NA value).

Layout on disk:

    <data-dir>/raw/<project_id>/clinical_supplement/
        nationwidechildrens.org_clinical_patient_<proj>.txt
        nationwidechildrens.org_clinical_follow_up_v*_<proj>.txt
        nationwidechildrens.org_clinical_nte_<proj>.txt
        nationwidechildrens.org_clinical_drug_<proj>.txt
        nationwidechildrens.org_clinical_radiation_<proj>.txt
        manifest.json
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Any

from tcga2hf_pipeline.gdc import GDCClient, and_, eq

# Ordered list (not a set) so `_form_kind` is deterministic. Order encodes
# priority for filenames that match more than one form: e.g. BLCA's
# `nationwidechildrens.org_clinical_follow_up_v4.0_nte_blca.txt` carries
# both follow-up structure (one row per encounter) and NTE fields populated
# when relevant — putting `clinical_follow_up` first means we route it to
# the follow_up table (where its row-per-encounter shape is the natural fit).
# `clinical_omf` is "Other Mutation Files" (germline mutation records);
# we ship it for completeness with the other Clinical Supplement biotabs.
WANTED_FORMS: list[str] = [
    "clinical_follow_up",
    "clinical_patient",
    "clinical_nte",
    "clinical_drug",
    "clinical_radiation",
    "clinical_ablation",
    "clinical_omf",
]

# `[Not Available]` / `[Unknown]` etc. all encode "no data" in BCR biotab.
_NA_TOKENS: set[str] = {
    "[Not Available]",
    "[Not Applicable]",
    "[Unknown]",
    "[Discrepancy]",
    "[Not Evaluated]",
    "[Completed]",
    "[Pending]",
    "",
}

FILE_FIELDS: list[str] = [
    "file_id",
    "file_name",
    "file_size",
    "md5sum",
    "access",
    "data_format",
    # Carried so the `files` table can describe these rows like every other
    # modality's. Without them a biotab shows up with a null `data_type`,
    # which reads as missing data rather than as a project-level form.
    "data_type",
    "data_category",
]


def _form_kind(file_name: str) -> str | None:
    """Classify a biotab filename by form prefix; None if it's not one we want."""
    for kind in WANTED_FORMS:
        # Match e.g. "nationwidechildrens.org_clinical_patient_blca.txt" by checking
        # for the form name preceded by `_` (so `clinical_drug` doesn't match
        # `clinical_drug_addendum` accidentally; the GDC doesn't have that today
        # but stays robust as new forms are added).
        if f"_{kind}_" in file_name or f"_{kind}." in file_name:
            return kind
    return None


def list_clinical_supplement_files(
    client: GDCClient,
    project_id: str,
) -> list[dict[str, Any]]:
    """Return all open-access biotab Clinical Supplement files for a project."""
    clauses = [
        eq("cases.project.project_id", project_id),
        eq("access", "open"),
        eq("data_type", "Clinical Supplement"),
        eq("data_format", "bcr biotab"),
    ]
    return client.files(filters=and_(*clauses), fields=FILE_FIELDS, page_size=200)


def fetch_clinical_supplements(
    client: GDCClient,
    project_id: str,
    out_dir: Path,
    skip_existing: bool = True,
    batch_size: int = 50,
) -> list[dict[str, Any]]:
    """Download relevant biotab files to <out_dir>/<file_name> and write manifest.json.

    Filters to the forms in `WANTED_FORMS`; skips biospecimen / ssf files.
    Existing files (matching size) are skipped unless skip_existing=False.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    hits = list_clinical_supplement_files(client, project_id)

    eligible: list[dict[str, Any]] = []
    skipped_other: list[dict[str, Any]] = []
    for hit in hits:
        kind = _form_kind(hit["file_name"])
        if kind is None:
            skipped_other.append(hit)
            continue
        hit["_form_kind"] = kind
        eligible.append(hit)

    cached_ids: set[str] = set()
    to_download_ids: list[str] = []
    for hit in eligible:
        target = out_dir / hit["file_name"]
        if skip_existing and target.exists() and target.stat().st_size == hit.get("file_size"):
            cached_ids.add(hit["file_id"])
        else:
            to_download_ids.append(hit["file_id"])

    if len(to_download_ids) >= 2:
        client.bulk_download(to_download_ids, out_dir, batch_size=batch_size)
    elif len(to_download_ids) == 1:
        only = next(h for h in eligible if h["file_id"] == to_download_ids[0])
        client.download(only["file_id"], out_dir / only["file_name"])

    # Per-file version provenance, same as the genomic modalities — see
    # `GDCClient.versions`. Biotabs are revised across releases more often
    # than most GDC files, so `gdc_superseded` is worth carrying here.
    versions = client.versions([h["file_id"] for h in eligible]) if eligible else {}

    manifest: list[dict[str, Any]] = []
    for hit in eligible:
        status = "cached" if hit["file_id"] in cached_ids else "downloaded"
        version = versions.get(hit["file_id"], {})
        manifest.append({
            "file_id": hit["file_id"],
            "file_name": hit["file_name"],
            "file_size": hit.get("file_size"),
            "md5sum": hit.get("md5sum"),
            "data_format": hit.get("data_format"),
            "data_type": hit.get("data_type"),
            "data_category": hit.get("data_category"),
            "access": hit.get("access"),
            # Deliberately no `cases`: a biotab is a *project*-level form
            # covering every case, so exploding it per case would multiply
            # one file into 51 rows that each imply a per-patient file.
            "form_kind": hit["_form_kind"],
            "gdc_version": version.get("version"),
            "gdc_first_release": version.get("release"),
            "gdc_superseded": bool(
                version and version.get("latest_id") not in (None, hit["file_id"])
            ),
            "_status": status,
        })
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    return manifest


def md5_of(path: Path, chunk: int = 1 << 16) -> str:
    h = hashlib.md5()
    with path.open("rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------


def parse_biotab(path: Path) -> list[dict[str, str]]:
    """Parse a BCR biotab TSV into a list of row dicts keyed by column name.

    BCR biotabs have a 3-row header:
      row 0: column names (the actual schema)
      row 1: aliased preferred names
      row 2: CDE_ID:<id> tags
    Data rows start at row 3.
    """
    # BCR biotabs use Windows-1252 / Latin-1, not UTF-8 (encountered \x96 = en-dash
    # in the wild). latin-1 round-trips every byte so we never crash on encoding.
    with path.open(newline="", encoding="latin-1") as f:
        reader = csv.reader(f, delimiter="\t")
        rows = list(reader)
    if len(rows) < 4:
        return []
    header = rows[0]
    out: list[dict[str, str]] = []
    for raw in rows[3:]:
        if not raw or all(c == "" for c in raw):
            continue
        rec = {header[i]: raw[i] if i < len(raw) else "" for i in range(len(header))}
        out.append(rec)
    return out


def _is_populated(v: str | None) -> bool:
    return v is not None and v not in _NA_TOKENS


def load_supplements_for_project(supp_dir: Path) -> dict[str, dict[str, Any]]:
    """Parse every biotab in supp_dir and group records by bcr_patient_barcode.

    Returns a dict:
        {case_submitter_id: {
            "patient": dict | None,
            "follow_ups": list[dict],
            "ntes": list[dict],
            "drugs": list[dict],
            "radiations": list[dict],
            "ablations": list[dict],
            "omfs": list[dict],
        }}

    The `patient` slot is a single dict (one initial form per patient);
    every other slot is a list (a patient may have N follow-ups, N drugs, etc.).
    """
    if not supp_dir.exists():
        return {}

    by_case: dict[str, dict[str, Any]] = {}

    def _slot(case_id: str) -> dict[str, Any]:
        return by_case.setdefault(case_id, {
            "patient": None,
            "follow_ups": [],
            "ntes": [],
            "drugs": [],
            "radiations": [],
            "ablations": [],
            "omfs": [],
        })

    for path in sorted(supp_dir.glob("*.txt")):
        kind = _form_kind(path.name)
        if kind is None:
            continue
        records = parse_biotab(path)
        for rec in records:
            case_id = rec.get("bcr_patient_barcode")
            if not case_id or not case_id.startswith("TCGA-"):
                continue
            slot = _slot(case_id)
            if kind == "clinical_patient":
                slot["patient"] = rec
            elif kind == "clinical_follow_up":
                slot["follow_ups"].append(rec)
            elif kind == "clinical_nte":
                slot["ntes"].append(rec)
            elif kind == "clinical_drug":
                slot["drugs"].append(rec)
            elif kind == "clinical_radiation":
                slot["radiations"].append(rec)
            elif kind == "clinical_ablation":
                slot["ablations"].append(rec)
            elif kind == "clinical_omf":
                slot["omfs"].append(rec)
    return by_case


def iter_values(supp: dict[str, Any], field: str) -> list[str]:
    """All non-NA values of `field` across patient + follow-up forms, in form order.

    Liu's "first course" outcome can land on the patient form OR any later
    follow-up form (e.g. TCGA-DK-A2I6's value "Complete Remission/Response"
    only appears on a follow-up; the patient form is "[Not Available]").
    A patient may also have *multiple* follow-up forms with different
    values — TCGA-2G-AAGA has both "Complete Remission/Response" and
    "No Measureable Tumor or Tumor Markers" on different follow-ups.

    We return every populated value in encounter order so callers decide
    the right collapse (Liu's disease-free check wants "did any form ever
    record CR/CRR?" — see `any_disease_free_signal`).
    """
    out: list[str] = []
    patient = supp.get("patient") or {}
    candidate = patient.get(field)
    if _is_populated(candidate):
        out.append(candidate)
    follow_ups = supp.get("follow_ups") or []
    for fu in sorted(follow_ups, key=lambda r: r.get("bcr_followup_barcode") or ""):
        candidate = fu.get(field)
        if _is_populated(candidate):
            out.append(candidate)
    return out


def any_disease_free_signal(supp: dict[str, Any]) -> bool:
    """True if any patient/follow-up form recorded a Liu disease-free outcome.

    Liu's algorithm: `treatment_outcome_first_course == "Complete Remission/Response"`
    counts as disease-free at end of first course. Modern data also uses
    `"Complete Response"` (no slash) for the same concept.
    """
    cr_set = {"Complete Response", "Complete Remission/Response"}
    return any(v in cr_set for v in iter_values(supp, "treatment_outcome_first_course"))


def first_value(supp: dict[str, Any], field: str) -> str | None:
    """First populated value (patient form, then earliest follow-up). Convenience helper."""
    vals = iter_values(supp, field)
    return vals[0] if vals else None


# ---------------------------------------------------------------------------
# Pipeline integration
# ---------------------------------------------------------------------------


def attach_supplements(
    rows: list[dict[str, Any]],
    supplements_by_case: dict[str, dict[str, Any]],
) -> list[dict[str, Any]]:
    """Attach raw supplement records to each patient row under `clinical_supplement`.

    Survival re-derivation can then read `row["clinical_supplement"]` to
    pull `treatment_outcome_first_course` (and other fields) from the
    biotab data when the harmonized API value is missing.
    """
    for row in rows:
        case_id = row.get("case_submitter_id")
        row["clinical_supplement"] = supplements_by_case.get(case_id)
    return rows


# Map our tabular table-name suffix -> form_kind in `WANTED_FORMS`. The
# tabular dataset will emit one parquet per (project, kind) pair under
# `<processed>/<project>/clinical_supplement_<suffix>/data.parquet`. Schema
# is inferred per project (no cross-project union) — different cancer types
# carry different BCR fields and we ship each project's set faithfully.
TABULAR_FORM_KINDS: list[str] = [
    "patient",
    "follow_up",
    "nte",
    "drug",
    "radiation",
    "ablation",
    "omf",
]


def _kind_to_form(suffix: str) -> str:
    return f"clinical_{suffix}"


def build_tabular_rows(
    supp_dir: Path,
    case_id_to_submitter: dict[str, str] | None = None,
) -> dict[str, list[dict[str, Any]]]:
    """Per-form list of biotab records for tabular emission.

    Returns `{tabular_table_suffix: [row, ...]}` where each row is the raw
    biotab record (one per data row in the source file) prefixed with
    `case_submitter_id`. Rows whose `bcr_patient_barcode` doesn't start
    with `TCGA-` are skipped (occasional placeholder/header artifacts).
    """
    out: dict[str, list[dict[str, Any]]] = {suffix: [] for suffix in TABULAR_FORM_KINDS}
    if not supp_dir.exists():
        return out

    for path in sorted(supp_dir.glob("*.txt")):
        kind = _form_kind(path.name)
        if kind is None:
            continue
        suffix = kind.removeprefix("clinical_")
        if suffix not in out:
            continue
        for rec in parse_biotab(path):
            barcode = rec.get("bcr_patient_barcode")
            if not barcode or not barcode.startswith("TCGA-"):
                continue
            # `case_submitter_id` == bcr_patient_barcode in the biotab; we
            # surface it as our canonical FK column for HF joinability.
            row = {"case_submitter_id": barcode, **rec}
            out[suffix].append(row)
    return out
