"""Biospecimen Supplement biotab fetcher + parser.

The sibling of `clinical_supplement`. Where the clinical biotabs describe
the *patient* (diagnosis, treatment, follow-up), the biospecimen biotabs
describe the *specimen chain* — how a tumour got from the operating room to
a sequencer — and the pathologist's read on each slide along the way.

Some of this duplicates what the harmonized `/cases` API already nests into
our `cases` table (sample / portion / analyte / aliquot ids and types). Much
of it does not. The fields with no `/cases` equivalent are the reason to
ship these:

  - `biospecimen_slide` — per-slide `percent_tumor_nuclei`,
    `percent_necrosis`, `percent_stromal_cells`, `section_location`. This
    is the QC layer every "is this sample actually tumour?" question
    bottoms out in, and it is the standard covariate for deconvolution and
    purity work.
  - `biospecimen_analyte` — `a260_a280_ratio`, `concentration`, `spectro-
    photometer_method`. Nucleic-acid quality, which drives batch effects.
  - `biospecimen_protocol` / `biospecimen_shipment_portion` — the plate,
    shipment and centre a specimen moved through, i.e. the raw material
    for batch-effect analysis.
  - `ssf_tumor_samples` / `ssf_normal_controls` — site-specific factors,
    the disease-specific pathology fields (e.g. Gleason components,
    Breslow depth) that the pan-cancer clinical schema has no column for.
  - `biospecimen_cqcf` — the submitting centre's clinical quality control
    form, present only for the `genome.wustl.edu`-submitted LUAD set.

Like the clinical supplements these are **flex-schema**: the column set
differs by project and by submitting centre, so each (project, form) pair
gets its own parquet with its own inferred schema rather than a padded
pan-cancer union.

Two submitters ship TCGA biospecimen biotabs — `nationwidechildrens.org`
for 334 of the 340 files and `genome.wustl.edu` for 6 (all TCGA-LUAD). The
form classifier therefore matches on the form name embedded in the
filename and never on the submitter prefix.

Layout on disk:

    <data-dir>/raw/<project_id>/biospecimen_supplement/
        nationwidechildrens.org_biospecimen_sample_chol.txt
        nationwidechildrens.org_biospecimen_slide_chol.txt
        ...
        manifest.json
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from tcga2hf_pipeline.clinical_supplement import parse_biotab
from tcga2hf_pipeline.gdc import GDCClient, and_, eq

# Ordered so `_form_kind` is deterministic, and ordered longest-first within
# each family so a specific form wins over a prefix of itself. The GDC has
# no such collision today (`biospecimen_shipment_portion` does not contain
# the substring `_biospecimen_portion_`, and `biospecimen_diagnostic_slides`
# does not contain `_biospecimen_slide_`), but the ordering makes that a
# property of this list rather than a lucky accident of GDC naming.
WANTED_FORMS: list[str] = [
    "biospecimen_diagnostic_slides",
    "biospecimen_shipment_portion",
    "biospecimen_slide",
    "biospecimen_portion",
    "biospecimen_aliquot",
    "biospecimen_analyte",
    "biospecimen_protocol",
    "biospecimen_sample",
    "biospecimen_cqcf",
    "ssf_tumor_samples",
    "ssf_normal_controls",
    "auxiliary",
]

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
    """Classify a biospecimen biotab filename by form name; None if unwanted."""
    for kind in WANTED_FORMS:
        if f"_{kind}_" in file_name or f"_{kind}." in file_name:
            return kind
    return None


def list_biospecimen_supplement_files(
    client: GDCClient,
    project_id: str,
) -> list[dict[str, Any]]:
    """Return all open-access biotab Biospecimen Supplement files for a project."""
    clauses = [
        eq("cases.project.project_id", project_id),
        eq("access", "open"),
        eq("data_type", "Biospecimen Supplement"),
        eq("data_format", "bcr biotab"),
    ]
    return client.files(filters=and_(*clauses), fields=FILE_FIELDS, page_size=200)


def fetch_biospecimen_supplements(
    client: GDCClient,
    project_id: str,
    out_dir: Path,
    skip_existing: bool = True,
    batch_size: int = 50,
) -> list[dict[str, Any]]:
    """Download biotab files to <out_dir>/<file_name> and write manifest.json.

    Structurally identical to `clinical_supplement.fetch_clinical_supplements`
    — the difference is only which `data_type` and which `WANTED_FORMS` are
    in play. Kept as a separate function rather than parameterising the
    clinical one because the two modalities land in separate raw directories
    and carry separate GDC release stamps.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    hits = list_biospecimen_supplement_files(client, project_id)

    eligible: list[dict[str, Any]] = []
    for hit in hits:
        kind = _form_kind(hit["file_name"])
        if kind is None:
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

    versions = client.versions([h["file_id"] for h in eligible]) if eligible else {}

    manifest: list[dict[str, Any]] = []
    for hit in eligible:
        status = "cached" if hit["file_id"] in cached_ids else "downloaded"
        version = versions.get(hit["file_id"], {})
        manifest.append(
            {
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
            }
        )
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    return manifest


# ---------------------------------------------------------------------------
# Tabular emission
# ---------------------------------------------------------------------------

# Suffixes used in the tabular table names (`biospecimen_supplement_<suffix>`).
# Derived from WANTED_FORMS by stripping the redundant `biospecimen_` prefix
# — the table name already says it.
TABULAR_FORM_KINDS: list[str] = [
    "sample",
    "portion",
    "analyte",
    "aliquot",
    "slide",
    "diagnostic_slides",
    "shipment_portion",
    "protocol",
    "cqcf",
    "ssf_tumor_samples",
    "ssf_normal_controls",
    "auxiliary",
]

# Barcode columns that identify a row's patient, in preference order. Only
# `bcr_patient_barcode` is guaranteed present; the specimen-level forms are
# keyed on their own entity and some omit the patient column, in which case
# the patient barcode is the first three hyphen groups of the entity barcode
# (`TCGA-3X-AAV9-01A-...` -> `TCGA-3X-AAV9`) — a stable property of the TCGA
# barcode grammar, not a heuristic.
_BARCODE_COLUMNS: list[str] = [
    "bcr_patient_barcode",
    "bcr_sample_barcode",
    "bcr_portion_barcode",
    "bcr_analyte_barcode",
    "bcr_aliquot_barcode",
    "bcr_slide_barcode",
    "bcr_shipment_portion_barcode",
]

_PATIENT_BARCODE_GROUPS = 3


def _case_submitter_id(rec: dict[str, str]) -> str | None:
    """Recover the patient barcode from whichever barcode column a form carries."""
    for col in _BARCODE_COLUMNS:
        value = (rec.get(col) or "").strip()
        if not value.startswith("TCGA-"):
            continue
        groups = value.split("-")
        if len(groups) < _PATIENT_BARCODE_GROUPS:
            continue
        return "-".join(groups[:_PATIENT_BARCODE_GROUPS])
    return None


def _suffix_for(kind: str) -> str:
    return kind.removeprefix("biospecimen_")


def load_supplements_for_project(supp_dir: Path) -> dict[str, dict[str, list[dict]]]:
    """Parse every biotab in supp_dir and group records by patient barcode.

    Returns `{case_submitter_id: {form_suffix: [record, ...]}}`, with one key
    per entry in `TABULAR_FORM_KINDS`.

    Every slot is a list, unlike `clinical_supplement.load_supplements_for_project`
    where `patient` is a single dict: these forms are all one-row-per-entity
    (a patient has N samples, N portions, N aliquots, N slides), so there is
    no once-per-patient form among them.
    """
    if not supp_dir.exists():
        return {}

    by_case: dict[str, dict[str, list[dict]]] = {}
    for path in sorted(supp_dir.glob("*.txt")):
        kind = _form_kind(path.name)
        if kind is None:
            continue
        suffix = _suffix_for(kind)
        if suffix not in TABULAR_FORM_KINDS:
            continue
        for rec in parse_biotab(path):
            case_submitter_id = _case_submitter_id(rec)
            if case_submitter_id is None:
                continue
            slot = by_case.setdefault(
                case_submitter_id, {k: [] for k in TABULAR_FORM_KINDS}
            )
            slot[suffix].append(rec)
    return by_case


def attach_supplements(
    rows: list[dict[str, Any]],
    supplements_by_case: dict[str, dict[str, list[dict]]],
) -> list[dict[str, Any]]:
    """Attach raw biospecimen records to each patient row under `biospecimen_supplement`.

    Keyed on `case_submitter_id` (the BCR patient barcode) rather than
    `case_id`, since the biotabs carry barcodes and never UUIDs for the case.
    """
    for row in rows:
        row["biospecimen_supplement"] = supplements_by_case.get(row.get("case_submitter_id"))
    return rows


def build_tabular_rows(supp_dir: Path) -> dict[str, list[dict[str, Any]]]:
    """Per-form list of biotab records for tabular emission.

    Returns `{tabular_table_suffix: [row, ...]}`, each row the raw biotab
    record prefixed with `case_submitter_id`. Rows whose patient barcode
    can't be recovered from any barcode column are skipped — that catches
    the placeholder rows BCR leaves in some forms without silently
    attributing them to a patient.

    When two submitters ship the same form for one project (TCGA-LUAD has
    both `nationwidechildrens.org` and `genome.wustl.edu` versions), their
    rows are concatenated into one table. The union of their columns is what
    the parquet schema infers, so a row from one submitter reports null for
    columns only the other defines.
    """
    out: dict[str, list[dict[str, Any]]] = {suffix: [] for suffix in TABULAR_FORM_KINDS}
    if not supp_dir.exists():
        return out

    for path in sorted(supp_dir.glob("*.txt")):
        kind = _form_kind(path.name)
        if kind is None:
            continue
        suffix = _suffix_for(kind)
        if suffix not in out:
            continue
        for rec in parse_biotab(path):
            case_submitter_id = _case_submitter_id(rec)
            if case_submitter_id is None:
                continue
            out[suffix].append({"case_submitter_id": case_submitter_id, **rec})
    return out
