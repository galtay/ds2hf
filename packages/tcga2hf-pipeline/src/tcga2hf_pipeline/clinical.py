from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq
from tcga2hf.schema import (
    ALIQUOT_FIELDS,
    ANALYTE_FIELDS,
    ANNOTATION_FIELDS,
    CENTER_FIELDS,
    DEMOGRAPHIC_FIELDS,
    DIAGNOSIS_FIELDS,
    EXPOSURE_FIELDS,
    FAMILY_HISTORY_FIELDS,
    FOLLOW_UP_FIELDS,
    MOLECULAR_TEST_FIELDS,
    OTHER_CLINICAL_ATTRIBUTE_FIELDS,
    PATHOLOGY_DETAIL_FIELDS,
    PATIENTS,
    PORTION_FIELDS,
    PROGRAM_FIELDS,
    SAMPLE_FIELDS,
    SLIDE_FIELDS,
    SSGSEA_COLLECTIONS,
    TISSUE_SOURCE_SITE_FIELDS,
    TREATMENT_FIELDS,
)

from tcga2hf_pipeline.gdc import GDCClient, in_

# GDC `/cases` exposes 41 expandable groups. We request every one of them
# except the `files.*` subtree — 29 in total.
#
# Two things force this into two requests rather than one:
#
#   1. **The API silently truncates a long `expand` list.** Asking for all 29
#      at once returns HTTP 200 with `hits: []` and `total: null` — no error,
#      no warning. The ceiling measured against the live API is 21 groups.
#   2. GDC `expand` is per-level: every intermediate must be listed for its
#      own fields to come back, so the deep biospecimen chain burns entries
#      quickly (without `samples.portions.analytes`, `analyte_type` is
#      missing).
#
# So the groups are split along the natural seam — clinical vs biospecimen —
# and the two responses are merged on `case_id` in `fetch_clinical`.
#
# `files.*` is excluded deliberately rather than for lack of room: it costs
# 190 KB per case against 16 KB for everything here (~2.2 GB pan-TCGA), and
# it duplicates the `files` table, which is built from the `/files`
# endpoint's own manifests and carries better provenance (md5, gdc_version,
# supersession) than the nested copy does.
CLINICAL_EXPANSIONS = [
    "annotations",
    "demographic",
    "diagnoses",
    "diagnoses.annotations",
    "diagnoses.pathology_details",
    "diagnoses.treatments",
    "exposures",
    "family_histories",
    "follow_ups",
    "follow_ups.molecular_tests",
    "follow_ups.other_clinical_attributes",
    "project",
    "project.program",
    "summary",
    "summary.data_categories",
    "summary.experimental_strategies",
    "tissue_source_site",
]

BIOSPECIMEN_EXPANSIONS = [
    "samples",
    "samples.annotations",
    "samples.portions",
    "samples.portions.analytes",
    "samples.portions.analytes.aliquots",
    "samples.portions.analytes.aliquots.annotations",
    "samples.portions.analytes.aliquots.center",
    "samples.portions.analytes.annotations",
    "samples.portions.annotations",
    "samples.portions.center",
    "samples.portions.slides",
    "samples.portions.slides.annotations",
]

EXPANSIONS = CLINICAL_EXPANSIONS + BIOSPECIMEN_EXPANSIONS

# `fields` is required by the API but does not restrict expanded entities.
# Case-level scalar fields (per GDC dictionary): index_date names the timeline
# anchor (TCGA: "Diagnosis"); the days_to_* and consent/follow-up status all
# share that anchor. Have to request each explicitly — they don't come back
# with `expand`.
TOP_LEVEL_FIELDS = [
    "case_id",
    "submitter_id",
    "primary_site",
    "disease_type",
    "index_date",
    "consent_type",
    "days_to_consent",
    "days_to_lost_to_followup",
    "lost_to_followup",
]


def fetch_clinical(projects: list[str], client: GDCClient) -> list[dict[str, Any]]:
    """Fetch every requestable field of each case, merging two `/cases` calls.

    See `CLINICAL_EXPANSIONS` for why this is two requests. Each response is
    a complete case record for its half of the tree; they are merged on
    `case_id`, which both halves always carry.

    Raises if either half comes back empty while the other did not — that is
    the signature of the API's silent `expand`-too-long truncation, and
    treating it as "this project has no cases" would quietly ship a dataset
    missing its entire biospecimen tree.
    """
    filters = in_("project.project_id", projects)

    def _call(expansions: list[str]) -> list[dict[str, Any]]:
        return client.cases(
            filters=filters,
            fields=TOP_LEVEL_FIELDS,
            expand=expansions,
            page_size=200,
        )

    clinical_hits = _call(CLINICAL_EXPANSIONS)
    biospecimen_hits = _call(BIOSPECIMEN_EXPANSIONS)
    if bool(clinical_hits) != bool(biospecimen_hits):
        raise RuntimeError(
            f"GDC returned {len(clinical_hits)} clinical and "
            f"{len(biospecimen_hits)} biospecimen cases for {projects}. An "
            "empty half means the `expand` list was silently truncated; "
            "shorten CLINICAL_EXPANSIONS / BIOSPECIMEN_EXPANSIONS."
        )

    merged: dict[str, dict[str, Any]] = {}
    for hit in clinical_hits:
        merged[hit["case_id"]] = dict(hit)
    for hit in biospecimen_hits:
        merged.setdefault(hit["case_id"], {}).update(hit)
    return list(merged.values())


def _pick(obj: dict[str, Any] | None, fields: list[pa.Field]) -> dict[str, Any] | None:
    """Return a dict containing exactly the named fields (None for missing).

    Returns None if the source object itself is missing — pyarrow renders that as a
    null struct, distinct from a struct of all-null fields.
    """
    if obj is None:
        return None
    return {f.name: obj.get(f.name) for f in fields}


def _sort_temporal(
    items: list[dict[str, Any]], days_field: str, id_field: str
) -> list[dict[str, Any]]:
    """Sort by (days_field asc, id_field asc) with null `days_field` last.

    The `is None` sentinel keeps None values from being compared to ints — Python
    short-circuits tuple comparison on the first differing element, so a True/False
    bool prefix safely segregates the nullable group.
    """
    return sorted(
        items,
        key=lambda x: (x.get(days_field) is None, x.get(days_field), x.get(id_field) or ""),
    )


def _sort_by_id(items: list[dict[str, Any]], id_field: str) -> list[dict[str, Any]]:
    return sorted(items, key=lambda x: x.get(id_field) or "")


def _scalars(obj: dict[str, Any], fields: list[pa.Field], nested: set[str]) -> dict[str, Any]:
    """Every field of `fields` except the nested containers listed in `nested`."""
    return {f.name: obj.get(f.name) for f in fields if f.name not in nested}


def _annotations(obj: dict[str, Any]) -> list[dict[str, Any]]:
    """Curator annotations, sorted by id so the order is reproducible."""
    picked = [_pick(a, ANNOTATION_FIELDS) for a in (obj.get("annotations") or [])]
    return _sort_by_id(picked, "annotation_id")


def _diagnosis_dict(dx: dict[str, Any]) -> dict[str, Any]:
    out = _scalars(dx, DIAGNOSIS_FIELDS, {"treatments", "annotations", "pathology_details"})
    treatments = [_pick(tx, TREATMENT_FIELDS) for tx in (dx.get("treatments") or [])]
    out["treatments"] = _sort_temporal(treatments, "days_to_treatment_start", "treatment_id")
    out["annotations"] = _annotations(dx)
    details = [_pick(d, PATHOLOGY_DETAIL_FIELDS) for d in (dx.get("pathology_details") or [])]
    out["pathology_details"] = _sort_by_id(details, "pathology_detail_id")
    return out


def _follow_up_dict(fu: dict[str, Any]) -> dict[str, Any]:
    out = _scalars(fu, FOLLOW_UP_FIELDS, {"molecular_tests", "other_clinical_attributes"})
    tests = [_pick(t, MOLECULAR_TEST_FIELDS) for t in (fu.get("molecular_tests") or [])]
    out["molecular_tests"] = _sort_by_id(tests, "molecular_test_id")
    others = [
        _pick(o, OTHER_CLINICAL_ATTRIBUTE_FIELDS)
        for o in (fu.get("other_clinical_attributes") or [])
    ]
    out["other_clinical_attributes"] = _sort_by_id(others, "other_clinical_attribute_id")
    return out


def _aliquot_dict(aliquot: dict[str, Any]) -> dict[str, Any]:
    out = _scalars(aliquot, ALIQUOT_FIELDS, {"annotations", "center"})
    out["annotations"] = _annotations(aliquot)
    out["center"] = _pick(aliquot.get("center"), CENTER_FIELDS)
    return out


def _analyte_dict(analyte: dict[str, Any]) -> dict[str, Any]:
    out = _scalars(analyte, ANALYTE_FIELDS, {"aliquots", "annotations"})
    aliquots = [_aliquot_dict(a) for a in (analyte.get("aliquots") or [])]
    out["aliquots"] = _sort_by_id(aliquots, "aliquot_id")
    out["annotations"] = _annotations(analyte)
    return out


def _portion_dict(portion: dict[str, Any]) -> dict[str, Any]:
    out = _scalars(portion, PORTION_FIELDS, {"analytes", "annotations", "center", "slides"})
    analytes = [_analyte_dict(a) for a in (portion.get("analytes") or [])]
    out["analytes"] = _sort_by_id(analytes, "analyte_id")
    out["annotations"] = _annotations(portion)
    out["center"] = _pick(portion.get("center"), CENTER_FIELDS)
    slides = [_pick(sl, SLIDE_FIELDS) for sl in (portion.get("slides") or [])]
    out["slides"] = _sort_by_id(slides, "slide_id")
    return out


def _sample_dict(sample: dict[str, Any]) -> dict[str, Any]:
    out = _scalars(sample, SAMPLE_FIELDS, {"portions", "annotations"})
    portions = [_portion_dict(p) for p in (sample.get("portions") or [])]
    out["portions"] = _sort_by_id(portions, "portion_id")
    out["annotations"] = _annotations(sample)
    return out


def _summary_dict(case: dict[str, Any]) -> dict[str, Any] | None:
    """GDC's own file tallies for the case, counts sorted for reproducibility."""
    summary = case.get("summary")
    if summary is None:
        return None
    return {
        "file_count": summary.get("file_count"),
        "file_size": summary.get("file_size"),
        "data_categories": sorted(
            (
                {
                    "data_category": d.get("data_category"),
                    "file_count": d.get("file_count"),
                }
                for d in (summary.get("data_categories") or [])
            ),
            key=lambda d: d["data_category"] or "",
        ),
        "experimental_strategies": sorted(
            (
                {
                    "experimental_strategy": e.get("experimental_strategy"),
                    "file_count": e.get("file_count"),
                }
                for e in (summary.get("experimental_strategies") or [])
            ),
            key=lambda e: e["experimental_strategy"] or "",
        ),
    }


def _patient_row(case: dict[str, Any]) -> dict[str, Any]:
    project = case.get("project") or {}
    diagnoses = [_diagnosis_dict(dx) for dx in (case.get("diagnoses") or [])]
    follow_ups = [_follow_up_dict(fu) for fu in (case.get("follow_ups") or [])]
    exposures = [_pick(ex, EXPOSURE_FIELDS) for ex in (case.get("exposures") or [])]
    family_histories = [
        _pick(fh, FAMILY_HISTORY_FIELDS) for fh in (case.get("family_histories") or [])
    ]
    samples = [_sample_dict(s) for s in (case.get("samples") or [])]
    case_id = case.get("case_id")
    return {
        "case_id": case_id,
        "case_submitter_id": case.get("submitter_id"),
        "project_id": project.get("project_id"),
        "gdc_portal_url": (
            f"https://portal.gdc.cancer.gov/cases/{case_id}" if case_id else None
        ),
        "primary_site": case.get("primary_site"),
        "disease_type": case.get("disease_type"),
        "index_date": case.get("index_date"),
        "consent_type": case.get("consent_type"),
        "days_to_consent": case.get("days_to_consent"),
        "days_to_lost_to_followup": case.get("days_to_lost_to_followup"),
        "lost_to_followup": case.get("lost_to_followup"),
        "annotations": _annotations(case),
        "tissue_source_site": _pick(
            case.get("tissue_source_site"), TISSUE_SOURCE_SITE_FIELDS
        ),
        "program": _pick(project.get("program"), PROGRAM_FIELDS),
        "summary": _summary_dict(case),
        "demographic": _pick(case.get("demographic"), DEMOGRAPHIC_FIELDS),
        "diagnoses": _sort_temporal(diagnoses, "days_to_diagnosis", "diagnosis_id"),
        "follow_ups": _sort_temporal(follow_ups, "days_to_follow_up", "follow_up_id"),
        "exposures": _sort_by_id(exposures, "exposure_id"),
        "family_histories": _sort_by_id(family_histories, "family_history_id"),
        "samples": _sort_temporal(samples, "days_to_collection", "sample_id"),
        # Per-file modality columns default to []; build step fills them in
        # if the corresponding raw data has been fetched. Every one of these
        # must be listed: a project can legitimately have zero files for a
        # modality (TCGA-LAML has no RPPA and no pathology reports), in which
        # case `attach` is never called and this default is the only thing
        # that puts the column on the row.
        "samples_masked_somatic_mutation": [],
        "samples_gene_expression_quantification": [],
        "samples_pathology_report": [],
        "samples_allele_specific_copy_number_segment": [],
        "samples_masked_copy_number_segment": [],
        "samples_mirna_expression_quantification": [],
        "samples_protein_expression_quantification": [],
        **{f"samples_ssgsea_{_c}": [] for _c in SSGSEA_COLLECTIONS},
        # Re-derived survival endpoints (Liu et al. 2018 algorithm). Build
        # step fills the struct via `survival.attach_survival`. Initialized
        # here as an all-null struct so pyarrow's `from_pylist(rows, schema=PATIENTS)`
        # can lift the column even when survival hasn't been attached yet.
        "survival_derived": {
            "os_event": None, "os_time": None,
            "dss_event": None, "dss_time": None,
            "pfi_event": None, "pfi_time": None,
            "dfi_event": None, "dfi_time": None,
        },
    }


def to_patient_rows(cases: Iterable[dict[str, Any]]) -> list[dict[str, Any]]:
    return [_patient_row(case) for case in cases]


# Row group size targeting the HuggingFace Dataset Viewer's 300 MB scan limit
# (see https://huggingface.co/docs/hub/datasets-data-files-configuration). Each
# patient row is ~1-1.5 MB compressed, dominated by the 60k-gene expression
# arrays, so 50 rows per group ≈ 50-75 MB — comfortably below the limit and
# small enough that the Viewer's per-page random access is snappy.
_ROW_GROUP_SIZE = 50


_SUPPLEMENT_FORMS = ("patient", "follow_ups", "ntes", "drugs",
                     "radiations", "ablations", "omfs")


def _build_clinical_supplement_column(rows: list[dict[str, Any]]) -> pa.Array | None:
    """Build a per-project pyarrow column for the flex clinical_supplement struct.

    BCR biotab schemas vary by cancer type, so we don't enumerate sub-fields
    in the global PATIENTS schema; instead we let pyarrow infer the struct
    type from the actual data per project. Form keys with no data for any
    patient in the project are dropped (e.g. LAML has no follow_ups, so its
    `clinical_supplement` struct has no `follow_ups` key) — this keeps
    pyarrow from trying to infer a type from all-empty lists.

    Returns None if no patient in the project has any supplement data.
    """
    payloads: list[dict[str, Any] | None] = []
    keys_with_data: set[str] = set()
    for row in rows:
        supp = row.get("clinical_supplement")
        if not supp:
            payloads.append(None)
            continue
        keep: dict[str, Any] = {}
        for k in _SUPPLEMENT_FORMS:
            v = supp.get(k)
            if k == "patient":
                if v:
                    keep[k] = v
                    keys_with_data.add(k)
            else:
                # list-valued forms — keep the list (even if empty for this
                # patient) only if at least one patient in the project has
                # entries for it. Otherwise drop the key entirely.
                if v:
                    keys_with_data.add(k)
                keep[k] = v or []
        payloads.append(keep)

    if not keys_with_data:
        return None

    # Drop any form keys that turned out to be empty for every patient in
    # the project (pyarrow can't infer a struct type from all-empty lists).
    cleaned: list[dict[str, Any] | None] = []
    for p in payloads:
        if p is None:
            cleaned.append(None)
            continue
        cleaned.append({k: v for k, v in p.items() if k in keys_with_data})
    return pa.array(cleaned)


def _build_biospecimen_supplement_column(rows: list[dict[str, Any]]) -> pa.Array | None:
    """Build a per-project pyarrow column for the flex biospecimen_supplement struct.

    Same flex-schema reasoning as `_build_clinical_supplement_column`: BCR
    biotab columns vary by cancer type and by submitting centre, so the
    struct type is inferred per project rather than enumerated in the global
    PATIENTS schema.

    Simpler than its clinical counterpart in one way — every biospecimen form
    is list-valued (a patient has N samples, N aliquots, N slides), so there
    is no single-dict `patient` slot to special-case. Form keys empty for
    every patient in the project are dropped, since pyarrow cannot infer a
    struct type from all-empty lists (only TCGA-LUAD has `cqcf`, and only 9
    projects have `auxiliary`).

    Returns None if no patient in the project has any biospecimen data.
    """
    from tcga2hf_pipeline.biospecimen_supplement import TABULAR_FORM_KINDS

    payloads: list[dict[str, Any] | None] = []
    keys_with_data: set[str] = set()
    for row in rows:
        supp = row.get("biospecimen_supplement")
        if not supp:
            payloads.append(None)
            continue
        keep: dict[str, Any] = {}
        for k in TABULAR_FORM_KINDS:
            v = supp.get(k)
            if v:
                keys_with_data.add(k)
            keep[k] = v or []
        payloads.append(keep)

    if not keys_with_data:
        return None

    cleaned: list[dict[str, Any] | None] = []
    for p in payloads:
        if p is None:
            cleaned.append(None)
            continue
        cleaned.append({k: v for k, v in p.items() if k in keys_with_data})
    return pa.array(cleaned)


def write_patients(rows: list[dict[str, Any]], processed_dir: Path, project_id: str) -> Path:
    """Write a single project's patient rows to <project_id>/data.parquet.

    The per-project directory layout matches the HuggingFace `configs:` convention
    so each TCGA project surfaces as its own loadable subset.

    Row groups are sized for the HF Dataset Viewer (see `_ROW_GROUP_SIZE`) and
    a page index is written so the Viewer can read only the bytes it needs to
    render a page of rows rather than scanning the whole row group.

    The optional `clinical_supplement` column is built separately from the
    fixed PATIENTS schema because BCR biotab fields vary by cancer type;
    its struct shape is inferred per project from the row data.
    """
    out_path = processed_dir / project_id / "data.parquet"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    # Strip clinical_supplement before applying the strict PATIENTS schema
    # (pyarrow's from_pylist silently drops extra dict keys, but being
    # explicit is clearer); attach it as a separate flex-typed column after.
    supp_col = _build_clinical_supplement_column(rows)
    bio_col = _build_biospecimen_supplement_column(rows)
    table = pa.Table.from_pylist(rows, schema=PATIENTS)
    if supp_col is not None:
        table = table.append_column("clinical_supplement", supp_col)
    if bio_col is not None:
        table = table.append_column("biospecimen_supplement", bio_col)
    pq.write_table(
        table,
        out_path,
        row_group_size=_ROW_GROUP_SIZE,
        write_page_index=True,
    )
    return out_path
