"""Generic open-access GDC file fetcher for genomic / molecular data.

Mutations and expression (and miRNA, copy number, etc. later) all share the same
shape: query /files with a project + data_type filter, download each hit's bytes,
and write a manifest mapping file_id -> case_submitter_id + sample info. The
patient-row build step can join on those FKs later — we keep raw downloads and
manifest separate so iteration on the schema doesn't trigger re-downloads.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from tcga2hf_pipeline.gdc import GDCClient, and_, eq

# Fields we request on each /files hit. The GDC `/files` endpoint
# *replaces* the default response with whatever you list here (it does
# NOT append to a default), so this enumerates every field our pipeline
# consumes — file metadata we save into the manifest, plus the case and
# aliquot linkages we need for downstream join logic. Trimmed to only
# what's actually read; cf. `expression._file_aliquot_and_case` for the
# reason `cases.samples.portions.analytes.aliquots.aliquot_id` is in
# here despite looking deep.
#
# `analysis.workflow_type` is the canonical GDC path — the unprefixed
# `workflow_type` comes back as an unrecognized field warning. Flattened
# to a top-level `workflow_type` column at manifest-write time.
FILE_FIELDS: list[str] = [
    "file_id",
    "file_name",
    "file_size",
    "md5sum",
    "access",
    "data_category",
    "data_type",
    "data_format",
    "experimental_strategy",
    "analysis.workflow_type",
    # Needed to tell the methylation array generations apart: TCGA ships
    # Illumina 450k (9,812 files), 27k (2,662) and EPIC v2 (53), and their
    # probe sets differ, so a beta value is only comparable within a
    # platform. Harmless for every other modality.
    "platform",
    "cases.case_id",
    "cases.submitter_id",
    "cases.samples.portions.analytes.aliquots.aliquot_id",
    # GDC's uniform "what biospecimen entity is this file about" pointer.
    # The nested `cases.samples...aliquots` path above only resolves for
    # aliquot-level modalities; `associated_entities` names the entity at
    # whatever level the file actually attaches (aliquot for MAF / STAR,
    # sample for Pathology Report, portion for RPPA), so modalities that
    # aren't aliquot-scoped can still resolve an FK.
    "associated_entities.entity_id",
    "associated_entities.entity_type",
    "associated_entities.entity_submitter_id",
]


# Per-modality filter values that lock the `/files` request to the exact
# tooling we expect. Filtering on `data_type` alone is technically enough
# for TCGA today, but adding `data_format` + `experimental_strategy` +
# `analysis.workflow_type` guards against future GDC additions silently
# shipping different workflows under the same data_type. Values verified
# against live `/files` responses on 2026-05; if any of them changes
# upstream the filter will return zero hits and surface the regression
# loudly.
MODALITY_FILTERS: dict[str, dict[str, str]] = {
    "Masked Somatic Mutation": {
        "data_format": "MAF",
        "data_category": "Simple Nucleotide Variation",
        "experimental_strategy": "WXS",
        "analysis.workflow_type": "Aliquot Ensemble Somatic Variant Merging and Masking",
    },
    "Gene Expression Quantification": {
        "data_format": "TSV",
        "data_category": "Transcriptome Profiling",
        "experimental_strategy": "RNA-Seq",
        "analysis.workflow_type": "STAR - Counts",
    },
    # Pathology reports are scanned documents, not pipeline output, so they
    # carry no experimental_strategy or workflow_type to pin — data_format
    # and data_category are the whole lock.
    "Pathology Report": {
        "data_format": "PDF",
        "data_category": "Clinical",
    },
    # miRNA-Seq deliberately does NOT pin `data_format`: 11,082 files are
    # TXT and 359 are TSV (GBM / OV / LUSC, named `*.mirnaseq.*` rather than
    # `*.mirbase21.*` — the older Genome Analyzer platform). Both carry the
    # identical four-column layout, so pinning format would silently drop
    # 359 aliquots.
    "miRNA Expression Quantification": {
        "data_category": "Transcriptome Profiling",
        "experimental_strategy": "miRNA-Seq",
        "analysis.workflow_type": "BCGSC miRNA Profiling",
    },
    "Protein Expression Quantification": {
        "data_format": "TSV",
        "data_category": "Proteome Profiling",
        "experimental_strategy": "Reverse Phase Protein Array",
    },
    # The two copy-number segment modalities pin neither
    # `experimental_strategy` nor `analysis.workflow_type`, and that is
    # deliberate rather than an oversight:
    #
    #   - Allele-specific segments come from three callers (ASCAT2 and
    #     ASCAT3 on Genotyping Array, AscatNGS on WGS). All three are real
    #     TCGA content and all three ship. `workflow_type` is carried as a
    #     column instead, because ASCAT2 and ASCAT3 genuinely disagree —
    #     they fit ploidy independently, so the same aliquot can be called
    #     modal CN 2 by one and 3-4 by the other. A consumer must choose.
    #   - Masked segments are DNAcopy-only today, but the mask (germline
    #     CNV removal) is the meaningful lock and it's already implied by
    #     the data_type.
    "Allele-specific Copy Number Segment": {
        "data_format": "TXT",
        "data_category": "Copy Number Variation",
    },
    "Masked Copy Number Segment": {
        "data_format": "TXT",
        "data_category": "Copy Number Variation",
    },
    # The unmasked counterpart of the two above, and the one data_type whose
    # workflows are not interchangeable at the parser level:
    #
    #   - DNAcopy (22,629 files) writes the same six columns as the masked
    #     type over the very same aliquots — the only difference is that
    #     germline CNV segments are still present.
    #   - GATK4 CNV (10,722 files) is WGS rather than Genotyping Array, and
    #     writes `GDC_Aliquot_ID` holding a *barcode* where every other seg
    #     file writes `GDC_Aliquot` holding a UUID. Its chromosomes are
    #     `chr`-prefixed and it is a paired caller (two aliquot entities).
    #
    # So `analysis.workflow_type` is pinned by the caller via `extra_filters`
    # and each workflow lands in its own raw directory.
    "Copy Number Segment": {
        "data_format": "TXT",
        "data_category": "Copy Number Variation",
    },
    # All four workflows (ASCAT2, ASCAT3, AscatNGS, ABSOLUTE LiftOver) write
    # an identical eight-column projection onto the same GENCODE v36 gene
    # model, so unlike `Copy Number Segment` they share one directory and one
    # manifest, with `workflow_type` carried as a column — the same treatment
    # `Allele-specific Copy Number Segment` gives its three callers.
    #
    # ABSOLUTE LiftOver is the reason this data_type is fetched at all: it
    # ships no segment file anywhere in the GDC, so it is the only copy
    # number source that cannot be reconstructed from data we already hold.
    "Gene Level Copy Number": {
        "data_format": "TSV",
        "data_category": "Copy Number Variation",
    },
    # SeSAMe level-3 betas. One headerless two-column TXT per aliquot:
    # probe id, beta. `platform` is deliberately unpinned and carried as a
    # column instead — the three array generations have different probe
    # sets, so pinning one would silently drop the others.
    "Methylation Beta Value": {
        "data_format": "TXT",
        "data_category": "DNA Methylation",
        "experimental_strategy": "Methylation Array",
        "analysis.workflow_type": "SeSAMe Methylation Beta Estimation",
    },
    # miRNA isoform-level quantification — the per-isoform companion to
    # `miRNA Expression Quantification`, from the same BCGSC run. Like that
    # modality it does NOT pin `data_format`: 11,082 files are TXT and 359
    # are TSV (the older Genome Analyzer platform), with identical columns.
    # The two BCR supplement data types span several serializations of the
    # same underlying submission — a project-level `bcr biotab` TSV plus
    # per-case XML. `data_format` is therefore NOT pinned here; the caller
    # passes it, so each serialization lands in its own raw directory.
    # `data_category` is the only safe lock.
    "Clinical Supplement": {
        "data_category": "Clinical",
    },
    "Biospecimen Supplement": {
        "data_category": "Biospecimen",
    },
    "Isoform Expression Quantification": {
        "data_category": "Transcriptome Profiling",
        "experimental_strategy": "miRNA-Seq",
        "analysis.workflow_type": "BCGSC miRNA Profiling",
    },
}


def list_open_files(
    client: GDCClient,
    project_id: str,
    data_type: str,
    extra_filters: list[dict[str, Any]] | None = None,
) -> list[dict[str, Any]]:
    """Return all open-access /files hits for the given project + data_type.

    Filter clauses:
      - cases.project.project_id = <project_id>
      - access = "open"
      - data_type = <data_type>
      - …plus every key/value pair in `MODALITY_FILTERS[data_type]` (locks
        format, experimental_strategy, workflow_type for the modalities we
        explicitly support).
    """
    clauses = [
        eq("cases.project.project_id", project_id),
        eq("access", "open"),
        eq("data_type", data_type),
    ]
    for field, value in MODALITY_FILTERS.get(data_type, {}).items():
        clauses.append(eq(field, value))
    if extra_filters:
        clauses.extend(extra_filters)
    return client.files(filters=and_(*clauses), fields=FILE_FIELDS, page_size=500)


def fetch_files(
    client: GDCClient,
    project_id: str,
    data_type: str,
    out_dir: Path,
    extra_filters: list[dict[str, Any]] | None = None,
    skip_existing: bool = True,
    batch_size: int = 50,
    max_files: int | None = None,
) -> list[dict[str, Any]]:
    """Download hits' bytes to <out_dir>/<file_name> and write manifest.json.

    Returns the manifest list. Existing files (matching size) are skipped
    unless skip_existing=False. Downloads are batched via POST /data
    (`batch_size` per request) so wall time scales with total bytes rather
    than per-file latency.

    `max_files` caps the *total* number of files we'll have on disk after
    this call (cached + freshly downloaded). Set to 0 to discover and
    write the full manifest without downloading any new bytes — useful for
    populating the `files` table across projects where you only want a
    sampler of the actual molecular content. The manifest always lists
    every hit returned by `/files`, regardless of how many were
    downloaded; entries trimmed off by `max_files` are tagged
    `_status="manifest_only"`.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    hits = list_open_files(client, project_id, data_type, extra_filters)

    cached_ids: set[str] = set()
    eligible_for_download: list[str] = []
    for hit in hits:
        target = out_dir / hit["file_name"]
        if skip_existing and target.exists() and target.stat().st_size == hit.get("file_size"):
            cached_ids.add(hit["file_id"])
        else:
            eligible_for_download.append(hit["file_id"])

    if max_files is None:
        to_download_ids = eligible_for_download
    else:
        # `max_files` is a target *total* count on disk, so subtract what's
        # already cached. If we already have ≥ max_files cached, download
        # nothing new.
        budget = max(0, max_files - len(cached_ids))
        to_download_ids = eligible_for_download[:budget]
    skipped_ids = set(eligible_for_download) - set(to_download_ids)

    if len(to_download_ids) >= 2:
        client.bulk_download(to_download_ids, out_dir, batch_size=batch_size)
    elif len(to_download_ids) == 1:
        only = next(h for h in hits if h["file_id"] == to_download_ids[0])
        client.download(only["file_id"], out_dir / only["file_name"])

    # Per-file version provenance. The `/files` + `/data` endpoints only
    # ever serve the current GDC release, so the release we fetch at is a
    # timestamp, not a description of the bytes. `gdc_version` +
    # `gdc_first_release` say what the bytes actually are, and
    # `gdc_superseded` flags the one case that matters — a file GDC has
    # since replaced with a newer version under a different id.
    versions = client.versions([h["file_id"] for h in hits]) if hits else {}

    manifest: list[dict[str, Any]] = []
    for hit in hits:
        version = versions.get(hit["file_id"], {})
        # GDC nests workflow info under `analysis.*`; flatten it into the
        # manifest so downstream consumers see a single `workflow_type`
        # column matching the dataset's `files` schema.
        analysis = hit.get("analysis") or {}
        if hit["file_id"] in cached_ids:
            status = "cached"
        elif hit["file_id"] in skipped_ids:
            status = "manifest_only"
        else:
            status = "downloaded"
        manifest.append(
            {
                "file_id": hit["file_id"],
                "file_name": hit["file_name"],
                "file_size": hit.get("file_size"),
                "md5sum": hit.get("md5sum"),
                "data_category": hit.get("data_category"),
                "data_type": hit.get("data_type"),
                "data_format": hit.get("data_format"),
                "experimental_strategy": hit.get("experimental_strategy"),
                "platform": hit.get("platform"),
                "workflow_type": analysis.get("workflow_type"),
                "access": hit.get("access"),
                "gdc_version": version.get("version"),
                "gdc_first_release": version.get("release"),
                "gdc_superseded": bool(
                    version and version.get("latest_id") not in (None, hit["file_id"])
                ),
                "cases": hit.get("cases", []),
                "associated_entities": hit.get("associated_entities", []),
                "_status": status,
            }
        )

    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    return manifest


def fetch_file_index(client: GDCClient, project_id: str, out_path: Path) -> list[dict[str, Any]]:
    """Write the project's complete open-access file list to `out_path`.

    This is what lets the `files` table describe the project's whole GDC
    footprint rather than only the parts we downloaded. It is a listing, not
    a download: no bytes are fetched, so it costs one paged `/files` query
    plus the version lookups regardless of how large the excluded modalities
    are (TCGA-CHOL's 110 whole-slide images are 88.6 GiB and appear here as
    110 rows).

    Controlled-access files are deliberately not listed. A row carrying a
    download URL nobody reading an open dataset can use would be noise, and
    `cases.summary.data_categories` already reports that controlled data
    exists for a case.
    """
    clauses = [
        eq("cases.project.project_id", project_id),
        eq("access", "open"),
    ]
    hits = client.files(filters=and_(*clauses), fields=FILE_FIELDS, page_size=500)
    versions = client.versions([h["file_id"] for h in hits]) if hits else {}

    index: list[dict[str, Any]] = []
    for hit in hits:
        version = versions.get(hit["file_id"], {})
        analysis = hit.get("analysis") or {}
        index.append(
            {
                "file_id": hit["file_id"],
                "file_name": hit["file_name"],
                "file_size": hit.get("file_size"),
                "md5sum": hit.get("md5sum"),
                "data_category": hit.get("data_category"),
                "data_type": hit.get("data_type"),
                "data_format": hit.get("data_format"),
                "experimental_strategy": hit.get("experimental_strategy"),
                "platform": hit.get("platform"),
                "workflow_type": analysis.get("workflow_type"),
                "access": hit.get("access"),
                "gdc_version": version.get("version"),
                "gdc_first_release": version.get("release"),
                "gdc_superseded": bool(
                    version and version.get("latest_id") not in (None, hit["file_id"])
                ),
                "cases": hit.get("cases", []),
            }
        )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(index, indent=2))
    return index
