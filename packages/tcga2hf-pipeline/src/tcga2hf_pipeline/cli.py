from __future__ import annotations

import json
import os
import shutil
from collections.abc import Callable
from pathlib import Path
from typing import Annotated

import typer
from dotenv import load_dotenv
from tcga2hf.schema import TABULAR_TABLES

from tcga2hf_pipeline import (
    biospecimen_supplement,
    cdr,
    clinical,
    clinical_supplement,
    copy_number,
    dataset_card,
    expression,
    expression_dataset,
    genomic,
    hf_upload,
    mirna,
    msigdb,
    mutations,
    pathology,
    protein_expression,
    ssgsea,
    survival,
    tabular,
    verify,
)
from tcga2hf_pipeline import webdataset as wds_mod
from tcga2hf_pipeline.gdc import GDC_BASE_URL, GDCClient, write_cases_json
from tcga2hf_pipeline.gdc import eq as gdc_eq

# Load .env from cwd (or any parent) on import. override=True so the project's
# .env wins over any inherited shell variable: the HF_TOKEN here is scoped to
# this project, and we don't want a stale global token to silently take over.
load_dotenv(override=True)

app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help="Download public TCGA data from the NCI GDC and stage it for the HF Hub.",
)

DEFAULT_DATA_DIR = Path.home() / "data" / "tcga2hf"


def _resolve_data_dir(data_dir: Path | None) -> Path:
    if data_dir is not None:
        return data_dir
    env = os.environ.get("TCGA2HF_DATA_DIR")
    if env:
        return Path(env).expanduser()
    return DEFAULT_DATA_DIR


DataDirOpt = Annotated[
    Path | None,
    typer.Option(
        "--data-dir",
        help="Root data dir. Defaults to $TCGA2HF_DATA_DIR or $HOME/data/tcga2hf.",
    ),
]

ProjectFilterOpt = Annotated[
    list[str] | None,
    typer.Option(
        "--project",
        help=(
            "TCGA project id (repeatable). If omitted, all projects under "
            "<data-dir>/raw/ are processed and the output tree is rebuilt from "
            "scratch. If given, only those projects' output dirs are replaced — "
            "other projects' parquets are left untouched."
        ),
    ),
]


def _clear_processed(processed_dir: Path, project: list[str] | None) -> None:
    """Remove prior output for the projects about to be rebuilt.

    Without `--project` the whole tree goes, so projects dropped from raw/
    and any legacy layout don't leave stale files behind. With `--project`
    only the named projects' directories go: a full wipe there would delete
    parquets we aren't regenerating, and (because HF uploads diff against
    what's on disk) would turn an incremental append into a full re-upload.
    """
    if project is None:
        if processed_dir.exists():
            shutil.rmtree(processed_dir)
        processed_dir.mkdir(parents=True)
        return
    processed_dir.mkdir(parents=True, exist_ok=True)
    for proj in project:
        target = processed_dir / proj
        if target.exists():
            shutil.rmtree(target)


def _card_inputs(
    raw_dir: Path, processed_dir: Path, built_marker: str
) -> tuple[list[str], dict[str, str]]:
    """Return (projects, gdc_releases) for every project present in the output tree.

    The dataset card describes the tree as a whole, not just the projects
    this invocation rebuilt — otherwise a `--project` build would drop the
    other projects from the `configs:` block and orphan their parquets on
    the Hub. `built_marker` is a glob, relative to a project dir, that
    identifies a completed build for the layout being written.

    Each project's GDC release is read back from its own raw
    `gdc_status.json`, so projects fetched at different releases report
    their own — the card renders a per-project list rather than a single
    release when they diverge.
    """
    projects: list[str] = []
    gdc_releases: dict[str, str] = {}
    if not processed_dir.exists():
        return projects, gdc_releases
    for project_dir in sorted(p for p in processed_dir.iterdir() if p.is_dir()):
        if not any(project_dir.glob(built_marker)):
            continue
        proj = project_dir.name
        projects.append(proj)
        status_path = raw_dir / proj / "gdc_status.json"
        if status_path.exists():
            gdc_releases[proj] = json.loads(status_path.read_text()).get(
                "data_release", "<unknown>"
            )
    return projects, gdc_releases


def _write_ssgsea_stats(processed_dir: Path, only: set[str] | None) -> None:
    """Build each ssGSEA collection's stats table from the written scores.

    Skipped for collections whose stats table wasn't requested, and a no-op
    when no project has scores on disk yet.
    """
    import pyarrow as pa
    import pyarrow.parquet as pq
    from tcga2hf.schema import SSGSEA_COLLECTIONS, TABULAR_TABLES, ssgsea_stats_table

    for coll in SSGSEA_COLLECTIONS:
        table_name = ssgsea_stats_table(coll)
        if only is not None and table_name not in only:
            continue
        by_project = tabular.ssgsea_stats_rows(processed_dir, coll)
        if not by_project:
            continue
        for proj, rows in by_project.items():
            out_path = processed_dir / proj / table_name / "data.parquet"
            out_path.parent.mkdir(parents=True, exist_ok=True)
            pq.write_table(
                pa.Table.from_pylist(rows, schema=TABULAR_TABLES[table_name]),
                out_path,
                write_page_index=True,
            )
        n = sum(len(r) for r in by_project.values())
        typer.echo(f"  {table_name}: {n:,} rows across {len(by_project)} projects")


def _select_projects(raw_dir: Path, project: list[str] | None) -> list[Path]:
    """Resolve `--project` to the cases.json paths to process (all if None)."""
    project_files = sorted(raw_dir.glob("*/cases.json"))
    if project:
        wanted = set(project)
        project_files = [p for p in project_files if p.parent.name in wanted]
        missing = wanted - {p.parent.name for p in project_files}
        if missing:
            raise typer.BadParameter(f"requested projects missing from raw/: {sorted(missing)}")
    if not project_files:
        raise typer.BadParameter(f"no cases.json files found under {raw_dir}.")
    return project_files


@app.command("fetch-clinical")
def fetch_clinical_cmd(
    project: Annotated[
        list[str],
        typer.Option(
            "--project",
            help="TCGA project id (repeatable), e.g. --project TCGA-CHOL --project TCGA-DLBC.",
        ),
    ],
    data_dir: DataDirOpt = None,
) -> None:
    """Fetch clinical case JSON for one or more TCGA projects from the GDC."""
    import hashlib

    root = _resolve_data_dir(data_dir)
    raw_dir = root / "raw"
    typer.echo(f"raw dir: {raw_dir}")

    with GDCClient() as client:
        # Snapshot the GDC server's data release once at the top so every project
        # fetched in this invocation pairs with the same release tag.
        status = client.status()
        release_str = status.get("data_release", "<unknown>")
        typer.echo(f"GDC: {release_str} (tag {status.get('tag', '?')})")

        # Capture the dictionary alongside, tagged by release. Reuses the cached
        # copy if a fetch with the same release already happened locally.
        version = status.get("data_release_version") or {}
        major = version.get("major", "x")
        minor = version.get("minor", "y")
        dict_path = raw_dir / f"gdc_dictionary.{major}.{minor}.json"
        if not dict_path.exists():
            typer.echo(f"capturing GDC dictionary -> {dict_path.name} ...")
            dict_path.parent.mkdir(parents=True, exist_ok=True)
            dictionary = client.dictionary()
            dict_text = json.dumps(dictionary, indent=2, sort_keys=True)
            dict_path.write_text(dict_text)
        else:
            dict_text = dict_path.read_text()
        # SHA-256 of the (sorted) dictionary JSON pins provenance precisely
        # without needing to diff 10MB blobs across runs.
        status["_gdc_dictionary_sha256"] = hashlib.sha256(dict_text.encode()).hexdigest()
        status["_gdc_dictionary_path"] = dict_path.name

        for proj in project:
            typer.echo(f"fetching {proj} ...")
            cases = clinical.fetch_clinical([proj], client)
            out_path = raw_dir / proj / "cases.json"
            write_cases_json(cases, out_path)
            (raw_dir / proj / "gdc_status.json").write_text(json.dumps(status, indent=2))
            typer.echo(f"  wrote {len(cases)} cases -> {out_path}")


def _fetch_modality(
    project: list[str],
    data_dir: Path | None,
    data_type: str,
    modality_dir: str,
    max_files: int | None = None,
    workflow_type: str | None = None,
    data_format: str | None = None,
) -> None:
    """Shared body for fetch-mutations / fetch-expression / future modalities.

    The GDC release is recorded per modality, in the modality's own
    directory — not on the project. Modalities are fetched at different
    times (a new one added years later pulls whatever release GDC is
    serving then), so one status file per project would report only the
    most recent fetch and silently overwrite the release + dictionary hash
    that `fetch-clinical` recorded for `cases.json`.

    `workflow_type` narrows the fetch to one GDC `analysis.workflow_type`
    within the data type. It exists for `Copy Number Segment`, whose two
    workflows write incompatible headers and so need one raw directory
    each; a data type whose workflows share a parser (allele-specific
    segments, gene-level copy number) leaves it None and keeps them
    together under a single manifest.

    `data_format` does the same for the BCR supplements, where one
    `data_type` spans a project-level biotab and several per-case XML
    serializations that have nothing in common but their name.
    """
    root = _resolve_data_dir(data_dir)
    raw_dir = root / "raw"
    typer.echo(f"raw dir: {raw_dir}")

    with GDCClient() as client:
        status = client.status()
        typer.echo(f"GDC: {status.get('data_release', '<unknown>')} (tag {status.get('tag', '?')})")
        for proj in project:
            out_dir = raw_dir / proj / modality_dir
            label = f"{data_type!r}" + (f" [{workflow_type}]" if workflow_type else "")
            typer.echo(f"fetching {label} for {proj} -> {out_dir}")
            extra_filters = [
                clause
                for field, value in (
                    ("analysis.workflow_type", workflow_type),
                    ("data_format", data_format),
                )
                if value
                for clause in [gdc_eq(field, value)]
            ] or None
            manifest = genomic.fetch_files(
                client,
                proj,
                data_type,
                out_dir,
                extra_filters=extra_filters,
                max_files=max_files,
            )
            n_dl = sum(1 for m in manifest if m["_status"] == "downloaded")
            n_cache = sum(1 for m in manifest if m["_status"] == "cached")
            n_skip = sum(1 for m in manifest if m["_status"] == "manifest_only")
            total_mb = sum((m.get("file_size") or 0) for m in manifest) / 1e6
            extra = f", {n_skip} manifest-only" if n_skip else ""
            typer.echo(
                f"  {len(manifest):>4} files ({n_dl} downloaded, {n_cache} cached{extra}), "
                f"{total_mb:.1f} MB total in manifest"
            )
            (out_dir / "gdc_status.json").write_text(json.dumps(status, indent=2))


MaxFilesOpt = Annotated[
    int | None,
    typer.Option(
        "--max-files",
        help=(
            "Cap on total files on disk (cached + new) per project. Set to 1 to "
            "sample one file per project, or 0 to populate the manifest without "
            "downloading any new bytes. Manifest always lists every discovered file."
        ),
    ),
]


@app.command("fetch-file-index")
def fetch_file_index_cmd(
    project: Annotated[
        list[str],
        typer.Option("--project", help="TCGA project id (repeatable)."),
    ],
    data_dir: DataDirOpt = None,
) -> None:
    """List every open-access GDC file for a project (no bytes downloaded).

    Writes `<data-dir>/raw/<project>/files_index.json`, which lets the
    `files` table carry a row for **every** open file the GDC holds — not
    only the ones this pipeline downloads — each with a `gdc_download_url`
    and an `in_dataset` flag.

    That is the compromise between completeness and size: TCGA-CHOL's 110
    whole-slide images are 88.6 GiB of bytes but 110 rows of index, so the
    dataset can describe its own scope honestly without carrying them.
    """
    root = _resolve_data_dir(data_dir)
    raw_dir = root / "raw"
    with GDCClient() as client:
        status = client.status()
        typer.echo(f"GDC: {status.get('data_release', '<unknown>')}")
        for proj in project:
            out_path = raw_dir / proj / "files_index.json"
            index = genomic.fetch_file_index(client, proj, out_path)
            by_type: dict[str, int] = {}
            for entry in index:
                by_type[entry["data_type"]] = by_type.get(entry["data_type"], 0) + 1
            total_gb = sum((e.get("file_size") or 0) for e in index) / 1e9
            typer.echo(
                f"  {proj:<12} {len(index):>5} open files, {len(by_type)} data types, "
                f"{total_gb:.1f} GB at GDC -> {out_path.name}"
            )


@app.command("fetch-mutations")
def fetch_mutations_cmd(
    project: Annotated[
        list[str],
        typer.Option(
            "--project",
            help="TCGA project id (repeatable). Downloads open-access Masked Somatic Mutation MAFs.",  # noqa: E501
        ),
    ],
    data_dir: DataDirOpt = None,
    max_files: MaxFilesOpt = None,
) -> None:
    """Download open-access Masked Somatic Mutation MAF files (DNA)."""
    _fetch_modality(project, data_dir, "Masked Somatic Mutation", "mutations", max_files=max_files)


@app.command("fetch-expression")
def fetch_expression_cmd(
    project: Annotated[
        list[str],
        typer.Option(
            "--project",
            help="TCGA project id (repeatable). Downloads open-access RNA-Seq STAR gene-count TSVs.",  # noqa: E501
        ),
    ],
    data_dir: DataDirOpt = None,
    max_files: MaxFilesOpt = None,
) -> None:
    """Download open-access Gene Expression Quantification TSVs (RNA-Seq, STAR counts)."""
    _fetch_modality(
        project, data_dir, "Gene Expression Quantification", "expression", max_files=max_files
    )


@app.command("fetch-pathology-reports")
def fetch_pathology_reports_cmd(
    project: Annotated[
        list[str],
        typer.Option(
            "--project",
            help="TCGA project id (repeatable). Downloads open-access Pathology Report PDFs.",
        ),
    ],
    data_dir: DataDirOpt = None,
    max_files: MaxFilesOpt = None,
) -> None:
    """Download open-access Pathology Report PDFs (scanned BCR documents).

    One report per case's tumor sample, ~11,200 across TCGA (~2.6 GB; LAML
    has none). The PDFs are shipped into the dataset verbatim — see
    `tcga2hf_pipeline.pathology` for why no text extraction happens here.
    """
    _fetch_modality(project, data_dir, "Pathology Report", "pathology_reports", max_files=max_files)


@app.command("fetch-mirna")
def fetch_mirna_cmd(
    project: Annotated[
        list[str],
        typer.Option(
            "--project",
            help="TCGA project id (repeatable). Downloads open-access miRNA-Seq quantifications.",
        ),
    ],
    data_dir: DataDirOpt = None,
    max_files: MaxFilesOpt = None,
) -> None:
    """Download open-access miRNA Expression Quantification files (miRNA-Seq).

    One file per aliquot, ~1,881 miRBase v21 mature miRNAs each; 11,441
    across TCGA (~580 MB). Isoform-level quantification is a separate,
    much larger modality (4 GB) that we don't fetch.
    """
    _fetch_modality(
        project, data_dir, "miRNA Expression Quantification", "mirna", max_files=max_files
    )


@app.command("fetch-protein-expression")
def fetch_protein_expression_cmd(
    project: Annotated[
        list[str],
        typer.Option(
            "--project",
            help="TCGA project id (repeatable). Downloads open-access RPPA protein expression.",
        ),
    ],
    data_dir: DataDirOpt = None,
    max_files: MaxFilesOpt = None,
) -> None:
    """Download open-access Protein Expression Quantification files (RPPA).

    Reverse Phase Protein Array, ~487 antibodies per portion; 7,906 files
    across TCGA (~172 MB) covering 7,827 cases — the narrowest of the
    molecular modalities, since RPPA was only run on a subset.
    """
    _fetch_modality(
        project,
        data_dir,
        "Protein Expression Quantification",
        "protein_expression",
        max_files=max_files,
    )


@app.command("fetch-methylation")
def fetch_methylation_cmd(
    project: Annotated[
        list[str],
        typer.Option(
            "--project",
            help="TCGA project id (repeatable). Downloads SeSAMe methylation beta values.",
        ),
    ],
    data_dir: DataDirOpt = None,
    max_files: MaxFilesOpt = None,
) -> None:
    """Download open-access Methylation Beta Value files (12,527 pan-TCGA).

    SeSAMe level-3 betas from the Illumina methylation arrays: one
    headerless two-column TXT per aliquot, probe id and beta in [0, 1].
    A 450k file is 486,427 probes, of which ~15% are masked and written
    `NA` — a real "not trustworthy" rather than zero.

    Three array generations ship and their probe sets differ (450k: 9,812
    files, 27k: 2,662, EPIC v2: 53), so `platform` is carried as a column
    and betas are only comparable within one.

    Not fetched: `Masked Intensities`, the raw two-channel IDATs these betas
    are computed from (25,054 files pan-TCGA, binary).
    """
    _fetch_modality(
        project, data_dir, "Methylation Beta Value", "methylation", max_files=max_files
    )


@app.command("fetch-mirna-isoform")
def fetch_mirna_isoform_cmd(
    project: Annotated[
        list[str],
        typer.Option(
            "--project",
            help="TCGA project id (repeatable). Downloads miRNA isoform quantification.",
        ),
    ],
    data_dir: DataDirOpt = None,
    max_files: MaxFilesOpt = None,
) -> None:
    """Download open-access Isoform Expression Quantification files.

    The per-isoform companion to `fetch-mirna`, from the same BCGSC run over
    the same aliquots: where that modality gives one number per mature
    miRNA, this splits it across the distinct read pileups collapsed into
    it (~4,500 isoforms per aliquot against ~1,881 mature miRNAs).
    """
    _fetch_modality(
        project,
        data_dir,
        "Isoform Expression Quantification",
        "mirna_isoform",
        max_files=max_files,
    )


@app.command("fetch-copy-number")
def fetch_copy_number_cmd(
    project: Annotated[
        list[str],
        typer.Option(
            "--project",
            help="TCGA project id (repeatable). Downloads open-access copy number files.",
        ),
    ],
    data_dir: DataDirOpt = None,
    max_files: MaxFilesOpt = None,
    skip_gene_level: Annotated[
        bool,
        typer.Option(
            "--skip-gene-level",
            help=(
                "Fetch only the four segment-level sources. Gene-level copy "
                "number is ~97% of the bytes (10.6 GiB for TCGA-BRCA)."
            ),
        ),
    ] = False,
) -> None:
    """Download every open-access copy number file GDC serves for TCGA.

    The GDC exposes copy number under four open `data_type`s spanning six
    workflows, and this command covers all of them — a one-to-one mapping
    between the source and the published tables:

      - Allele-specific Copy Number Segment (23,225 files, ~215 MB) —
        integer total / major / minor copy number per segment, from ASCAT2,
        ASCAT3 and AscatNGS.
      - Masked Copy Number Segment (22,629 files, ~344 MB) — DNAcopy log2
        ratio segment means with germline CNVs masked out.
      - Copy Number Segment (33,351 files, ~2.3 GB) — the unmasked calls.
        DNAcopy (22,629 files) covers the very same aliquots as the masked
        type with the germline segments still in; GATK4 CNV (10,722 files)
        is WGS coverage over aliquots the genotyping arrays never touched.
      - Gene Level Copy Number (33,902 files, ~109 GB) — the calls projected
        onto GENCODE v36, from all four callers.

    Nothing open-access is left behind: `Raw Intensities` and `Intermediate
    Analysis Archive` are the only other copy number data types, and both
    are entirely controlled-access.

    On gene-level bytes. The ASCAT2 / ASCAT3 / AscatNGS projections are
    exactly reproducible from the allele-specific segments (verified on
    TCGA-CHOL: zero mismatches across 3 aliquots / 180k genes, with
    min/max copy number for boundary-spanning genes recovered as the
    min/max over overlapping segments), so they are redundant with data we
    already hold. ABSOLUTE LiftOver is not: it ships no segment file
    anywhere in the GDC, so its purity- and ploidy-corrected absolute copy
    number exists only at gene level. All four are fetched anyway, because
    the goal is source parity rather than a minimal spanning set — and the
    34 GB per workflow is TSV bloat, not information. The gene model is
    byte-identical across all 33,902 files and copy number runs in long
    constant stretches, so the same calls re-encode to roughly 1.3 KB per
    aliquot as parallel arrays, against GDC's 3,350 KB of TSV.
    """
    _fetch_modality(
        project,
        data_dir,
        "Allele-specific Copy Number Segment",
        "copy_number_allele_specific",
        max_files=max_files,
    )
    _fetch_modality(
        project,
        data_dir,
        "Masked Copy Number Segment",
        "copy_number_masked",
        max_files=max_files,
    )
    # One raw dir per workflow: DNAcopy writes `GDC_Aliquot` (a UUID) and
    # bare chromosome names, GATK4 CNV writes `GDC_Aliquot_ID` (a barcode)
    # and `chr`-prefixed names, so they cannot share a parser or a manifest.
    _fetch_modality(
        project,
        data_dir,
        "Copy Number Segment",
        "copy_number_segment_dnacopy",
        max_files=max_files,
        workflow_type="DNAcopy",
    )
    _fetch_modality(
        project,
        data_dir,
        "Copy Number Segment",
        "copy_number_segment_gatk4",
        max_files=max_files,
        workflow_type="GATK4 CNV",
    )
    if skip_gene_level:
        typer.echo("skipping Gene Level Copy Number (--skip-gene-level)")
        return
    # All four gene-level workflows share one directory: identical columns
    # over an identical gene model, distinguished by `workflow_type`.
    _fetch_modality(
        project,
        data_dir,
        "Gene Level Copy Number",
        "gene_level_copy_number",
        max_files=max_files,
    )


@app.command("fetch-bcr-xml")
def fetch_bcr_xml_cmd(
    project: Annotated[
        list[str],
        typer.Option(
            "--project",
            help="TCGA project id (repeatable). Downloads the per-case BCR XML supplements.",
        ),
    ],
    data_dir: DataDirOpt = None,
    max_files: MaxFilesOpt = None,
) -> None:
    """Download the per-case BCR XML Clinical / Biospecimen Supplements.

    The `Clinical Supplement` and `Biospecimen Supplement` data types each
    ship in several serializations, and `fetch-clinical-supplements` /
    `fetch-biospecimen-supplements` take only the `bcr biotab` ones — the
    project-level TSV forms this pipeline parses into columns. That leaves
    the majority of the files behind: for TCGA-CHOL, 7 of 65 clinical and
    10 of 112 biospecimen supplement files.

    The rest are **per-case XML**, one file per patient, and this fetches
    all four kinds:

      - Clinical `bcr xml`      — the patient's full BCR clinical record
      - Clinical `bcr omf xml`  — the "other malignancy" form
      - Biospecimen `bcr xml`   — the full specimen chain
      - Biospecimen `bcr ssf xml` — site-specific factors

    They are small (5.4 MB for all 160 CHOL files) and are shipped as the
    XML text verbatim rather than parsed: the biotab tables already give a
    parsed view of the same underlying BCR data, so parsing again would be
    deriving a second time rather than recording what TCGA holds.
    """
    for label, dtype, fmt, out in (
        ("clinical", "Clinical Supplement", "bcr xml", "clinical_supplement_xml"),
        ("clinical omf", "Clinical Supplement", "bcr omf xml", "clinical_supplement_omf_xml"),
        ("biospecimen", "Biospecimen Supplement", "bcr xml", "biospecimen_supplement_xml"),
        (
            "biospecimen ssf",
            "Biospecimen Supplement",
            "bcr ssf xml",
            "biospecimen_supplement_ssf_xml",
        ),
    ):
        typer.echo(f"--- {label} ({fmt}) ---")
        _fetch_modality(
            project, data_dir, dtype, out, max_files=max_files, data_format=fmt
        )


@app.command("fetch-clinical-supplements")
def fetch_clinical_supplements_cmd(
    project: Annotated[
        list[str],
        typer.Option(
            "--project",
            help="TCGA project id (repeatable). Downloads BCR biotab Clinical Supplement files.",
        ),
    ],
    data_dir: DataDirOpt = None,
) -> None:
    """Download BCR biotab Clinical Supplement files (the original BCR forms).

    These TSV files carry clinical fields the harmonized `/cases?expand=...`
    API drops or under-populates — most importantly Liu et al. 2018's
    `treatment_outcome_first_course`, the disease-free signal that drives
    DFI re-derivation. Files land at `<data-dir>/raw/<project>/clinical_supplement/`.
    """
    root = _resolve_data_dir(data_dir)
    raw_dir = root / "raw"
    typer.echo(f"raw dir: {raw_dir}")

    with GDCClient() as client:
        status = client.status()
        typer.echo(f"GDC: {status.get('data_release', '<unknown>')} (tag {status.get('tag', '?')})")
        for proj in project:
            out_dir = raw_dir / proj / "clinical_supplement"
            typer.echo(f"fetching Clinical Supplement biotabs for {proj} -> {out_dir}")
            manifest = clinical_supplement.fetch_clinical_supplements(client, proj, out_dir)
            n_dl = sum(1 for m in manifest if m["_status"] == "downloaded")
            n_cache = sum(1 for m in manifest if m["_status"] == "cached")
            kinds = sorted({m["form_kind"] for m in manifest})
            typer.echo(
                f"  {len(manifest):>3} files ({n_dl} downloaded, {n_cache} cached)  "
                f"forms: {', '.join(kinds)}"
            )


@app.command("fetch-biospecimen-supplements")
def fetch_biospecimen_supplements_cmd(
    project: Annotated[
        list[str],
        typer.Option(
            "--project",
            help="TCGA project id (repeatable). Downloads BCR biotab Biospecimen Supplement files.",
        ),
    ],
    data_dir: DataDirOpt = None,
) -> None:
    """Download BCR biotab Biospecimen Supplement files (the specimen chain).

    The counterpart to the clinical biotabs: slide-level tumour-nuclei and
    necrosis percentages, analyte quality ratios, plate / shipment / centre
    provenance, and the disease-specific site-specific-factor forms. 340
    files across TCGA (~76 MB). Files land at
    `<data-dir>/raw/<project>/biospecimen_supplement/`.
    """
    root = _resolve_data_dir(data_dir)
    raw_dir = root / "raw"
    typer.echo(f"raw dir: {raw_dir}")

    with GDCClient() as client:
        status = client.status()
        typer.echo(f"GDC: {status.get('data_release', '<unknown>')} (tag {status.get('tag', '?')})")
        for proj in project:
            out_dir = raw_dir / proj / "biospecimen_supplement"
            typer.echo(f"fetching Biospecimen Supplement biotabs for {proj} -> {out_dir}")
            manifest = biospecimen_supplement.fetch_biospecimen_supplements(client, proj, out_dir)
            n_dl = sum(1 for m in manifest if m["_status"] == "downloaded")
            n_cache = sum(1 for m in manifest if m["_status"] == "cached")
            kinds = sorted({m["form_kind"] for m in manifest})
            typer.echo(
                f"  {len(manifest):>3} files ({n_dl} downloaded, {n_cache} cached)  "
                f"forms: {', '.join(kinds)}"
            )


@app.command("fetch-project")
def fetch_project_cmd(
    project: Annotated[
        list[str],
        typer.Option("--project", help="TCGA project id (repeatable)."),
    ],
    data_dir: DataDirOpt = None,
    skip_gene_level: Annotated[
        bool,
        typer.Option("--skip-gene-level", help="Skip Gene Level Copy Number (~97% of CNV bytes)."),
    ] = False,
) -> None:
    """Fetch every modality this pipeline uses, for one or more projects.

    The individual `fetch-*` commands exist for re-fetching one thing; this
    runs all of them in the order a fresh project needs, so working through
    the cohort is one command per project rather than thirteen.

    Every step skips files already on disk, so re-running is cheap and this
    is also the way to bring an older project tree up to date after a new
    modality is added.

    Order matters in one place: `fetch-clinical` writes `cases.json`, which
    every build step joins against, and `fetch-file-index` records what the
    GDC holds so the `files` table can describe what we chose not to carry.
    Both are cheap and come first.
    """
    root = _resolve_data_dir(data_dir)
    for proj in project:
        typer.echo(f"\n{'=' * 62}\n{proj}\n{'=' * 62}")
        steps: list[tuple[str, Callable[[], None]]] = [
            (
                "clinical case tree",
                lambda p=proj: fetch_clinical_cmd(project=[p], data_dir=data_dir),
            ),
            ("file index", lambda p=proj: fetch_file_index_cmd(project=[p], data_dir=data_dir)),
            (
                "somatic mutations",
                lambda p=proj: fetch_mutations_cmd(project=[p], data_dir=data_dir, max_files=None),
            ),
            (
                "gene expression",
                lambda p=proj: fetch_expression_cmd(project=[p], data_dir=data_dir, max_files=None),
            ),
            (
                "miRNA",
                lambda p=proj: fetch_mirna_cmd(project=[p], data_dir=data_dir, max_files=None),
            ),
            (
                "miRNA isoforms",
                lambda p=proj: fetch_mirna_isoform_cmd(
                    project=[p], data_dir=data_dir, max_files=None
                ),
            ),
            (
                "protein expression",
                lambda p=proj: fetch_protein_expression_cmd(
                    project=[p], data_dir=data_dir, max_files=None
                ),
            ),
            (
                "methylation",
                lambda p=proj: fetch_methylation_cmd(
                    project=[p], data_dir=data_dir, max_files=None
                ),
            ),
            (
                "copy number",
                lambda p=proj: fetch_copy_number_cmd(
                    project=[p],
                    data_dir=data_dir,
                    max_files=None,
                    skip_gene_level=skip_gene_level,
                ),
            ),
            (
                "pathology reports",
                lambda p=proj: fetch_pathology_reports_cmd(
                    project=[p], data_dir=data_dir, max_files=None
                ),
            ),
            (
                "clinical supplements",
                lambda p=proj: fetch_clinical_supplements_cmd(project=[p], data_dir=data_dir),
            ),
            (
                "biospecimen supplements",
                lambda p=proj: fetch_biospecimen_supplements_cmd(project=[p], data_dir=data_dir),
            ),
            (
                "BCR XML supplements",
                lambda p=proj: fetch_bcr_xml_cmd(project=[p], data_dir=data_dir, max_files=None),
            ),
        ]
        for label, step in steps:
            typer.echo(f"\n--- {label} ---")
            step()
        raw = root / "raw" / proj
        total = sum(f.stat().st_size for f in raw.rglob("*") if f.is_file())
        typer.echo(f"\n{proj}: raw tree now {total / 2**30:.2f} GiB at {raw}")


@app.command("fetch-cdr")
def fetch_cdr_cmd(
    data_dir: DataDirOpt = None,
) -> None:
    """Fetch the Liu et al. 2018 PanCanAtlas CDR workbook (curated TCGA survival).

    Downloads the GDC's PanCanAtlas auxiliary file (UUID
    1b5f413e-...) into <data-dir>/raw/cdr/. md5-pinned; safe to re-run.
    The file is frozen at the 2018 data freeze, so coverage stops there
    — post-2018 cases get `cdr_matched=False` at build time.
    """
    root = _resolve_data_dir(data_dir)
    raw_dir = root / "raw"
    typer.echo(f"raw dir: {raw_dir}")
    out_path = cdr.fetch_cdr_workbook(raw_dir)
    typer.echo(f"  wrote -> {out_path}")
    index = cdr.load_cdr_index(raw_dir)
    typer.echo(f"  CDR rows indexed: {len(index)}")


@app.command("fetch-msigdb")
def fetch_msigdb_cmd(
    collection: Annotated[
        list[str] | None,
        typer.Option(
            "--collection",
            help=(
                "MSigDB collection (repeatable). Defaults to hallmark. "
                f"Known: {', '.join(sorted(msigdb.COLLECTIONS))}."
            ),
        ),
    ] = None,
    data_dir: DataDirOpt = None,
) -> None:
    """Fetch md5-pinned MSigDB gene-set collections (GMT) for ssGSEA scoring.

    Lands in `<data-dir>/raw/msigdb/`. Gene-set membership changes between
    MSigDB releases and feeds straight into every score, so the version is
    pinned in `msigdb.MSIGDB_VERSION` and each file's md5 is verified.
    """
    root = _resolve_data_dir(data_dir)
    raw_dir = root / "raw"
    typer.echo(f"raw dir: {raw_dir}")
    typer.echo(f"MSigDB:  {msigdb.MSIGDB_VERSION}")
    for key in collection or ["hallmark"]:
        path = msigdb.fetch_collection(key, raw_dir)
        sets = ssgsea.load_gmt(path)
        sizes = sorted(len(v) for v in sets.values())
        typer.echo(f"  {key:<10} {len(sets):>5} sets  sizes {sizes[0]}-{sizes[-1]}  -> {path.name}")


@app.command("build")
def build_cmd(
    project: ProjectFilterOpt = None,
    data_dir: DataDirOpt = None,
) -> None:
    """Flatten raw clinical JSON into per-project patient Parquets + dataset card."""
    root = _resolve_data_dir(data_dir)
    raw_dir = root / "raw"
    processed_dir = root / "processed_patient"
    typer.echo(f"raw dir:       {raw_dir}")
    typer.echo(f"processed dir: {processed_dir}")

    if not raw_dir.exists():
        raise typer.BadParameter(f"raw dir does not exist: {raw_dir}. Run fetch-clinical first.")

    project_files = _select_projects(raw_dir, project)
    _clear_processed(processed_dir, project)

    for path in project_files:
        proj = path.parent.name
        cases = json.loads(path.read_text())
        rows = clinical.to_patient_rows(cases)

        # Attach molecular modalities if their raw data has been fetched.
        mut_by_case = mutations.load_for_project(path.parent)
        expr_by_case = expression.load_for_project(path.parent)
        path_by_case = pathology.load_for_project(path.parent)
        ascn_by_case = copy_number.load_allele_specific_for_project(path.parent)
        mcn_by_case = copy_number.load_masked_for_project(path.parent)
        mirna_by_case = mirna.load_for_project(path.parent)
        rppa_by_case = protein_expression.load_for_project(path.parent)
        if mut_by_case:
            mutations.attach(rows, mut_by_case)
        if expr_by_case:
            expression.attach(rows, expr_by_case)
        if path_by_case:
            pathology.attach(rows, path_by_case)
        if ascn_by_case:
            copy_number.attach(rows, ascn_by_case, "samples_allele_specific_copy_number_segment")
        if mcn_by_case:
            copy_number.attach(rows, mcn_by_case, "samples_masked_copy_number_segment")
        if mirna_by_case:
            mirna.attach(rows, mirna_by_case)
        if rppa_by_case:
            protein_expression.attach(rows, rppa_by_case)
        # ssGSEA pathway activity, one column per MSigDB collection. The TPM
        # matrix is memoized across collections: reading every STAR TSV
        # dominates the cost and the matrix is identical for all of them.
        from tcga2hf.schema import SSGSEA_COLLECTIONS, ssgsea_patient_column

        msigdb_dir = raw_dir / "msigdb"

        def _memoized_matrix(project_dir: Path):
            """One TPM matrix per project, shared across collections.

            Built by a factory rather than a closure over a loop variable so
            the memo cannot outlive the project it belongs to.
            """
            cache: list = []

            def _load():
                if not cache:
                    cache.append(expression.tpm_matrix_for_project(project_dir))
                return cache[0]

            return _load

        _matrix = _memoized_matrix(path.parent)
        n_ssgsea = 0
        for _coll in SSGSEA_COLLECTIONS:
            by_case = ssgsea.load_for_project(path.parent, _coll, msigdb_dir, _matrix)
            if by_case:
                ssgsea.attach(rows, by_case, ssgsea_patient_column(_coll))
                n_ssgsea += sum(len(v) for v in by_case.values())
        # Attach BCR biotab Clinical Supplement data if it's been fetched.
        # survival.attach_survival reads `row["clinical_supplement"]` for the
        # ~2.4×-better-populated `treatment_outcome_first_course` field,
        # which drives DFI re-derivation. Supplement data is consumed
        # in-memory only; not serialized to the patients dataset.
        supp_dir = path.parent / "clinical_supplement"
        supps = clinical_supplement.load_supplements_for_project(supp_dir)
        if supps:
            clinical_supplement.attach_supplements(rows, supps)
        # BCR biotab Biospecimen Supplements — the specimen chain (slide
        # tumour-nuclei percentages, analyte quality, plate/shipment
        # provenance, site-specific factors). Flex-schema like the clinical
        # supplements; serialized as its own per-project inferred struct.
        bio_supps = biospecimen_supplement.load_supplements_for_project(
            path.parent / "biospecimen_supplement"
        )
        if bio_supps:
            biospecimen_supplement.attach_supplements(rows, bio_supps)
        # Re-derive OS / DSS / PFI / DFI from current GDC data using Liu
        # et al. 2018's algorithm; results land in `survival_derived` struct.
        # See `dev_research/liu_2018/report.html` for validation against
        # Liu's curated 2018 CDR.
        survival.attach_survival(rows)

        n_variants = sum(len(r["samples_masked_somatic_mutation"]) for r in rows)
        n_expr = sum(len(r["samples_gene_expression_quantification"]) for r in rows)
        n_path = sum(len(r["samples_pathology_report"]) for r in rows)
        n_ascn = sum(len(r["samples_allele_specific_copy_number_segment"]) for r in rows)
        n_mcn = sum(len(r["samples_masked_copy_number_segment"]) for r in rows)
        n_mirna = sum(len(r["samples_mirna_expression_quantification"]) for r in rows)
        n_rppa = sum(len(r["samples_protein_expression_quantification"]) for r in rows)
        n_bio = sum(1 for r in rows if r.get("biospecimen_supplement"))
        n_os = sum(1 for r in rows if (r.get("survival_derived") or {}).get("os_event") is not None)
        typer.echo(
            f"  {proj:<12} {len(rows):>4} patients  "
            f"mutations={n_variants:>5} ({len(mut_by_case)} MAFs)  "
            f"expression={n_expr:>4} aliquots ({len(expr_by_case)} TSVs)  "
            f"path_reports={n_path:>4}  ssgsea={n_ssgsea:>4}  "
            f"cnv={n_ascn:>4}/{n_mcn:<4} mirna={n_mirna:>4} rppa={n_rppa:>4} "
            f"biospec={n_bio:>4}  os={n_os}/{len(rows)}"
        )

        out = clinical.write_patients(rows, processed_dir, proj)
        typer.echo(f"             -> {out}")

    projects, gdc_releases = _card_inputs(raw_dir, processed_dir, "data.parquet")
    card = dataset_card.write_card(processed_dir, projects, gdc_releases=gdc_releases)
    typer.echo(f"wrote dataset card -> {card}")
    typer.echo(f"done. processed_patient tree at: {processed_dir}")


@app.command("build-tabular")
def build_tabular_cmd(
    project: ProjectFilterOpt = None,
    table: Annotated[
        list[str] | None,
        typer.Option(
            "--table",
            help=(
                "Table name (repeatable). If omitted, every table is built. If "
                "given, only those tables are re-derived and written; the "
                "project's other parquets are left untouched. Use this to append "
                "a newly-added modality without re-deriving the whole tree."
            ),
        ),
    ] = None,
    data_dir: DataDirOpt = None,
) -> None:
    """Flatten raw TCGA data into per-(project, table) Parquets + dataset card.

    Companion to `build`: same raw inputs, different output shape — one HF
    subset per (project, table) under <data-dir>/processed_tabular/.
    """
    root = _resolve_data_dir(data_dir)
    raw_dir = root / "raw"
    processed_dir = root / "processed_tabular"
    typer.echo(f"raw dir:       {raw_dir}")
    typer.echo(f"processed dir: {processed_dir}")

    if not raw_dir.exists():
        raise typer.BadParameter(f"raw dir does not exist: {raw_dir}. Run fetch-clinical first.")

    project_files = _select_projects(raw_dir, project)

    only: set[str] | None = None
    if table:
        known = (
            set(TABULAR_TABLES)
            | {f"clinical_supplement_{kind}" for kind in clinical_supplement.TABULAR_FORM_KINDS}
            | {
                f"biospecimen_supplement_{kind}"
                for kind in biospecimen_supplement.TABULAR_FORM_KINDS
            }
        )
        unknown = set(table) - known
        if unknown:
            raise typer.BadParameter(f"unknown table(s): {sorted(unknown)}. Known: {sorted(known)}")
        only = set(table)
        typer.echo(f"tables:        {', '.join(sorted(only))} (other tables left as-is)")

    # Clear prior output so removed projects/tables don't leave stale data
    # behind. Scope the wipe to what we're about to rebuild: with --project
    # only those project dirs go, leaving every other project's parquets in
    # place so a single project can be re-derived (or a new table appended)
    # without re-running — or re-uploading — the whole cohort.
    # A --table build must not wipe the tables it isn't rebuilding, so the
    # project-level wipe is skipped entirely; write_tables overwrites each
    # named table's parquet in place.
    if only is None:
        _clear_processed(processed_dir, project)

    for path in project_files:
        proj = path.parent.name
        cases = json.loads(path.read_text())
        tables = tabular.build_tables(cases, path.parent, only=only)
        # Re-derive survival endpoints (Liu et al. 2018 algorithm) — same
        # in-memory enrichment the consolidated `build` does. Lands on
        # `tables["cases"][i]["survival_derived"]`; tabular.write_tables
        # then projects it into the standalone `survival_derived` table.
        # Skipped when neither table was requested — nothing would read it.
        os_note = ""
        if "cases" in tables:
            supp_dir = path.parent / "clinical_supplement"
            supps = clinical_supplement.load_supplements_for_project(supp_dir)
            if supps:
                clinical_supplement.attach_supplements(tables["cases"], supps)
            survival.attach_survival(tables["cases"])
            n_cases = len(tables["cases"])
            n_os = sum(
                1
                for r in tables["cases"]
                if (r.get("survival_derived") or {}).get("os_event") is not None
            )
            os_note = f"  os={n_os}/{n_cases}"
            if only is None or "survival_derived" in only:
                tables["survival_derived"] = tabular.derived_survival_rows(tables["cases"])
            elif "cases" not in only:
                # `cases` was only pulled in because survival_derived is
                # projected off it — drop it so a --table build doesn't
                # rewrite a table nobody asked for. It stays when the
                # caller named `cases` explicitly.
                tables.pop("cases")

        # Row counts come back from write_tables rather than being taken
        # before it: the per-gene tables are streamed as Arrow batches and
        # are never materialised, so their size isn't known until written.
        sizes: dict[str, int] = {}
        out_paths = tabular.write_tables(tables, processed_dir, proj, counts=sizes)
        # Compact one-line summary of row counts per table — easier to spot
        # cardinality regressions than scrolling through 14 lines per project.
        summary = " ".join(f"{name}={n}" for name, n in sizes.items())
        typer.echo(f"  {proj:<12} {summary}{os_note}")
        # Use any one table's path to print the project's output dir.
        # write_tables returns {} when every requested table was a
        # flex-schema one with no rows (nothing is written for those).
        if out_paths:
            typer.echo(f"             -> {next(iter(out_paths.values())).parent.parent}")

    # ssGSEA reference distributions are cohort-level, so they can only be
    # built once every project's scores are on disk. Reading them back from
    # the written parquets (rather than from memory) is what guarantees the
    # published stats describe the published scores.
    _write_ssgsea_stats(processed_dir, only)

    projects, gdc_releases = _card_inputs(raw_dir, processed_dir, "*/data.parquet")

    # Tables list combines the fixed-schema TABULAR_TABLES with the
    # flex-schema clinical_supplement_* tables (one per BCR biotab form);
    # _tabular_configs_yaml filters to (project, table) pairs that
    # actually exist on disk.
    all_tables = (
        list(TABULAR_TABLES)
        + [f"clinical_supplement_{kind}" for kind in clinical_supplement.TABULAR_FORM_KINDS]
        + [f"biospecimen_supplement_{kind}" for kind in biospecimen_supplement.TABULAR_FORM_KINDS]
    )
    card = dataset_card.write_tabular_card(
        processed_dir,
        projects,
        all_tables,
        gdc_releases=gdc_releases,
    )
    typer.echo(f"wrote dataset card -> {card}")
    typer.echo(f"done. processed_tabular tree at: {processed_dir}")


@app.command("build-project-tabular")
def build_project_tabular_cmd(
    project: Annotated[
        str,
        typer.Option("--project", help="TCGA project id, e.g. --project TCGA-BRCA."),
    ],
    data_dir: DataDirOpt = None,
    table: Annotated[
        list[str] | None,
        typer.Option(
            "--table",
            help="Build only these tables (repeatable). Others are left as-is on disk.",
        ),
    ] = None,
    msigdb_dir: Annotated[
        Path | None,
        typer.Option("--msigdb-dir", help="MSigDB GMT dir; defaults to <data-dir>/raw/msigdb."),
    ] = None,
) -> None:
    """Build one project's standalone tabular dataset.

    Writes `<data-dir>/processed_project_tabular/<PROJECT>/`, whose contents
    are the repo root for `gabrielaltay/tcga-<slug>-tabular-open`: one
    `<table>/data.parquet` per config, plus the card.

    This is the per-project counterpart to `build-tabular`. The reason it
    exists is the HF dataset viewer: the pan-cancer repo declares 1,188
    configs (33 projects x 36 tables) and its splits never get scheduled —
    every one sits `pending` and the viewer is off. One repo per project
    declares ~40 configs, well inside the range that works.

    The published tree is standalone — no table here joins against another
    dataset. The *build* does read the pan-cancer tree once, for ssGSEA
    reference distributions: those are cohort-level by construction (a
    pathway's pan-cancer median depends on every project), so they are
    computed from `processed_tabular/` and this project's rows copied in.
    Skipped with a warning when that tree isn't present.
    """
    root = _resolve_data_dir(data_dir)
    raw_dir = root / "raw"
    processed_dir = root / "processed_project_tabular"
    project_dir = processed_dir / project
    pancancer_dir = root / "processed_tabular"

    cases_path = raw_dir / project / "cases.json"
    if not cases_path.exists():
        raise typer.BadParameter(f"no cases.json for {project} under {raw_dir}.")

    only: set[str] | None = None
    if table:
        known = (
            set(TABULAR_TABLES)
            | {f"clinical_supplement_{kind}" for kind in clinical_supplement.TABULAR_FORM_KINDS}
            | {
                f"biospecimen_supplement_{kind}"
                for kind in biospecimen_supplement.TABULAR_FORM_KINDS
            }
        )
        unknown = set(table) - known
        if unknown:
            raise typer.BadParameter(f"unknown table(s): {sorted(unknown)}. Known: {sorted(known)}")
        only = set(table)
        typer.echo(f"tables:        {', '.join(sorted(only))} (other tables left as-is)")
    elif project_dir.exists():
        # Full rebuild: drop prior output so a table removed upstream doesn't
        # linger as an orphaned parquet the card would still declare.
        shutil.rmtree(project_dir)

    typer.echo(f"raw dir:       {raw_dir}")
    typer.echo(f"output:        {project_dir}")

    cases = json.loads(cases_path.read_text())
    tables = tabular.build_tables(
        cases,
        cases_path.parent,
        only=only,
        msigdb_dir=msigdb_dir or (raw_dir / "msigdb"),
    )

    os_note = ""
    if "cases" in tables:
        supps = clinical_supplement.load_supplements_for_project(
            cases_path.parent / "clinical_supplement"
        )
        if supps:
            clinical_supplement.attach_supplements(tables["cases"], supps)
        survival.attach_survival(tables["cases"])
        n_os = sum(
            1
            for r in tables["cases"]
            if (r.get("survival_derived") or {}).get("os_event") is not None
        )
        os_note = f"  os={n_os}/{len(tables['cases'])}"
        if only is None or "survival_derived" in only:
            tables["survival_derived"] = tabular.derived_survival_rows(tables["cases"])
        elif "cases" not in only:
            # `cases` was only pulled in because survival_derived is projected
            # off it — drop it so a --table build doesn't rewrite a table
            # nobody asked for. It stays when the caller named it explicitly.
            tables.pop("cases")

    # The ssGSEA stats tables are cohort-level and are written separately,
    # below, from the pan-cancer tree. `build_tables` returns them empty, and
    # writing that would leave a 0-row parquet the card then declares as a
    # config — worse than absent — if the pan-cancer tree turns out to be
    # missing. Drop them here so the only writer is the one with real rows.
    for name in [t for t in tables if t.startswith("ssgsea_stats_")]:
        tables.pop(name)

    sizes: dict[str, int] = {}
    tabular.write_tables(tables, processed_dir, project, counts=sizes)
    for name, n in sorted(sizes.items()):
        typer.echo(f"  {name:<44}{n:>12,} rows")
    typer.echo(f"  {'':<44}{'':>12}{os_note}")

    _write_project_ssgsea_stats(pancancer_dir, project_dir, project, only)

    all_tables = (
        list(TABULAR_TABLES)
        + [f"clinical_supplement_{kind}" for kind in clinical_supplement.TABULAR_FORM_KINDS]
        + [f"biospecimen_supplement_{kind}" for kind in biospecimen_supplement.TABULAR_FORM_KINDS]
    )
    status_path = raw_dir / project / "gdc_status.json"
    gdc_release = (
        json.loads(status_path.read_text()).get("data_release") if status_path.exists() else None
    )
    card = dataset_card.write_project_tabular_card(
        project_dir, project, all_tables, gdc_release=gdc_release
    )
    typer.echo(f"wrote dataset card -> {card}")
    slug = project.lower().replace("tcga-", "")
    typer.echo(f"done. upload with: --repo-id gabrielaltay/tcga-{slug}-tabular-open")


def _write_project_ssgsea_stats(
    pancancer_dir: Path,
    project_dir: Path,
    project: str,
    only: set[str] | None,
) -> None:
    """Copy this project's ssGSEA reference rows out of the pan-cancer tree.

    The stats carry both the project's own distributions and the pan-cancer
    ones, and the latter are only correct when computed over every project —
    hence reading `processed_tabular/` rather than the single-project tree.
    """
    import pyarrow as pa
    import pyarrow.parquet as pq
    from tcga2hf.schema import SSGSEA_COLLECTIONS, TABULAR_TABLES, ssgsea_stats_table

    if not pancancer_dir.exists():
        typer.echo(f"  ssgsea stats skipped: no pan-cancer tree at {pancancer_dir}")
        return
    for coll in SSGSEA_COLLECTIONS:
        table_name = ssgsea_stats_table(coll)
        if only is not None and table_name not in only:
            continue
        rows = tabular.ssgsea_stats_rows(pancancer_dir, coll).get(project)
        if not rows:
            continue
        out_path = project_dir / table_name / "data.parquet"
        out_path.parent.mkdir(parents=True, exist_ok=True)
        pq.write_table(
            pa.Table.from_pylist(rows, schema=TABULAR_TABLES[table_name]),
            out_path,
            compression="zstd",
            write_page_index=True,
        )
        typer.echo(f"  {table_name:<44}{len(rows):>12,} rows (pan-cancer reference)")


@app.command("verify-project")
def verify_project_cmd(
    project: Annotated[
        str,
        typer.Option("--project", help="TCGA project id, e.g. --project TCGA-CHOL."),
    ],
    data_dir: DataDirOpt = None,
    sample: Annotated[
        int,
        typer.Option("--sample", help="Raw files per modality to re-hash against GDC's md5."),
    ] = 3,
) -> None:
    """Check a built project dataset against the live GDC.

    The unit tests assert the pipeline does what its author intended; this
    asks whether the published tree agrees with what the GDC serves today.
    It hits the API and re-hashes bytes on disk rather than trusting any
    manifest we wrote.

    Exits non-zero if any check fails, so it can gate an upload.
    """
    root = _resolve_data_dir(data_dir)
    checks = verify.verify_project(
        project,
        root / "raw",
        root / "processed_project_tabular",
        sample=sample,
    )
    typer.echo(f"verifying {project} against {GDC_BASE_URL}\n")
    for check in checks:
        mark = "PASS" if check.passed else "FAIL"
        typer.echo(f"[{mark}] {check.name}: {check.summary}")
        for line in check.details:
            typer.echo(line)
    failed = [c.name for c in checks if not c.passed]
    typer.echo("")
    if failed:
        typer.echo(f"{len(failed)} check(s) failed: {', '.join(failed)}")
        raise typer.Exit(code=1)
    typer.echo(f"all {len(checks)} checks passed")


@app.command("upload")
def upload_cmd(
    repo_id: Annotated[
        str,
        typer.Option(
            "--repo-id",
            help="HF dataset repo id, e.g. galtay/tcga-patients-open.",
        ),
    ],
    private: Annotated[
        bool,
        typer.Option(
            "--private/--public",
            help="Upload as private (default) or public. Recommend private for first pushes.",
        ),
    ] = True,
    commit_message: Annotated[
        str | None,
        typer.Option(
            "--commit-message",
            "-m",
            help="Commit message for this upload.",
        ),
    ] = None,
    data_dir: DataDirOpt = None,
) -> None:
    """Push <data-dir>/processed_patient/ to a HuggingFace dataset repo."""
    root = _resolve_data_dir(data_dir)
    processed_dir = root / "processed_patient"
    visibility = "private" if private else "PUBLIC"
    typer.echo(f"processed dir: {processed_dir}")
    typer.echo(f"repo_id:       {repo_id} ({visibility})")
    if not private:
        typer.confirm(
            f"This will publish {repo_id} as PUBLIC. Continue?",
            abort=True,
        )

    url = hf_upload.upload_dataset(
        processed_dir=processed_dir,
        repo_id=repo_id,
        private=private,
        commit_message=commit_message,
    )
    typer.echo(f"\nuploaded -> {url}")


@app.command("upload-tabular")
def upload_tabular_cmd(
    repo_id: Annotated[
        str,
        typer.Option(
            "--repo-id",
            help="HF dataset repo id, e.g. gabrielaltay/tcga-tabular-open.",
        ),
    ],
    private: Annotated[
        bool,
        typer.Option(
            "--private/--public",
            help="Upload as private (default) or public.",
        ),
    ] = True,
    commit_message: Annotated[
        str | None,
        typer.Option("--commit-message", "-m", help="Commit message for this upload."),
    ] = None,
    data_dir: DataDirOpt = None,
) -> None:
    """Push <data-dir>/processed_tabular/ to a HuggingFace dataset repo.

    Companion to `upload`: same shape, just points at the tabular tree.
    """
    root = _resolve_data_dir(data_dir)
    processed_dir = root / "processed_tabular"
    visibility = "private" if private else "PUBLIC"
    typer.echo(f"processed dir: {processed_dir}")
    typer.echo(f"repo_id:       {repo_id} ({visibility})")
    if not private:
        typer.confirm(
            f"This will publish {repo_id} as PUBLIC. Continue?",
            abort=True,
        )

    url = hf_upload.upload_dataset(
        processed_dir=processed_dir,
        repo_id=repo_id,
        private=private,
        commit_message=commit_message,
    )
    typer.echo(f"\nuploaded -> {url}")


@app.command("upload-project-tabular")
def upload_project_tabular_cmd(
    project: Annotated[
        str,
        typer.Option("--project", help="TCGA project id, e.g. --project TCGA-BRCA."),
    ],
    repo_id: Annotated[
        str | None,
        typer.Option(
            "--repo-id",
            help="Defaults to gabrielaltay/tcga-<slug>-tabular-open for the project.",
        ),
    ] = None,
    private: Annotated[
        bool,
        typer.Option(
            "--private/--public",
            help="Upload as private, or public (the default for these datasets).",
        ),
    ] = False,
    commit_message: Annotated[
        str | None,
        typer.Option("--commit-message", "-m", help="Commit message for this upload."),
    ] = None,
    skip_verify: Annotated[
        bool,
        typer.Option(
            "--skip-verify",
            help="Upload without running verify-project first. Rarely what you want.",
        ),
    ] = False,
    data_dir: DataDirOpt = None,
) -> None:
    """Verify against the GDC, then push one project's tabular dataset.

    Uploads `<data-dir>/processed_project_tabular/<PROJECT>/`, which is the
    repo root: `README.md` plus one `<table>/data.parquet` per config.

    **Every publish costs a full dataset-viewer re-index.** The Hub
    invalidates all of the repo's cached splits and re-queues them, and the
    viewer is dark until that finishes — tens of minutes for a project with
    ~40 configs. So an upload is a deliberate act, not something to do after
    each edit: build and verify as often as you like, and push once you are
    ready to wait for the rebuild.

    That is also why verification runs *first* rather than after. Finding a
    problem post-upload means a second push and a second re-index; finding
    it here costs about ten seconds. `--skip-verify` exists for the case
    where the GDC API is unreachable.

    Uploads are **public by default**. These are open-access TCGA data with
    no redistribution restriction, and a private dataset is treated as low
    priority by the Hub's viewer queue — the thing that makes the repo
    usable. `--private` is there for a staging push.

    `upload_folder` publishes everything under the directory, so this also
    refuses to run while stray non-dataset files sit in it — a
    `.pytest_cache` has reached a public repo this way before.
    """
    root = _resolve_data_dir(data_dir)
    processed_dir = root / "processed_project_tabular" / project
    slug = project.lower().replace("tcga-", "")
    repo_id = repo_id or f"gabrielaltay/tcga-{slug}-tabular-open"
    visibility = "private" if private else "PUBLIC"
    if not processed_dir.exists():
        raise typer.BadParameter(
            f"{processed_dir} does not exist. Run `build-project-tabular --project {project}`."
        )
    typer.echo(f"processed dir: {processed_dir}")
    typer.echo(f"repo_id:       {repo_id} ({visibility})")

    strays = [
        p.relative_to(processed_dir)
        for p in processed_dir.rglob("*")
        if p.is_file() and p.name not in {"README.md", "data.parquet"}
    ]
    if strays:
        listed = ", ".join(str(s) for s in sorted(strays)[:10])
        raise typer.BadParameter(
            f"{len(strays)} unexpected file(s) under {processed_dir} would be "
            f"published: {listed}. Remove them and re-run."
        )

    parquets = sorted(processed_dir.glob("*/data.parquet"))
    total = sum(q.stat().st_size for q in parquets)
    typer.echo(f"{len(parquets)} table(s), {total / 1e9:.2f} GB")

    if skip_verify:
        typer.echo("skipping verification (--skip-verify)")
    else:
        typer.echo("\nverifying against the GDC before publishing ...")
        checks = verify.verify_project(
            project, root / "raw", root / "processed_project_tabular"
        )
        for check in checks:
            typer.echo(f"  [{'PASS' if check.passed else 'FAIL'}] {check.name}: {check.summary}")
            if not check.passed:
                for line in check.details:
                    typer.echo(f"  {line}")
        failed = [c.name for c in checks if not c.passed]
        if failed:
            typer.echo("")
            raise typer.BadParameter(
                f"not uploading: {', '.join(failed)} failed. Fix and rebuild, or pass "
                "--skip-verify if the GDC API is unreachable."
            )
        typer.echo("")

    url = hf_upload.upload_dataset(
        processed_dir=processed_dir,
        repo_id=repo_id,
        private=private,
        commit_message=commit_message or f"Update {project} tabular (open access) dataset",
        # A table dropped locally must not keep being served from the repo.
        delete_patterns=["*/data.parquet"],
    )
    typer.echo(f"\nuploaded -> {url}")


@app.command("build-webdataset")
def build_webdataset_cmd(
    project: ProjectFilterOpt = None,
    data_dir: DataDirOpt = None,
    no_gzip: Annotated[
        bool,
        typer.Option(
            "--no-gzip",
            help=(
                "Store members as verbatim GDC bytes instead of gzipping them. "
                "Roughly triples the shard size; md5sums then verify in place "
                "without a decompression step."
            ),
        ),
    ] = False,
    shard_bytes: Annotated[
        int,
        typer.Option("--shard-bytes", help="Approximate target size per tar shard."),
    ] = wds_mod.SHARD_TARGET_BYTES,
) -> None:
    """Pack raw GDC files into per-patient WebDataset shards + Parquet indexes.

    One sample per patient, members named for GDC's own `data_type` and
    `file_id`. Unlike the other two builds this one re-derives nothing: the
    shard members are the bytes GDC serves, with the BCR biotab supplements
    the single documented exception (row subsets of a project-scoped form).

    With `--project`, only those projects' shard directories are replaced and
    the index is merged with the rows already on disk, so a project can be
    appended without repacking the cohort.
    """
    root = _resolve_data_dir(data_dir)
    raw_dir = root / "raw"
    processed_dir = root / "processed_webdataset"
    typer.echo(f"raw dir:       {raw_dir}")
    typer.echo(f"processed dir: {processed_dir}")

    if not raw_dir.exists():
        raise typer.BadParameter(f"raw dir does not exist: {raw_dir}. Run fetch-clinical first.")

    project_files = _select_projects(raw_dir, project)
    processed_dir.mkdir(parents=True, exist_ok=True)

    # Carry forward index rows for projects we're not rebuilding, so
    # `--project` stays an append rather than a truncation.
    keep = {p.parent.name for p in project_files}
    case_rows: list[dict] = []
    file_rows: list[dict] = []
    import pyarrow.parquet as pq

    for name, sink in (("cases.parquet", case_rows), ("files.parquet", file_rows)):
        existing = processed_dir / name
        if existing.exists():
            sink.extend(
                r for r in pq.read_table(existing).to_pylist() if r["project_id"] not in keep
            )

    for path in project_files:
        proj = path.parent.name
        out_dir = processed_dir / "data" / proj
        if out_dir.exists():
            shutil.rmtree(out_dir)
        cases = json.loads(path.read_text())
        rows, members = wds_mod.build_project(
            cases,
            path.parent,
            out_dir,
            proj,
            gzip_members=not no_gzip,
            shard_target_bytes=shard_bytes,
        )
        case_rows.extend(rows)
        file_rows.extend(members)
        shards = sorted(out_dir.glob("*.tar"))
        total = sum(s.stat().st_size for s in shards)
        n_files = sum(r["n_files"] for r in rows)
        n_supp = sum(1 for m in members if m["subset_of_gdc_file"])
        typer.echo(
            f"  {proj:<12} {len(rows):>4} patients  {n_files:>5} GDC files  "
            f"{n_supp:>5} supplement slices  "
            f"{len(shards)} shard(s)  {total / 1e9:.2f} GB"
        )

    out = wds_mod.write_cases_index(case_rows, processed_dir)
    typer.echo(f"wrote cases -> {out} ({len(case_rows)} patients)")
    files_out = wds_mod.write_files_index(file_rows, processed_dir)
    n_gdc = sum(1 for r in file_rows if not r["subset_of_gdc_file"])
    typer.echo(
        f"wrote files -> {files_out} ({len(file_rows)} members: "
        f"{n_gdc} GDC files + {len(file_rows) - n_gdc} supplement slices)"
    )
    projects = sorted({r["project_id"] for r in case_rows})
    card = dataset_card.write_webdataset_card(processed_dir, projects, case_rows)
    typer.echo(f"wrote dataset card -> {card}")
    typer.echo(f"done. processed_webdataset tree at: {processed_dir}")


@app.command("upload-webdataset")
def upload_webdataset_cmd(
    repo_id: Annotated[
        str,
        typer.Option(
            "--repo-id",
            help="HF dataset repo id, e.g. gabrielaltay/tcga-wds-open.",
        ),
    ],
    private: Annotated[
        bool,
        typer.Option(
            "--private/--public",
            help="Upload as private (default) or public.",
        ),
    ] = True,
    commit_message: Annotated[
        str | None,
        typer.Option("--commit-message", "-m", help="Commit message for this upload."),
    ] = None,
    data_dir: DataDirOpt = None,
) -> None:
    """Push <data-dir>/processed_webdataset/ to a HuggingFace dataset repo.

    Companion to `upload` / `upload-tabular`, pointed at the shard tree.
    `upload_folder` publishes everything under the directory, so this refuses
    to run while stray non-dataset files are sitting in it.
    """
    root = _resolve_data_dir(data_dir)
    processed_dir = root / "processed_webdataset"
    visibility = "private" if private else "PUBLIC"
    typer.echo(f"processed dir: {processed_dir}")
    typer.echo(f"repo_id:       {repo_id} ({visibility})")

    expected = {"README.md", "cases.parquet", "files.parquet"}
    strays = [
        p.relative_to(processed_dir)
        for p in processed_dir.rglob("*")
        if p.is_file()
        and p.name not in expected
        and not (p.suffix == ".tar" and p.parent.parent.name == "data")
    ]
    if strays:
        listed = ", ".join(str(s) for s in sorted(strays)[:10])
        raise typer.BadParameter(
            f"{len(strays)} unexpected file(s) under {processed_dir} would be "
            f"published: {listed}. Remove them and re-run."
        )

    shards = sorted(processed_dir.glob("data/*/*.tar"))
    total = sum(s.stat().st_size for s in shards)
    projects = sorted({s.parent.name for s in shards})
    typer.echo(f"{len(shards)} shard(s) across {len(projects)} project(s), {total / 1e9:.2f} GB")
    if not private:
        typer.confirm(
            f"This will publish {repo_id} as PUBLIC. Continue?",
            abort=True,
        )

    url = hf_upload.upload_dataset(
        processed_dir=processed_dir,
        repo_id=repo_id,
        private=private,
        commit_message=commit_message or "Update TCGA WebDataset (open access) dataset",
        # Root Parquets have been renamed before (index.parquet -> cases.parquet);
        # sweep any that are no longer part of the local tree so the repo never
        # serves a stale, undeclared table. Shards and .gitattributes are untouched.
        delete_patterns=["*.parquet"],
    )
    typer.echo(f"\nuploaded -> {url}")


@app.command("build-expression")
def build_expression_cmd(
    data_dir: DataDirOpt = None,
    project: Annotated[
        list[str] | None,
        typer.Option("--project", help="Limit to these projects (repeatable). Default: all built."),
    ] = None,
) -> None:
    """Build the cohort-wide expression matrix dataset.

    Writes `<data-dir>/processed_expression/`, the repo root for
    `gabrielaltay/tcga-expression-open`: one row per sample, the gene axis
    positional, one config per GDC quantification.

    Reads the per-project `gene_expression_quantification` tables rather
    than the raw TSVs, so a project must be built before it can be included
    — which also means every value here is one already published in a
    per-project dataset, reshaped and not recomputed.

    This is a reshape of a single modality, not a replacement for the
    per-project datasets: nothing clinical, and no other molecular data,
    lands here.
    """
    root = _resolve_data_dir(data_dir)
    processed_project_dir = root / "processed_project_tabular"
    out_dir = root / "processed_expression"
    typer.echo(f"source:  {processed_project_dir}")
    typer.echo(f"output:  {out_dir}")

    if out_dir.exists():
        # A project dropped from the selection must not linger as a stale
        # row block inside an otherwise-rebuilt matrix.
        shutil.rmtree(out_dir)

    counts, strand = expression_dataset.build(
        processed_project_dir,
        root / "raw",
        out_dir,
        projects=project,
        echo=typer.echo,
    )
    typer.echo("")
    for name, n in sorted(counts.items()):
        size = (out_dir / name / "data.parquet").stat().st_size
        typer.echo(f"  {name:<24}{n:>10,} rows{size / 1e6:>10,.0f} MB")

    projects_built = sorted(
        {p.name for p in processed_project_dir.glob("TCGA-*")} if project is None else set(project)
    )
    status_path = root / "raw" / projects_built[0] / "gdc_status.json"
    gdc_release = (
        json.loads(status_path.read_text()).get("data_release") if status_path.exists() else None
    )
    typer.echo(
        f"\nstrand balance first/(first+second): median {strand['balance_median']:.4f} "
        f"[{strand['balance_p01']:.3f}, {strand['balance_p99']:.3f}]; "
        f"{strand['n_strand_specific']} of {strand['n_samples']:,} samples outside 0.4-0.6"
    )
    card = dataset_card.write_expression_card(
        out_dir, counts, projects_built, strand, gdc_release=gdc_release
    )
    typer.echo(f"\nwrote dataset card -> {card}")
    typer.echo("upload with: tcga2hf-pipeline upload-expression")


@app.command("upload-expression")
def upload_expression_cmd(
    repo_id: Annotated[
        str,
        typer.Option("--repo-id", help="HF dataset repo id."),
    ] = "gabrielaltay/tcga-expression-open",
    private: Annotated[
        bool,
        typer.Option("--private/--public", help="Upload as private, or public."),
    ] = False,
    commit_message: Annotated[
        str | None,
        typer.Option("--commit-message", "-m", help="Commit message for this upload."),
    ] = None,
    data_dir: DataDirOpt = None,
) -> None:
    """Push `<data-dir>/processed_expression/` to HF Hub.

    Same discipline as `upload-project-tabular`: every push costs a full
    dataset-viewer re-index, so build as often as you like and upload once.

    Refuses to run while stray non-dataset files sit in the tree, since
    `upload_folder` would publish them.
    """
    root = _resolve_data_dir(data_dir)
    processed_dir = root / "processed_expression"
    if not processed_dir.exists():
        raise typer.BadParameter(f"{processed_dir} does not exist. Run `build-expression`.")

    strays = [
        p.relative_to(processed_dir)
        for p in processed_dir.rglob("*")
        if p.is_file() and p.name not in {"README.md", "data.parquet"}
    ]
    if strays:
        listed = ", ".join(str(s) for s in sorted(strays)[:10])
        raise typer.BadParameter(
            f"{len(strays)} unexpected file(s) under {processed_dir} would be "
            f"published: {listed}. Remove them and re-run."
        )

    parquets = sorted(processed_dir.glob("*/data.parquet"))
    total = sum(q.stat().st_size for q in parquets)
    visibility = "private" if private else "PUBLIC"
    typer.echo(f"processed dir: {processed_dir}")
    typer.echo(f"repo_id:       {repo_id} ({visibility})")
    typer.echo(f"{len(parquets)} config(s), {total / 1e9:.2f} GB")

    url = hf_upload.upload_dataset(
        processed_dir=processed_dir,
        repo_id=repo_id,
        private=private,
        commit_message=commit_message or "Update TCGA expression (open access) dataset",
        delete_patterns=["*/data.parquet"],
    )
    typer.echo(f"\nuploaded -> {url}")


if __name__ == "__main__":
    app()
