# tcga2hf-pipeline

Build pipeline for the
[`gabrielaltay/tcga-patients-open`][patients] and
[`gabrielaltay/tcga-tabular-open`][tabular] HuggingFace datasets. Provides
the `tcga2hf-pipeline` CLI that fetches public, open-access TCGA data from
the National Cancer Institute (NCI) Genomic Data Commons (GDC), flattens
it into per-project Parquets, and pushes the result to HF Hub.

For the read-side typed models / pyarrow schemas alone, install the
companion [`tcga2hf`](../tcga2hf/) package instead — it has only `pydantic`
+ `pyarrow` as dependencies.

## Install

```bash
uv sync   # at the repo root, inside the workspace
```

(Or `pip install tcga2hf-pipeline` once published.)

## Auth

Only the `upload` and `upload-tabular` commands need auth. Put a
write-scoped HF token (from https://huggingface.co/settings/tokens) into a
local `.env` file in the repo root — gitignored, auto-loaded by the CLI:

```
HF_TOKEN=hf_xxxxxxxxxxxxxxxxxxxx
```

The CLI calls `load_dotenv(override=True)` so this project's `.env` always
wins over any inherited shell `HF_TOKEN`.

## Data location

All data lives outside the repo, under `$TCGA2HF_DATA_DIR` (default
`$HOME/data/tcga2hf`), overridable per command via `--data-dir`. The GDC
fetch step is open-access and needs no auth.

## Commands

These hit the GDC for open-access data. Acronyms used below: MAF =
Mutation Annotation Format; WXS = Whole Exome Sequencing.

```
tcga2hf-pipeline fetch-clinical           # GDC -> raw cases.json per project
tcga2hf-pipeline fetch-mutations          # GDC -> Masked Somatic Mutation MAFs from WXS (DNA)
tcga2hf-pipeline fetch-expression         # GDC -> Gene Expression Quantification TSVs from RNA-Seq (STAR counts)
tcga2hf-pipeline fetch-pathology-reports  # GDC -> Pathology Report PDFs (scanned BCR documents)
tcga2hf-pipeline fetch-msigdb             # Broad -> md5-pinned MSigDB gene sets (GMT) for ssGSEA
tcga2hf-pipeline fetch-copy-number        # GDC -> every open-access copy number file (4 data types, 6 workflows)
tcga2hf-pipeline fetch-methylation        # GDC -> SeSAMe methylation beta values (Illumina arrays)
tcga2hf-pipeline fetch-mirna-isoform      # GDC -> miRNA isoform-level quantification
tcga2hf-pipeline build                    # raw -> consolidated per-project parquets + dataset card
tcga2hf-pipeline build-tabular            # raw -> per-(project, table) parquets + dataset card
tcga2hf-pipeline build-project-tabular    # raw -> ONE project's standalone tabular dataset + card
tcga2hf-pipeline build-webdataset         # raw -> per-patient WebDataset tar shards + index.parquet
tcga2hf-pipeline upload                   # push processed/ to the consolidated HF dataset repo
tcga2hf-pipeline upload-tabular           # push processed_tabular/ to the tabular HF dataset repo
tcga2hf-pipeline upload-project-tabular   # push one project's tree to tcga-<slug>-tabular-open
```

### One dataset per project

`build-project-tabular` / `upload-project-tabular` are the current target;
`build-tabular` / `upload-tabular` build the older pan-cancer repo, which is
kept as-is but no longer extended.

The reason is the HF dataset viewer. `tcga-tabular-open` declares 33 projects
x 36 tables = **1,188 configs**, and although HF's documented ceiling is
3,000, none of those splits ever get scheduled — every one reports `pending`
and `is-valid` returns `viewer: false`. A per-project repo declares ~40
configs, comfortably inside the range `tcga-patients-open` (33 configs) is
fully green in. Config names lose the project prefix too, since the project
*is* the repo: `load_dataset("gabrielaltay/tcga-brca-tabular-open", "cases")`.

With those, a project dataset covers **every open-access molecular
`data_type` GDC serves** except two raw-instrument formats: `Slide Image`
(whole-slide `.svs`, which outweighs everything else by an order of
magnitude) and `Masked Intensities` (the `.idat` files the methylation betas
are computed from). The `files` table indexes whatever was actually fetched.

Each project dataset is standalone — no table joins against another dataset.
The build reads the pan-cancer tree once, for ssGSEA reference distributions
(a pathway's pan-cancer median depends on every project, so it cannot be
computed from one), and copies that project's rows in.

Two tables there are stored **normalized against `gene_model`** rather than
as verbatim copies of their source TSVs — `gene_expression_quantification`
and `gene_level_copy_number` carry only `gene_id`, and the GENCODE v36
metadata lives once in the `gene_model` config. Every GDC per-gene file
repeats the same 60,660-row model, which accounted for 51% of one project's
expression table. The source file is exactly reconstructible with a join on
`gene_id`; `gene_model` is replicated into each project repo (~1.25 MiB) so
nothing has to reach outside its own dataset.

`build-webdataset` is the odd one out: it re-derives nothing. Each patient
becomes one WebDataset sample whose members are the open-access files GDC
serves for them, named for GDC's own `data_type` and `file_id`. The single
exception is the BCR biotab supplements, which GDC ships per *project* — those
members are per-patient row subsets, flagged `subset_of_gdc_file` in the
sample's `files.jsonl`. Members are gzipped unless GDC already ships them
compressed (`--no-gzip` to store verbatim bytes at ~3x the size).

Only `index.parquet` is declared as a config in that dataset's card. HF
`datasets` resolves one builder module per repo, so declaring the `.tar`
shards too would make it read them as Parquet and fail; leaving them
undeclared lets the viewer render the index while `webdataset` streams the
shards by URL.

`tcga2hf-pipeline <cmd> --help` for arguments. Notable flags:

- `--project TCGA-LUAD` (repeatable) on every fetch command and on both
  build commands. Omit on the build commands to process every project
  under `<data-dir>/raw/` and rebuild the output tree from scratch.
  **With** `--project`, only those projects' output directories are
  replaced — everything else is left byte-identical, so a new modality can
  be appended without re-deriving (or re-uploading) the whole cohort.
- `--table <name>` (repeatable) on `build-tabular` /
  `build-project-tabular`: build only those tables, leaving the rest on
  disk untouched. Adding one modality doesn't mean re-reading every MAF.
- `--skip-gene-level` on `fetch-copy-number`: fetch only the four
  segment-level sources. Gene-level copy number is ~97% of the modality's
  bytes (10.6 GiB for TCGA-BRCA, ~109 GB pan-TCGA).
- `--max-files N` on `fetch-mutations` / `fetch-expression`: cap the
  total number of files on disk per project (cached + freshly downloaded).
  Set to 1 to sample one file per project, or 0 to populate the manifest
  without downloading any new bytes. The manifest always lists every
  discovered file; entries trimmed off get `_status="manifest_only"`.

## Provenance

The GDC API only ever serves the current data release, so *when* a
modality was fetched says nothing about whether its bytes changed. Two
records pin that down:

- **Per modality**: each `raw/<project>/<modality>/` directory keeps its
  own `gdc_status.json` with the release it was fetched at, so modalities
  added years apart each report their own. `raw/<project>/gdc_status.json`
  belongs to `fetch-clinical` and describes `cases.json` alone.
- **Per file**: manifests record `gdc_version`, `gdc_first_release`, and
  `gdc_superseded` from `POST /files/versions`. A file with
  `gdc_superseded=false` is byte-identical to what every release since
  `gdc_first_release` served — which is the guarantee that actually
  matters when modalities are fetched at different times. These surface on
  the tabular `files` table.

## Tests

Run from the workspace root:

```bash
uv run pytest -m "not network and not integration"  # offline
uv run pytest                                        # all (incl. network)
```

## Source

This package ships from the [`galtay/tcga2hf`][repo] monorepo. The
generated dataset cards live in `dataset_card.py`; build them via
`tcga2hf-pipeline build` / `build-tabular`.

[patients]: https://huggingface.co/datasets/gabrielaltay/tcga-patients-open
[tabular]: https://huggingface.co/datasets/gabrielaltay/tcga-tabular-open
[repo]: https://github.com/galtay/tcga2hf
