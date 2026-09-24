"""The dataset card (README.md) for the MSK-CHORD mirror.

The card is also where the license's attribution terms are met (section
3(a)): it names the creator, links the license and the source, notes the
warranty disclaimer, and lists every change made to the files.
"""

from __future__ import annotations

from datetime import UTC, datetime
from pathlib import Path

from mskchord2hf_pipeline import source
from mskchord2hf_pipeline.build import ConfigInfo

REPO_ID = "gabrielaltay/msk-chord-2024"
GITHUB = "https://github.com/galtay/tcga2hf"
LOADER = f"{GITHUB}/tree/main/packages/mskchord2hf"

# Parquet files shipped but not declared as configs. The CNA matrix has one
# column per sample; `datasets` spends many minutes inferring features for a
# table that wide, and the HF viewer builds every declared config. Undeclared,
# the file is still in the repo for pyarrow, DuckDB and the loader. Making it
# narrow would mean reshaping it, which the license does not allow us to share.
FILE_ONLY = frozenset({"cna"})

# What one row of each config is.
_ROW_IS = {
    "clinical_patient": "a patient",
    "clinical_sample": "a sequenced tumor sample",
    "gene_panel_matrix": "a sample: the gene panel behind each data type",
    "mutations": "a somatic mutation in a sample",
    "sv": "a structural variant in a sample",
    "cna": "a gene; one column per sample",
    "cna_hg19_seg": "a copy-number segment in a sample (hg19)",
}


def _row_is(config: str) -> str:
    if config.startswith("timeline_"):
        return "an event in a patient's timeline"
    return _ROW_IS[config]


def _meta(members: dict[str, bytes], name: str) -> dict[str, str]:
    """A cBioPortal `meta_*.txt` file as its `key: value` pairs."""
    pairs = (line.split(":", 1) for line in members[name].decode().splitlines() if ":" in line)
    return {k.strip(): v.strip() for k, v in pairs}


def _cell(text: str) -> str:
    return " ".join(text.replace("|", "\\|").split())


def _dictionary_table(path: Path) -> str:
    import pyarrow.parquet as pq

    rows = ["| column | display name | type | description |", "|---|---|---|---|"]
    for f in pq.read_schema(path):
        m = {k.decode(): v.decode() for k, v in (f.metadata or {}).items()}
        rows.append(
            f"| `{f.name}` | {_cell(m['display_name'])} | {m['datatype']} | {_cell(m['description'])} |"
        )
    return "\n".join(rows)


def write_card(out_dir: Path, infos: list[ConfigInfo]) -> Path:
    members = source.read_members(out_dir / "source" / source.TARBALL_NAME)
    study = _meta(members, "meta_study.txt")
    cna_values = _meta(members, "meta_cna.txt")["profile_description"]
    rows = {i.name: i.rows for i in infos}
    timestamp = datetime.now(UTC).strftime("%Y-%m-%d %H:%M:%S UTC")

    configs_yaml = "\n".join(
        f"  - config_name: {i.name}\n    data_files:\n      - split: train\n        path: {i.name}/data.parquet"
        for i in infos
        if i.name not in FILE_ONLY
    )
    frontmatter = f"""---
license: cc-by-nc-nd-4.0
pretty_name: MSK-CHORD (MSK, Nature 2024)
tags:
  - cancer
  - clinical
  - genomics
  - real-world-data
  - cbioportal
configs:
{configs_yaml}
---
"""
    config_table = "\n".join(
        [
            "| config | source file | rows | one row is |",
            "|---|---|---:|---|",
            *(
                f"| `{i.name}` | `{i.member}` | {i.rows:,} | {_row_is(i.name)}{' (file only)' if i.name in FILE_ONLY else ''} |"
                for i in infos
            ),
        ]
    )

    body = f"""\
# MSK-CHORD

MSK-CHORD, the clinicogenomic dataset of [Jee et al. (2024)][paper], as Memorial Sloan Kettering Cancer Center (MSK) published it through [cBioPortal][study]. It pairs MSK-IMPACT tumor sequencing with clinical timelines, some derived by natural language processing (NLP).

- **Size:** {rows["clinical_patient"]:,} patients, {rows["clinical_sample"]:,} samples
- **Source:** [`{source.TARBALL_NAME}`][tarball], sha256 `{source.TARBALL_SHA256}`; byte-identical to [`{source.DATAHUB_REPO}`][datahub] at commit `{source.DATAHUB_COMMIT[:7]}`
- **Built:** {timestamp}

## License

**Non-commercial use only.** MSK licenses MSK-CHORD under [Creative Commons BY-NC-ND 4.0][license]. For commercial use, contact MSK at datarequests@mskcc.org.

**Share it only unmodified.** The license forbids sharing modified versions, including joined, reshaped or annotated copies. It does let you make them for your own non-commercial use.

**This repository changes only the file format.** The license's section 2(a)(4) permits format changes and states they do not create a modified version. Every value is as MSK published it; [Differences from the source files](#differences-from-the-source-files) lists each change.

**Neither MSK nor cBioPortal made or endorses this repository.**

## Configs

A *config* is one table, loaded by name. Each is one data file of the cBioPortal study, same rows, same columns, same order:

{config_table}

**`cna` is a file, not a config.** It has one column per sample, too wide for `load_dataset` and the dataset viewer. Read `cna/data.parquet` with pyarrow or DuckDB, or use the loader's `cna_long`; see [Copy number](#copy-number).

The study's other files, its `meta_*.txt` descriptions and `case_lists/`, are in `source/{source.TARBALL_NAME}`.

## Load it

```python
from datasets import load_dataset

REPO = "{REPO_ID}"
patients = load_dataset(REPO, "clinical_patient", split="train").to_pandas()
treatment = load_dataset(REPO, "timeline_treatment", split="train").to_pandas()
```

**Query without downloading everything.** Every config is one parquet file:

```python
import duckdb

duckdb.sql(\"\"\"
    SELECT CANCER_TYPE, count(*) AS samples
    FROM 'hf://datasets/{REPO_ID}/clinical_sample/data.parquet'
    GROUP BY CANCER_TYPE ORDER BY samples DESC
\"\"\")
```

## One row per patient, on your machine

**The [`mskchord2hf`][loader] package nests every config under its patient locally.** Sharing that nested table would share a modified version, so this repository ships the code to build it, not the table:

```python
# pip install "mskchord2hf @ git+{GITHUB}#subdirectory=packages/mskchord2hf"
from mskchord2hf import Chord

chord = Chord.from_hub()
p = chord.patient("P-0000012")
for sample in p["samples"]:
    print(sample["SAMPLE_ID"], sample["CANCER_TYPE"], len(sample["mutations"]))
for event in p["timeline_treatment"]:
    print(event["START_DATE"], event["AGENT"])
```

**Keep what it returns to yourself.** The license lets you make it for non-commercial use, not share it.

## Keys

**`PATIENT_ID` identifies a patient** (`P-0000012`). Every clinical and timeline config has it.

**`SAMPLE_ID` identifies a sequenced sample** (`P-0000012-T03-IM3`): the patient ID, a tumor number and the panel version. Sample-level configs name the column differently, as cBioPortal's formats do:

| config | sample column |
|---|---|
| `clinical_sample`, `gene_panel_matrix` | `SAMPLE_ID` |
| `mutations` | `Tumor_Sample_Barcode` |
| `sv` | `Sample_Id` |
| `cna_hg19_seg` | `ID` |
| `cna` | one column per sample, named by its `SAMPLE_ID` |

## Timelines

**Dates are day offsets, not calendar dates.** `START_DATE` and `STOP_DATE` count days from a per-patient day 0. `STOP_DATE` is null for events without a duration.

**Day 0 is the patient's first sequencing.** Each patient's earliest `Sequencing` event in `timeline_specimen` falls on day 0, so events before sequencing, such as the diagnosis, have negative days:

```python
spec = load_dataset(REPO, "timeline_specimen", split="train").to_pandas()
assert (spec[spec.EVENT_TYPE == "Sequencing"].groupby("PATIENT_ID").START_DATE.min() == 0).all()
```

**Some columns come from NLP.** Their clinical display names end in `(NLP)`; timeline probabilities such as `NLP_PROGRESSION_PROBABILITY` are the model's outputs. Jee et al. describe the models and their validation.

## Copy number

**`cna` is a gene by sample matrix, stored as text.** Read it directly:

```python
import pyarrow.parquet as pq
from huggingface_hub import hf_hub_download

cna = pq.read_table(hf_hub_download(REPO, "cna/data.parquet", repo_type="dataset"))
```

Its values, as the study's `meta_cna.txt` defines them: {cna_values}

**Count the values before relying on that definition.** The matrix need not hold every listed level, and can hold others:

```python
from mskchord2hf import cna_long  # see "One row per patient"

cna_long(cna).column("cna").value_counts()
```

**`NA` means the sample's panel did not cover the gene.** It is not a copy-number call. The loader's `cna_long` turns the matrix into (gene, sample, value) rows and drops these cells.

**`cna_hg19_seg` holds the segments behind the calls,** on the hg19 genome build.

## Gene panels

**A gene absent from a sample's panel was not tested.** `gene_panel_matrix` names the MSK-IMPACT panel behind each sample's mutations, copy number and structural variants. Treat a mutation missing from an uncovered gene as unknown, not wild type.

**The panels' gene lists are not in this repository.** cBioPortal keeps them in [datahub's `reference_data/gene_panels`][panels].

## Differences from the source files

**Only the format changed.** Each tab-separated file became one parquet file. The rules:

- **Line endings:** a line ends at a line feed, with or without a preceding carriage return. The lab timelines end lines with both.
- **Empty cells** are null. Every other token, including `NA`, `Unknown` and `None`, is kept as text, as are leading and trailing spaces.
- **Numbers:** a column is integer or float where every non-empty cell is a number that type holds exactly; otherwise text. Formatting is not kept: `0` and `0.0` load the same. The CNA matrix is text throughout.
- **Declared types win.** Clinical files declare each column's type; a column declared `STRING` stays text even when it holds numbers.
- **Header lines:** the four `#` lines that open each clinical file (display name, description, type, priority) are each column's parquet field metadata. The data dictionary below renders them.

Row order, column order and column names are the source's.

## Check against the source

The tarball is in this repository at `source/{source.TARBALL_NAME}`, and cBioPortal serves it without a login:

```python
import hashlib, io, tarfile, urllib.request
import pandas as pd

raw = urllib.request.urlopen("{source.TARBALL_URL}").read()
assert hashlib.sha256(raw).hexdigest() == "{source.TARBALL_SHA256}"

with tarfile.open(fileobj=io.BytesIO(raw)) as tar:
    text = tar.extractfile("{source.STUDY_ID}/data_timeline_treatment.txt").read()
theirs = pd.read_csv(io.BytesIO(text), sep="\\t", dtype=str, keep_default_na=False)
ours = load_dataset(REPO, "timeline_treatment", split="train").to_pandas()
assert theirs.AGENT.tolist() == ours.AGENT.fillna("").tolist()
```

**A newer study version fails the hash.** cBioPortal corrects studies in place under the same URL; this repository stays on the pinned tarball until it is rebuilt.

## Clinical data dictionary

The four header lines of each clinical file, verbatim.

**`clinical_patient`**

{_dictionary_table(out_dir / "clinical_patient" / "data.parquet")}

**`clinical_sample`**

{_dictionary_table(out_dir / "clinical_sample" / "data.parquet")}

## Citation

Cite the paper when you use this dataset:

> Jee J, Fong C, Pichotta K, et al. Automated real-world data integration improves cancer outcome prediction. *Nature*. 2024;636(8043):728-736. doi:10.1038/s41586-024-08167-5

## Attribution

- **Title:** {study["name"]}
- **Creator:** Memorial Sloan Kettering Cancer Center
- **Source:** [cBioPortal study `{source.STUDY_ID}`][study]
- **License:** [CC BY-NC-ND 4.0][license]; the full text is in `LICENSE`, copied from the study.
- **Warranty:** the licensor offers the material as-is and makes no warranties; see the license's section 5.
- **Changes:** format only, as [listed above](#differences-from-the-source-files).

The pipeline that builds this repository is at [{GITHUB}][repo].

[paper]: https://doi.org/10.1038/s41586-024-08167-5
[study]: {source.STUDY_URL}
[tarball]: {source.TARBALL_URL}
[datahub]: https://github.com/{source.DATAHUB_REPO}/tree/{source.DATAHUB_COMMIT}/{source.DATAHUB_PATH}
[panels]: https://github.com/{source.DATAHUB_REPO}/tree/master/reference_data/gene_panels
[license]: {source.LICENSE_URL}
[loader]: {LOADER}
[repo]: {GITHUB}
"""
    out_path = out_dir / "README.md"
    out_path.write_text(frontmatter + body)
    return out_path
