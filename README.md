# ds2hf

Public cancer datasets, published to the HuggingFace Hub with their
provenance intact. One module per source, in two packages:

- **`packages/ds2hf/`**: the **read side** ([README](packages/ds2hf/README.md)),
  for people using the published datasets. Small dependency set (`pyarrow`,
  `pydantic`, `huggingface-hub`).
  - `ds2hf.tcga`: typed pydantic models and pyarrow schemas for the TCGA datasets.
  - `ds2hf.mskchord`: a local loader that builds per-patient views of the
    MSK-CHORD mirror on your own machine.
- **`packages/ds2hf-pipeline/`**: the **build side**
  ([README](packages/ds2hf-pipeline/README.md)), for maintainers. The
  `ds2hf-pipeline <source> <command>` CLI fetches, converts, verifies,
  writes dataset cards and pushes to the Hub.

Sources share the upload helper and nothing else. Each has its own
license, keys and verification.

## Develop

This is a [uv](https://docs.astral.sh/uv/) workspace; one `uv sync` at
the repo root brings up both packages and dev tooling:

```bash
uv sync
uv run pytest -m "not network"     # unit + integration tests
uv run ds2hf-pipeline --help       # the build CLI
```

## Datasets

Every dataset has its own card, written by its source's build command.

**TCGA**, open-access data from the NCI Genomic Data Commons (GDC):

- [`gabrielaltay/tcga-patients-open`][patients]: consolidated, one row
  per patient with the GDC entity tree fully nested, plus molecular vectors.
- [`gabrielaltay/tcga-tabular-open`][tabular]: the same source data
  as flat tables per project.

**MSK-CHORD**, the Memorial Sloan Kettering clinicogenomic dataset
(Jee et al., *Nature* 2024) from cBioPortal:

- `gabrielaltay/msk-chord-2024`: the study's files, format conversion
  only, as its CC BY-NC-ND 4.0 license requires. Not yet published.

[patients]: https://huggingface.co/datasets/gabrielaltay/tcga-patients-open
[tabular]: https://huggingface.co/datasets/gabrielaltay/tcga-tabular-open
