# mskchord2hf-pipeline

Build pipeline for the [`gabrielaltay/msk-chord-2024`][hf] HuggingFace
dataset: a mirror of the MSK-CHORD cBioPortal study (Jee et al., *Nature*
2024), converted from tab-separated text to parquet and otherwise unchanged.

For the per-patient view, install the companion [`mskchord2hf`](../mskchord2hf/)
loader instead.

## Why format conversion only

MSK licenses MSK-CHORD under [CC BY-NC-ND 4.0][license]. It may be shared,
for non-commercial use, only unmodified; section 2(a)(4) says format changes
do not count as modification. So the pipeline writes one parquet per study
data file, same rows, columns and order, and restructures nothing. Joins and
per-patient nesting live in the loader, which runs on the user's machine,
where the license (section 2(a)(1)(B)) allows them.

The conversion rules are in `convert.py`'s docstring, and the dataset card
lists them as the license's required notice of changes.

## Commands

```
mskchord2hf-pipeline fetch    # cBioPortal -> <data-dir>/raw/msk_chord_2024.tar.gz (sha256-pinned)
mskchord2hf-pipeline build    # tarball -> <data-dir>/processed/ (parquets, tarball, LICENSE, card)
mskchord2hf-pipeline verify   # every cell vs the tarball; the tarball vs cBioPortal and datahub
mskchord2hf-pipeline upload   # verify, then push processed/ (private unless --public)
```

Data lives under `$MSKCHORD2HF_DATA_DIR` (default `$HOME/data/mskchord2hf`).
`upload` needs a write-scoped `HF_TOKEN` in the repo-root `.env`.

## When the study changes

cBioPortal corrects studies in place under the same URL. `fetch` refuses
bytes other than the pinned tarball and keeps them beside it for review.
To move the pin, update `TARBALL_SHA256` and `DATAHUB_COMMIT` in
`source.py` to the new tarball and the datahub commit it matches, then
rebuild and verify.

[hf]: https://huggingface.co/datasets/gabrielaltay/msk-chord-2024
[license]: https://creativecommons.org/licenses/by-nc-nd/4.0/
