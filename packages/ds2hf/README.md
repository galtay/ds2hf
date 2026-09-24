# ds2hf

Read-side companions for the ds2hf HuggingFace datasets, one module per
source. To produce or publish the datasets yourself, see
[`ds2hf-pipeline`](../ds2hf-pipeline/).

## Install

```bash
pip install "ds2hf @ git+https://github.com/galtay/ds2hf#subdirectory=packages/ds2hf"
```

Dependencies are limited to `pyarrow`, `pydantic` and `huggingface-hub`.

## TCGA: `ds2hf.tcga`

Typed pydantic models and pyarrow schemas for the public
[`gabrielaltay/tcga-patients-open`][patients] and
[`gabrielaltay/tcga-tabular-open`][tabular] datasets. Use them to validate
patient rows or browse them with full type information.

```python
import pyarrow.parquet as pq
from ds2hf.tcga import TcgaHfPatient

# Read one project's parquet from the consolidated dataset.
t = pq.read_table("TCGA-LUAD/data.parquet")
for row in t.to_pylist():
    patient = TcgaHfPatient.model_validate(row)
    for tumor, normal in patient.tumor_normal_pairs():
        print(tumor.submitter_id, "vs", normal.submitter_id)
    for variant in patient.mutations_by_gene().get("TP53", []):
        print(variant.HGVSp_Short, variant.t_alt_count, "/", variant.t_depth)
    for event in patient.timeline():
        print(f"day {event.day:+.0f}: {event.category} — {event.label}")
```

`TcgaHfPatient` mirrors the consolidated parquet schema field-for-field
(`extra="forbid"` strict mode) and adds convenience joins:

- `aliquot_to_sample` / `aliquot_lookup` — flatten the GDC biospecimen tree
- `tumor_normal_pairs` — distinct tumor/normal sample pairs from MAF rows
- `mutations_by_gene`, `mutations_by_consequence`
- `expression_for_gene` — per-aliquot expression lookup by gene symbol
- `timeline` — every dated event (consent, diagnosis, treatments,
  follow-ups, sample procurement, BCR receipt, death) sorted on the
  case's `index_date` anchor
- `consistency_check` — count of GDC quirks worth surfacing for
  downstream consumers (e.g. samples with `days_to_collection >
  days_to_death`)

The pyarrow `*_FIELDS` lists in `ds2hf.tcga.schema` are the single source of
truth for the dataset shape (regenerated from gdcdictionary YAMLs); the
pydantic models are derived from those same lists, so the two stay in
sync by construction.

## MSK-CHORD: `ds2hf.mskchord`

A local loader for the `gabrielaltay/msk-chord-2024` dataset, a
format-only mirror of the MSK-CHORD cBioPortal study. It builds the views
the mirror cannot ship: one nested row per patient, and the copy number
matrix as long rows.

```python
from ds2hf.mskchord import Chord, cna_long

chord = Chord.from_hub()                 # downloads the parquets, not the source tarball
chord.table("mutations")                 # any config, as published (pyarrow.Table)
chord.dictionary("clinical_patient")     # display name, description, type, priority per column

patients = chord.patient_table()         # one row per patient, nested
p = chord.patient("P-0000012")           # the same row as a dict
p["samples"][0]["mutations"]             # that sample's mutation rows
p["timeline_treatment"]                  # treatment events, sorted by START_DATE

cna = cna_long(chord.table("cna"))       # (Hugo_Symbol, SAMPLE_ID, cna), unprofiled cells dropped
```

**Keep the results to yourself.** MSK-CHORD is licensed
[CC BY-NC-ND 4.0][cc]. The license lets you make restructured versions
like these for non-commercial use (section 2(a)(1)(B)) but not share them.
Share this code, or the unmodified mirror.

## Source

This package ships from the [`galtay/ds2hf`][repo] monorepo; see the
repo for the build pipeline, dataset cards, and full documentation.

[patients]: https://huggingface.co/datasets/gabrielaltay/tcga-patients-open
[tabular]: https://huggingface.co/datasets/gabrielaltay/tcga-tabular-open
[repo]: https://github.com/galtay/ds2hf
[cc]: https://creativecommons.org/licenses/by-nc-nd/4.0/
