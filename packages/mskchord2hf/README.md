# mskchord2hf

Local loader for the [`gabrielaltay/msk-chord-2024`][hf] HuggingFace
dataset, a format-only mirror of the MSK-CHORD cBioPortal study. It builds
the views the mirror cannot ship: one nested row per patient, and the copy
number matrix as long rows.

```bash
pip install "mskchord2hf @ git+https://github.com/galtay/tcga2hf#subdirectory=packages/mskchord2hf"
```

```python
from mskchord2hf import Chord, cna_long

chord = Chord.from_hub()                 # downloads the parquets, not the source tarball
chord.table("mutations")                 # any config, as published (pyarrow.Table)
chord.dictionary("clinical_patient")     # display name, description, type, priority per column

patients = chord.patient_table()         # one row per patient, nested
p = chord.patient("P-0000012")           # the same row as a dict
p["samples"][0]["mutations"]             # that sample's mutation rows
p["timeline_treatment"]                  # treatment events, sorted by START_DATE

cna = cna_long(chord.table("cna"))       # (Hugo_Symbol, SAMPLE_ID, cna), unprofiled cells dropped
```

## Keep the results to yourself

MSK-CHORD is licensed [CC BY-NC-ND 4.0][license]. The license lets you make
restructured versions like these for non-commercial use (section
2(a)(1)(B)) but not share them. Share this code, or the unmodified mirror.

[hf]: https://huggingface.co/datasets/gabrielaltay/msk-chord-2024
[license]: https://creativecommons.org/licenses/by-nc-nd/4.0/
