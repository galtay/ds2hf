## Section 10 — Extension beyond Liu (post-2018 coverage)

Liu's CDR is frozen at the 2018 data release. Our re-derivation runs against the modern GDC and includes patients who entered TCGA after Liu's freeze.

### Headline numbers

| | count |
|---|---|
| Total cohort (all 33 TCGA projects) | 11428 |
| In Liu's CDR (matched) | 11160 |
| **Post-Liu's freeze (extension)** | **268** |

### Re-derived endpoint coverage on the extension cohort

For the 268 post-freeze patients, how many have a populated value for each re-derived endpoint:

| endpoint | post-freeze patients with re-derived value | of total post-freeze |
| --- | --- | --- |
| OS | 116 | 268 |
| DSS | 116 | 268 |
| PFI | 116 | 268 |
| DFI | 73 | 268 |

### Top projects by extension count

Most of the +268 patients land in a few projects (TGCT and LUAD lead by a wide margin). The rest of TCGA was essentially closed by 2018; only a handful of projects continued accruing.

| project | post-freeze patients |
| --- | --- |
| TGCT | 129 |
| LUAD | 63 |
| GBM | 21 |
| OV | 21 |
| UCEC | 12 |
| DLBC | 10 |
| CHOL | 6 |
| COAD | 2 |
| READ | 2 |
| BRCA | 1 |
| LGG | 1 |

### Why TGCT / LUAD dominate the extension

Two structural reasons:

- **TGCT** had a small published cohort in Liu's 2018 freeze (134 patients) but BCR continued accruing testicular germ cell cases substantially through 2019-2020. The 129 new TGCT patients more than double the 2018 cohort.
- **LUAD** had ongoing tissue acquisitions through 2019 from the original TCGA centers; +63 new patients reflects late-stage data freezes.

Other projects (GBM, OV, UCEC, DLBC) added handfuls of new patients, mostly cases that were in-flight at the 2018 freeze and finalized later.

### Implication for downstream analysis

If you're reproducing Liu's results: use `gabrielaltay/tcga-pancanatlas-cdr`. You'll work on the 11,160 patients of the 2018 freeze, with Liu's values verbatim.

If you want maximum cohort size or current vital status: use `{os,dss,pfi,dfi}_event/_time` on all 11,428 patients. The +268 post-freeze patients are most useful for TGCT-specific analyses (where they're a meaningful fraction of the cohort). For most other projects the addition is small enough that conclusions wouldn't change either way.
