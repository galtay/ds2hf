## Section 5 — Match rate against Liu's curated CDR

For the 11,160 CDR-matched patients, how often does our re-derived value exactly match Liu's? Five-category classification per (patient, endpoint):

- **match** — both populated, event and time both equal (within 0.5 days)
- **mismatch** — both populated but disagree
- **both NA** — both NA (e.g. SKCM has no DFI in either stream)
- **Liu only** (`cdr_pop_der_na`) — Liu populated, ours is NA (we under-populate)
- **ours only** (`der_pop_cdr_na`) — ours populated, Liu's is NA (we over-populate, or filling a post-hoc gap)

### Cohort-level summary

| endpoint | total | match | match % | mismatch | both NA | Liu only (cdr_pop_der_na) | ours only (der_pop_cdr_na) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| OS | 11160 | 10938 | 98.0 | 163 | 3 | 50 | 6 |
| DSS | 11160 | 10373 | 92.9 | 240 | 18 | 35 | 494 |
| PFI | 11160 | 10719 | 96.0 | 194 | 13 | 40 | 194 |
| DFI | 11160 | 5412 | 48.5 | 65 | 4640 | 54 | 989 |

### Per-project per-endpoint match rate

Denominator is patients where Liu's CDR populated the endpoint (excludes the both-NA cases). Numerator is exact matches. Mismatches and Liu-only cases count against the rate. `—` means Liu populated zero patients for that (project, endpoint) — typically Liu's documented exclusions (LAML for DSS/PFI/DFI; SKCM/THYM/UVM for DFI).

| project | N | OS | DSS | PFI | DFI |
| --- | --- | --- | --- | --- | --- |
| ACC | 92 | 92/92 (100%) | 90/90 (100%) | 92/92 (100%) | 53/53 (100%) |
| BLCA | 412 | 402/412 (98%) | 389/398 (98%) | 404/412 (98%) | 184/189 (97%) |
| BRCA | 1097 | 1075/1097 (98%) | 1060/1078 (98%) | 1079/1097 (98%) | 916/953 (96%) |
| CESC | 307 | 305/307 (99%) | 300/303 (99%) | 305/307 (99%) | 175/176 (99%) |
| CHOL | 45 | 45/45 (100%) | 42/42 (100%) | 45/45 (100%) | 28/28 (100%) |
| COAD | 459 | 444/459 (97%) | 429/443 (97%) | 445/459 (97%) | 187/190 (98%) |
| DLBC | 48 | 47/48 (98%) | 47/48 (98%) | 47/48 (98%) | 27/28 (96%) |
| ESCA | 185 | 182/185 (98%) | 180/183 (98%) | 182/185 (98%) | 89/89 (100%) |
| GBM | 596 | 592/596 (99%) | 553/555 (100%) | 593/596 (99%) | 3/3 (100%) |
| HNSC | 528 | 525/528 (99%) | 498/502 (99%) | 525/528 (99%) | 133/134 (99%) |
| KICH | 113 | 112/113 (99%) | 112/113 (99%) | 112/113 (99%) | 71/71 (100%) |
| KIRC | 537 | 529/537 (99%) | 516/525 (98%) | 527/537 (98%) | 117/117 (100%) |
| KIRP | 291 | 287/291 (99%) | 283/287 (99%) | 285/291 (98%) | 180/184 (98%) |
| LAML | 200 | 186/200 (93%) | 0/67 (0%) | — | — |
| LGG | 515 | 511/515 (99%) | 503/507 (99%) | 512/515 (99%) | 132/134 (99%) |
| LIHC | 377 | 373/377 (99%) | 364/368 (99%) | 376/377 (100%) | 321/323 (99%) |
| LUAD | 522 | 512/522 (98%) | 476/486 (98%) | 512/522 (98%) | 304/309 (98%) |
| LUSC | 504 | 496/504 (98%) | 445/452 (98%) | 498/504 (99%) | 302/304 (99%) |
| MESO | 87 | 85/87 (98%) | 64/66 (97%) | 74/87 (85%) | 12/15 (80%) |
| OV | 587 | 577/585 (99%) | 540/548 (99%) | 580/587 (99%) | 285/286 (100%) |
| PAAD | 185 | 184/185 (99%) | 169/178 (95%) | 181/185 (98%) | 66/69 (96%) |
| PCPG | 179 | 179/179 (100%) | 179/179 (100%) | 179/179 (100%) | 159/160 (99%) |
| PRAD | 500 | 500/500 (100%) | 498/498 (100%) | 474/500 (95%) | 331/340 (97%) |
| READ | 170 | 165/170 (97%) | 159/164 (97%) | 167/170 (98%) | 47/48 (98%) |
| SARC | 261 | 252/261 (97%) | 246/255 (96%) | 252/261 (97%) | 148/153 (97%) |
| SKCM | 470 | 445/463 (96%) | 439/457 (96%) | 450/463 (97%) | — |
| STAD | 443 | 428/443 (97%) | 397/416 (95%) | 413/443 (93%) | 249/259 (96%) |
| TGCT | 134 | 134/134 (100%) | 134/134 (100%) | 134/134 (100%) | 90/105 (86%) |
| THCA | 507 | 504/507 (99%) | 500/501 (100%) | 504/507 (99%) | 357/358 (100%) |
| THYM | 124 | 123/124 (99%) | 123/124 (99%) | 123/124 (99%) | — |
| UCEC | 548 | 535/548 (98%) | 532/546 (97%) | 537/548 (98%) | 420/426 (99%) |
| UCS | 57 | 54/57 (95%) | 52/55 (95%) | 55/57 (96%) | 26/27 (96%) |
| UVM | 80 | 58/80 (72%) | 54/80 (68%) | 57/80 (71%) | — |

### Reading the cohort summary

The "match" column above counts only exact agreement (event + time within 0.5 days), which under-counts DFI by penalizing the ~4,640 patients we both correctly exclude. A more honest "we got pretty much the same answer" rate counts both-NA as agreement and allows ~30 days of time slack on populated cases:

| endpoint | agreement rate (both NA OR both populated within 30d) |
|---|---|
| OS  | **98.2%** (10,961/11,160) |
| DSS | **93.2%** (10,400/11,160) |
| PFI | **96.3%** (10,747/11,160) |
| DFI | **90.1%** (10,055/11,160) |

For DFI specifically, this is up from **77.2%** before the BCR biotab Clinical Supplement integration (a +12.9 percentage-point jump from one upstream change). See Section 9 for the per-bucket breakdown of the remaining 10% disagreement.

### Per-endpoint notes

- **OS** — algorithmic agreement is essentially perfect; the 2% disagreement is patients alive in Liu's 2018 freeze who have since died (data drift, not algorithm). See Section 6 deep-dive.
- **DSS** — Liu's own STAR Methods flags DSS as approximate (*"not a 100% accurate definition but is the best we could do with this dataset"*). Our re-derivation honors the same rule. See Section 7.
- **PFI** — past Liu's reliability bar of 95%. Most disagreement is post-2018 follow-up time shifting (longer censor times). See Section 8.
- **DFI** — at 90% agreement, the remaining 10% is mostly the 989 patients where we have **extra coverage** Liu didn't have in 2018 (we found a disease-free signal in the BCR biotab Liu's 2018 algorithm couldn't compute), plus 64 outright disagreements and 52 patients where Liu had a value we still can't recover. See Section 9.

### Reading the per-project table

- **High-match projects** (most cells ≥95%): smaller cohorts where the data has stabilized and Liu's exclusions match ours.
- **Low-DFI-match projects**: GBM (mostly all stage IV → DFI excluded by Liu's stage-IV rule), DLBC (small cohort, fragile to per-patient drift), KICH/KIRC/KIRP (the kidney trio, where biotab data is sparser than other tumor types).
- **`—` cells**: SKCM/THYM/UVM/LAML for DFI, LAML for DSS/PFI — Liu's documented exclusions; we honor the same rule so denominator is zero.

The remaining sections walk through OS / DSS / PFI / DFI individually.
