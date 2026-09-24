## Section 4 — Per-project survival landscape

Liu's Figure 1B-E shows Kaplan-Meier curves per project for each endpoint. We surface the same information as event counts plus survival probabilities at 1-year, 3-year, and 5-year timepoints — the standard clinical reporting milestones. The same patterns Liu's panels show visually appear here numerically: aggressive cancers (SKCM, OV, GBM, MESO) drop fast; indolent cancers (TGCT, PRAD, THCA, KICH) stay flat.

Computed against our re-derived `{os,dss,pfi,dfi}_event/_time` columns on the full current cohort (11,428 patients). `—` means insufficient follow-up to estimate the rate at that timepoint (curve hasn't reached it). Projects with zero populated patients for an endpoint (e.g. LAML / SKCM / THYM / UVM for DFI — Liu's exclusions) are omitted from that table.

### OS — Overall Survival

| project | N (populated) | events | 1-yr | 3-yr | 5-yr |
| --- | --- | --- | --- | --- | --- |
| ACC | 92 | 34 | 91% | 72% | 59% |
| BLCA | 412 | 182 | 78% | 48% | 40% |
| BRCA | 1097 | 152 | 98% | 90% | 81% |
| CESC | 307 | 72 | 92% | 72% | 63% |
| CHOL | 48 | 22 | 80% | 49% | 22% |
| COAD | 458 | 102 | 88% | 76% | 62% |
| DLBC | 48 | 9 | 90% | 80% | 64% |
| ESCA | 185 | 77 | 76% | 40% | 10% |
| GBM | 594 | 492 | 58% | 13% | 6% |
| HNSC | 527 | 224 | 82% | 57% | 46% |
| KICH | 112 | 12 | 96% | 89% | 85% |
| KIRC | 537 | 177 | 90% | 75% | 62% |
| KIRP | 290 | 44 | 94% | 86% | 74% |
| LAML | 188 | 122 | 53% | 27% | 21% |
| LGG | 515 | 126 | 94% | 77% | 62% |
| LIHC | 376 | 132 | 83% | 62% | 48% |
| LUAD | 513 | 184 | 88% | 62% | 40% |
| LUSC | 499 | 216 | 83% | 58% | 47% |
| MESO | 86 | 73 | 67% | 18% | 5% |
| OV | 583 | 349 | 90% | 62% | 34% |
| PAAD | 185 | 100 | 74% | 32% | 20% |
| PCPG | 179 | 6 | 97% | 77% | 77% |
| PRAD | 500 | 10 | 100% | 98% | 97% |
| READ | 170 | 28 | 93% | 81% | 48% |
| SARC | 261 | 99 | 88% | 67% | 54% |
| SKCM | 462 | 222 | 93% | 71% | 58% |
| STAD | 440 | 175 | 75% | 45% | 33% |
| TGCT | 247 | 7 | 99% | 98% | 96% |
| THCA | 507 | 16 | 99% | 96% | 92% |
| THYM | 123 | 9 | 98% | 93% | 89% |
| UCEC | 547 | 91 | 95% | 82% | 76% |
| UCS | 57 | 35 | 76% | 37% | 28% |
| UVM | 78 | 33 | 90% | 66% | 46% |

### PFI — Progression-Free Interval

| project | N (populated) | events | 1-yr | 3-yr | 5-yr |
| --- | --- | --- | --- | --- | --- |
| ACC | 92 | 49 | 67% | 49% | 42% |
| BLCA | 412 | 176 | 70% | 46% | 39% |
| BRCA | 1097 | 145 | 96% | 88% | 78% |
| CESC | 307 | 72 | 88% | 68% | 63% |
| CHOL | 48 | 23 | 54% | 45% | 45% |
| COAD | 458 | 123 | 84% | 65% | 58% |
| DLBC | 48 | 12 | 80% | 75% | 66% |
| ESCA | 185 | 87 | 62% | 31% | 26% |
| GBM | 594 | 504 | 29% | 8% | 2% |
| HNSC | 527 | 199 | 76% | 60% | 48% |
| KICH | 112 | 16 | 91% | 87% | 84% |
| KIRC | 537 | 163 | 84% | 73% | 64% |
| KIRP | 290 | 57 | 88% | 79% | 72% |
| LAML | 188 | 0 | 100% | 100% | 100% |
| LGG | 515 | 193 | 84% | 56% | 41% |
| LIHC | 376 | 185 | 64% | 38% | 24% |
| LUAD | 513 | 210 | 82% | 49% | 37% |
| LUSC | 499 | 148 | 84% | 64% | 54% |
| MESO | 86 | 64 | 48% | 7% | 7% |
| OV | 583 | 413 | 69% | 22% | 13% |
| PAAD | 185 | 114 | 61% | 21% | 15% |
| PCPG | 179 | 21 | 95% | 82% | 78% |
| PRAD | 500 | 112 | 89% | 75% | 65% |
| READ | 170 | 39 | 87% | 64% | 37% |
| SARC | 261 | 138 | 66% | 45% | 39% |
| SKCM | 462 | 314 | 75% | 50% | 34% |
| STAD | 440 | 162 | 72% | 48% | 37% |
| TGCT | 247 | 49 | 89% | 82% | 80% |
| THCA | 507 | 52 | 96% | 88% | 84% |
| THYM | 123 | 22 | 91% | 83% | 76% |
| UCEC | 547 | 124 | 91% | 74% | 68% |
| UCS | 57 | 37 | 52% | 28% | 28% |
| UVM | 78 | 36 | 77% | 60% | 35% |

### DFI — Disease-Free Interval

| project | N (populated) | events | 1-yr | 3-yr | 5-yr |
| --- | --- | --- | --- | --- | --- |
| ACC | 69 | 28 | 78% | 62% | 53% |
| BLCA | 190 | 31 | 93% | 76% | 61% |
| BRCA | 950 | 79 | 98% | 92% | 85% |
| CESC | 196 | 29 | 95% | 81% | 76% |
| CHOL | 31 | 10 | 68% | 64% | 64% |
| COAD | 327 | 49 | 93% | 78% | 65% |
| DLBC | 37 | 6 | 93% | 87% | 76% |
| ESCA | 143 | 48 | 71% | 49% | 40% |
| GBM | 5 | 3 | 75% | 0% | — |
| HNSC | 193 | 42 | 85% | 77% | 54% |
| KICH | 71 | 6 | 96% | 93% | 86% |
| KIRC | 117 | 15 | 95% | 90% | 83% |
| KIRP | 184 | 27 | 93% | 86% | 79% |
| LGG | 139 | 23 | 96% | 70% | 70% |
| LIHC | 322 | 148 | 67% | 44% | 29% |
| LUAD | 398 | 133 | 86% | 59% | 49% |
| LUSC | 435 | 102 | 89% | 71% | 62% |
| MESO | 17 | 10 | 63% | 0% | 0% |
| OV | 338 | 231 | 81% | 28% | 17% |
| PAAD | 111 | 52 | 72% | 32% | 25% |
| PCPG | 160 | 5 | 98% | 90% | 90% |
| PRAD | 423 | 75 | 92% | 81% | 70% |
| READ | 127 | 19 | 92% | 72% | 42% |
| SARC | 153 | 66 | 74% | 54% | 51% |
| STAD | 336 | 90 | 81% | 60% | 52% |
| TGCT | 166 | 35 | 88% | 81% | 80% |
| THCA | 390 | 29 | 97% | 90% | 89% |
| UCEC | 474 | 77 | 94% | 80% | 76% |
| UCS | 37 | 18 | 67% | 42% | 42% |

### DSS — Disease-Specific Survival

| project | N (populated) | events | 1-yr | 3-yr | 5-yr |
| --- | --- | --- | --- | --- | --- |
| ACC | 92 | 32 | 91% | 73% | 59% |
| BLCA | 412 | 138 | 83% | 56% | 51% |
| BRCA | 1097 | 102 | 98% | 92% | 87% |
| CESC | 307 | 57 | 93% | 74% | 70% |
| CHOL | 48 | 21 | 82% | 50% | 23% |
| COAD | 458 | 80 | 90% | 79% | 72% |
| DLBC | 48 | 4 | 95% | 87% | 87% |
| ESCA | 185 | 53 | 85% | 47% | 28% |
| GBM | 594 | 485 | 58% | 14% | 6% |
| HNSC | 527 | 158 | 87% | 66% | 56% |
| KICH | 112 | 9 | 97% | 92% | 88% |
| KIRC | 537 | 123 | 92% | 81% | 71% |
| KIRP | 290 | 32 | 94% | 89% | 80% |
| LAML | 188 | 122 | 53% | 27% | 21% |
| LGG | 515 | 122 | 94% | 78% | 63% |
| LIHC | 376 | 89 | 90% | 70% | 57% |
| LUAD | 513 | 149 | 90% | 67% | 49% |
| LUSC | 499 | 142 | 90% | 70% | 58% |
| MESO | 86 | 64 | 70% | 22% | 6% |
| OV | 583 | 339 | 90% | 62% | 34% |
| PAAD | 185 | 88 | 78% | 37% | 24% |
| PCPG | 179 | 4 | 98% | 78% | 78% |
| PRAD | 500 | 7 | 100% | 98% | 97% |
| READ | 170 | 23 | 94% | 84% | 53% |
| SARC | 261 | 88 | 90% | 70% | 57% |
| SKCM | 462 | 201 | 94% | 73% | 61% |
| STAD | 440 | 130 | 81% | 56% | 43% |
| TGCT | 247 | 6 | 99% | 98% | 96% |
| THCA | 507 | 13 | 99% | 96% | 94% |
| THYM | 123 | 4 | 99% | 93% | 93% |
| UCEC | 547 | 63 | 97% | 86% | 81% |
| UCS | 57 | 33 | 78% | 39% | 29% |
| UVM | 78 | 26 | 91% | 74% | 53% |

### Reading the tables

- **N (populated)** — patients with both event and time populated for that endpoint. For DFI especially this is much smaller than the project N (most patients lack a disease-free signal at end of first course; see Section 9 deep-dive).
- **events** — count of `event=1` patients in the populated subset.
- **1-yr / 3-yr / 5-yr** — Kaplan-Meier survival probability at that timepoint, accounting for right-censoring. `—` flags projects whose follow-up is shorter than the timepoint (e.g. mostly post-2018 cases).
- **OS vs DSS** — DSS censors cancer-unrelated deaths, so DSS rates are equal to or higher than OS rates at the same timepoint within a project.
- **OS vs PFI** — PFI events are generally earlier than OS events (Liu: *"PFI is generally considered a more informative endpoint"*); PFI rates at 1-year are typically lower than OS rates at 1-year.

### Cross-reference to Liu's reliability assessment

Liu's Table 3 marks each (project, endpoint) as recommended / acceptable / not recommended based on event count and assumption tests. The N-populated column above is the relevant input — projects with very low N for an endpoint don't pass Liu's reliability bar. See the original paper Table 3 for the per-cell recommendations; we don't reproduce that here because it's an editorial overlay rather than a reproducible computation.
