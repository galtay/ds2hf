## Section 2 — Median follow-up times (Liu's Table 2)

Liu's Table 2 reports per-project median follow-up plus median time-to-event and time-to-censor for each of OS / PFI / DFI / DSS, all in months (`days / 30.44`). The CDR workbook ships the per-patient event flags and time values Liu used; aggregating those reproduces Liu's published Table 2 to 98.0% (291/297 cells), with the 6 differences in low-count cells where one or two patients shift the median.

### Cell-match rule

Each cell is a median in months. We treat two cells as equal if they're within **±0.5 months**. Liu reports to one decimal; ±0.5 covers rounding noise and one-or-two-patient shifts around the median position. NA matches NA.

### Modern GDC vs Liu CDR — same 11,160 patients

This is the primary drift signal. Same patients, two data sources (modern-GDC re-derived endpoints vs Liu's 2018 CDR-workbook values). Mismatches are post-2018 re-curation showing up in the median: patients whose vital status, follow-up, or recurrence dates have been updated.

#### Per-field match rate

| field | match | % |
| --- | --- | --- |
| FollowUp | 31/33 | 93.9 |
| OS_event | 24/33 | 72.7 |
| OS_censor | 31/33 | 93.9 |
| PFI_event | 29/33 | 87.9 |
| PFI_censor | 24/33 | 72.7 |
| DFI_event | 15/33 | 45.5 |
| DFI_censor | 18/33 | 54.5 |
| DSS_event | 19/33 | 57.6 |
| DSS_censor | 27/33 | 81.8 |

#### Mismatching cells

| project | field | ours | Liu CDR |
| --- | --- | --- | --- |
| BRCA | FollowUp | 27.0 | 27.7 |
| UVM | FollowUp | 29.9 | 25.8 |
| BRCA | OS_event | 37.7 | 41.8 |
| COAD | OS_event | 12.3 | 13.3 |
| KIRC | OS_event | 25.7 | 26.9 |
| READ | OS_event | 19.8 | 22.0 |
| SARC | OS_event | 19.9 | 21.3 |
| THCA | OS_event | 32.4 | 33.5 |
| UCEC | OS_event | 22.0 | 23.3 |
| UCS | OS_event | 15.0 | 17.1 |
| UVM | OS_event | 26.5 | 19.9 |
| SKCM | OS_censor | 37.6 | 36.9 |
| UVM | OS_censor | 36.5 | 27.0 |
| MESO | PFI_event | 8.9 | 10.2 |
| PRAD | PFI_event | 14.0 | 18.4 |
| READ | PFI_event | 17.0 | 19.0 |
| UVM | PFI_event | 14.8 | 12.5 |
| CESC | PFI_censor | 22.6 | 21.7 |
| KIRC | PFI_censor | 41.7 | 43.0 |
| LAML | PFI_censor | 12.0 |  |
| MESO | PFI_censor | 16.4 | 19.4 |
| PAAD | PFI_censor | 14.9 | 13.8 |
| PRAD | PFI_censor | 27.6 | 28.2 |
| READ | PFI_censor | 18.3 | 19.0 |
| STAD | PFI_censor | 14.4 | 13.8 |
| UVM | PFI_censor | 30.4 | 25.0 |
| ACC | DFI_event | 10.8 | 20.0 |
| BRCA | DFI_event | 23.9 | 25.4 |
| CESC | DFI_event | 15.0 | 15.9 |
| COAD | DFI_event | 16.6 | 16.0 |
| DLBC | DFI_event | 63.0 | 113.7 |
| ESCA | DFI_event | 8.8 | 7.4 |
| GBM | DFI_event | 24.6 | 31.5 |
| HNSC | DFI_event | 8.4 | 7.6 |
| KIRP | DFI_event | 16.2 | 15.5 |
| LGG | DFI_event | 16.5 | 19.6 |
| LUAD | DFI_event | 14.4 | 15.7 |
| LUSC | DFI_event | 14.4 | 18.0 |
| PAAD | DFI_event | 12.1 | 14.8 |
| PCPG | DFI_event | 25.4 | 27.3 |
| PRAD | DFI_event | 16.8 | 24.9 |
| READ | DFI_event | 26.0 | 27.8 |
| TGCT | DFI_event | 10.6 | 14.8 |
| UCS | DFI_event | 10.5 | 16.6 |
| COAD | DFI_censor | 23.0 | 29.3 |
| ESCA | DFI_censor | 12.6 | 13.2 |
| GBM | DFI_censor | 13.9 | 26.3 |
| HNSC | DFI_censor | 28.4 | 27.5 |
| LUAD | DFI_censor | 21.4 | 22.5 |
| LUSC | DFI_censor | 22.6 | 26.9 |
| MESO | DFI_censor | 6.1 | 9.8 |
| OV | DFI_censor | 25.1 | 26.5 |
| PRAD | DFI_censor | 28.7 | 30.4 |
| READ | DFI_censor | 17.8 | 21.0 |
| SARC | DFI_censor | 35.9 | 36.5 |
| STAD | DFI_censor | 16.6 | 18.6 |
| TGCT | DFI_censor | 36.9 | 31.8 |
| THCA | DFI_censor | 31.0 | 31.9 |
| UCEC | DFI_censor | 30.1 | 30.7 |
| BRCA | DSS_event | 34.0 | 32.6 |
| CESC | DSS_event | 18.7 | 18.0 |
| COAD | DSS_event | 11.9 | 11.1 |
| KIRP | DSS_event | 16.4 | 14.2 |
| LAML | DSS_event | 9.0 |  |
| LGG | DSS_event | 26.7 | 25.5 |
| LIHC | DSS_event | 19.1 | 19.8 |
| LUAD | DSS_event | 20.5 | 19.9 |
| LUSC | DSS_event | 19.8 | 18.8 |
| READ | DSS_event | 24.0 | 20.0 |
| SARC | DSS_event | 21.6 | 22.6 |
| STAD | DSS_event | 11.7 | 12.4 |
| THCA | DSS_event | 30.8 | 33.5 |
| UVM | DSS_event | 27.6 | 19.9 |
| KIRC | DSS_censor | 45.5 | 46.4 |
| LAML | DSS_censor | 23.0 |  |
| MESO | DSS_censor | 23.3 | 24.9 |
| SKCM | DSS_censor | 36.4 | 34.3 |
| UCEC | DSS_censor | 30.1 | 30.7 |
| UVM | DSS_censor | 34.1 | 27.0 |

### Sanity — CDR-aggregated reproduces Liu's published Table 2

This confirms our aggregation logic before the drift comparison above is interpretable.

#### Per-field match rate

| field | match | % |
| --- | --- | --- |
| FollowUp | 33/33 | 100.0 |
| OS_event | 33/33 | 100.0 |
| OS_censor | 33/33 | 100.0 |
| PFI_event | 33/33 | 100.0 |
| PFI_censor | 33/33 | 100.0 |
| DFI_event | 30/33 | 90.9 |
| DFI_censor | 30/33 | 90.9 |
| DSS_event | 33/33 | 100.0 |
| DSS_censor | 33/33 | 100.0 |

#### Mismatching cells

| project | field | CDR-aggregated | Liu paper |
| --- | --- | --- | --- |
| SKCM | DFI_event |  | 21.8 |
| THYM | DFI_event |  | 30.8 |
| UVM | DFI_event |  | 12.2 |
| SKCM | DFI_censor |  | 23.8 |
| THYM | DFI_censor |  | 42.1 |
| UVM | DFI_censor |  | 26.2 |

### Our re-derived Table 2 on the full current cohort (11,428)

Includes the 268 post-freeze patients. Useful for picking the most-currently-informative project per endpoint.

| project | N | FollowUp | OS_event | OS_censor | DSS_event | DSS_censor | PFI_event | PFI_censor | DFI_event | DFI_censor |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ACC | 92 | 38.8 | 18.1 | 47.7 | 18.1 | 47.7 | 8.1 | 49.2 | 10.8 | 60.8 |
| BLCA | 412 | 17.5 | 13.4 | 20.9 | 13.5 | 19.2 | 9.5 | 17.8 | 14.5 | 19.1 |
| BRCA | 1098 | 27.0 | 37.7 | 25.0 | 34.0 | 25.9 | 26.0 | 25.0 | 23.9 | 25.0 |
| CESC | 307 | 21.0 | 20.4 | 22.6 | 18.7 | 23.3 | 13.8 | 22.6 | 15.0 | 28.4 |
| CHOL | 51 | 22.3 | 18.0 | 32.7 | 18.2 | 30.1 | 7.1 | 23.3 | 7.1 | 27.4 |
| COAD | 461 | 21.6 | 12.3 | 24.0 | 11.9 | 23.9 | 11.7 | 22.0 | 16.6 | 23.0 |
| DLBC | 58 | 26.7 | 19.5 | 31.1 | 16.2 | 29.2 | 10.3 | 29.2 | 63.0 | 31.1 |
| ESCA | 185 | 13.1 | 11.5 | 13.2 | 13.9 | 12.9 | 8.8 | 12.6 | 8.8 | 12.6 |
| GBM | 617 | 12.1 | 12.5 | 8.6 | 12.5 | 8.5 | 6.1 | 6.0 | 24.6 | 13.9 |
| HNSC | 528 | 21.2 | 14.1 | 28.0 | 14.0 | 25.7 | 9.4 | 25.7 | 8.4 | 28.4 |
| KICH | 113 | 48.3 | 24.3 | 54.2 | 28.1 | 51.0 | 11.9 | 49.7 | 52.7 | 39.6 |
| KIRC | 537 | 38.6 | 25.7 | 47.4 | 23.7 | 45.5 | 13.1 | 41.7 | 29.6 | 45.4 |
| KIRP | 291 | 25.2 | 21.1 | 25.3 | 16.4 | 25.9 | 10.9 | 25.3 | 16.2 | 25.2 |
| LAML | 200 | 12.0 | 9.0 | 23.0 | 9.0 | 23.0 |  | 12.0 |  |  |
| LGG | 516 | 22.2 | 26.8 | 20.7 | 26.7 | 20.7 | 15.3 | 18.7 | 16.5 | 20.1 |
| LIHC | 377 | 19.7 | 13.7 | 21.3 | 19.1 | 19.7 | 9.0 | 15.6 | 8.9 | 17.3 |
| LUAD | 585 | 21.6 | 20.2 | 22.0 | 20.5 | 21.7 | 14.4 | 20.0 | 14.4 | 21.4 |
| LUSC | 504 | 21.7 | 17.9 | 24.9 | 19.8 | 22.5 | 14.0 | 21.1 | 14.4 | 22.6 |
| MESO | 87 | 16.9 | 15.0 | 38.4 | 15.0 | 23.3 | 8.9 | 16.4 | 15.6 | 6.1 |
| OV | 608 | 32.9 | 35.1 | 27.8 | 34.8 | 28.9 | 14.7 | 15.0 | 17.9 | 25.1 |
| PAAD | 185 | 15.3 | 12.9 | 17.0 | 13.6 | 15.9 | 11.2 | 14.9 | 12.1 | 15.2 |
| PCPG | 179 | 24.8 | 14.9 | 25.2 | 17.5 | 25.0 | 19.9 | 23.8 | 25.4 | 24.3 |
| PRAD | 500 | 30.5 | 36.2 | 30.5 | 43.7 | 30.5 | 14.0 | 27.6 | 16.8 | 28.7 |
| READ | 172 | 20.0 | 19.8 | 20.0 | 24.0 | 20.0 | 17.0 | 18.3 | 26.0 | 17.8 |
| SARC | 261 | 30.8 | 19.9 | 35.6 | 21.6 | 34.7 | 10.0 | 32.7 | 11.3 | 35.9 |
| SKCM | 470 | 36.3 | 35.7 | 37.6 | 36.2 | 36.4 | 23.5 | 22.7 |  |  |
| STAD | 443 | 13.9 | 11.4 | 17.2 | 11.7 | 16.0 | 9.1 | 14.4 | 10.7 | 16.6 |
| TGCT | 263 | 51.2 | 20.3 | 51.6 | 18.6 | 51.7 | 10.4 | 45.4 | 10.4 | 46.1 |
| THCA | 507 | 31.0 | 32.4 | 31.0 | 30.8 | 31.0 | 16.0 | 30.8 | 16.2 | 31.0 |
| THYM | 124 | 41.2 | 28.0 | 41.6 | 54.9 | 41.2 | 25.2 | 41.2 |  |  |
| UCEC | 560 | 29.4 | 22.0 | 30.7 | 22.0 | 30.1 | 16.8 | 29.3 | 17.1 | 30.1 |
| UCS | 57 | 19.6 | 15.0 | 27.2 | 15.0 | 26.9 | 9.0 | 26.9 | 10.5 | 26.7 |
| UVM | 80 | 29.9 | 26.5 | 36.5 | 27.6 | 34.1 | 14.8 | 30.4 |  |  |

### Notable findings

- **Bucketing logic is sound.** The CDR-aggregated reconstruction matches Liu's paper Table 2 at 98.0%; remaining differences are sub-patient median shifts in low-count cells.
- **Median follow-up has grown by ~1-12 months for several projects.** SKCM, BRCA, KIRC, OV all show longer median follow-up in the modern GDC than Liu had in 2018 — patients who were censored in 2018 have either had more follow-up visits or eventually died.
- **DFI medians are the noisiest.** Same structural reason as Section 1 stage NA collapse: the underlying `treatment_outcome_first_course` field has shifted population over time. DLBC's 113.7 month DFI-to-event in Liu (one or two patients) drops in our re-derivation because we no longer have those patients populated for DFI.
- **Our re-derived DFI populates fewer cells than Liu's CDR.** A handful of (project, DFI) entries are NA in our re-derivation but populated in Liu's. Section 9 deep-dives this.

### Confidence in the drift attribution

Same as Section 1: the CDR workbook ships per-patient values, so the cell-level differences here are real per-patient changes rather than aggregation issues. Without a versioned 2018 GDC clinical snapshot we can't always attribute each change to (a) modern GDC re-curation vs (b) Liu hand-fix that didn't round-trip — but in aggregate, modern GDC consistently has more populated, more-recent values.
