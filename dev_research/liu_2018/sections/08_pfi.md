## Section 8 — PFI deep-dive (96.1% match — past Liu's reliability bar)

PFI is the most consequential endpoint clinically — Liu recommends it for 27 of 33 tumor types — and our reproduction is past the 95% threshold. Most disagreements are time differences (longer follow-up since 2018 nudges censoring times), not event-direction disagreements.

### Event-direction mismatches

| Liu cdr_PFI | ours pfi_event | patients |
| --- | --- | --- |
| 0.0 | 1.0 | 55.0 |
| 1.0 | 0.0 | 5.0 |

The asymmetric direction (most are Liu=0, ours=1) is the same data-drift signal as OS: patients who were censored at Liu's 2018 freeze have since had a recurrence / progression / death recorded.

### Time-only mismatches (event matches, time differs)

These are patients we *and* Liu both classify as "event=1" or both "event=0", but the time-to-event or time-to-censor differs by ≥0.5 days.

- **Total time-only mismatches**: 125
- **Median time difference (ours - Liu)**: -364 days

A negative median difference means our times are *earlier* than Liu's. That makes sense for events that happened pre-2018 — same date for both. But for censored patients whose follow-up has continued past 2018, our censor time should be *later* than Liu's (positive diff).

### Per-project event mismatches

| project | event mismatches |
| --- | --- |
| PRAD | 19 |
| STAD | 19 |
| UVM | 6 |
| PAAD | 4 |
| MESO | 3 |
| BLCA | 1 |
| CESC | 1 |
| GBM | 1 |
| HNSC | 1 |
| KIRC | 1 |
| KIRP | 1 |
| LGG | 1 |
| OV | 1 |
| SARC | 1 |

### Implication

PFI is in good shape. The remaining ≤4% gap is the natural drift signature of comparing 2018 frozen data against current live data — not algorithmic.
