"""Section 11 — Which survival values to use.

Pure-prose section with concrete code-level guidance for picking between
Liu's curated values (published verbatim in tcga-pancanatlas-cdr) and our
re-derived `<ep>_event`/`_time` (published in the `survival_derived`
struct/table).

Writes: sections/11_usage.md
"""

from __future__ import annotations

from pathlib import Path

HERE = Path(__file__).parent
OUT = HERE / "sections" / "11_usage.md"


REPORT = """\
## Section 11 — Which survival values to use

The two sources compared in this report are published separately:

- **Liu's values (curated, frozen)** — [`gabrielaltay/tcga-pancanatlas-cdr`][cdr], config `cdr`: the workbook's `TCGA-CDR` sheet verbatim, one row per patient keyed by `bcr_patient_barcode`, with columns `OS`, `OS.time`, `DSS`, …. Covers the 11,160 patients of the 2018 freeze.
- **Our re-derived values (live)** — `os_event`, `os_time`, … for all 11,428 current patients, keyed by `case_submitter_id`. Nested as the `survival_derived` struct in each [`gabrielaltay/tcga-patients-open`][patients] row, and flat as the `<project>_survival_derived` configs of [`gabrielaltay/tcga-tabular-open`][tabular].

A TCGA patient barcode is the same in both, so `case_submitter_id` joins to `bcr_patient_barcode`.

### Pick based on what you need

**For exact reproducibility against Liu et al.** Use Liu's values. Your results will line up with the paper.

```python
from datasets import load_dataset

liu = load_dataset("gabrielaltay/tcga-pancanatlas-cdr", "cdr", split="train").to_pandas()
liu = liu[(liu.type == "LUAD") & (liu.Redaction != "Redacted")]
# OS curve from OS / OS.time → reproduces Liu's published numbers.
```

**For maximum cohort size and current vital status.** Use the re-derived values. They include the post-freeze patients and reflect the GDC's current data — patients who died after 2018 show as dead.

```python
ours = load_dataset(
    "gabrielaltay/tcga-tabular-open", "TCGA_LUAD_survival_derived", split="train"
).to_pandas()
modern_cohort = ours.dropna(subset=["os_event"])
# OS curve uses current vital status; includes post-freeze patients.
```

**For audit / sanity-check.** Join the two and compare. Disagreement is a red flag worth investigating, especially for OS (where it's almost always data drift; see Section 6).

```python
both = ours.merge(liu, left_on="case_submitter_id", right_on="bcr_patient_barcode")
both = both.dropna(subset=["os_event", "OS"])
mismatches = both[both.os_event != both.OS]
# Each mismatch is a patient whose vital status changed since Liu's 2018 freeze.
```

### Rules of thumb per endpoint

Agreement rate = both correctly NA, or both populated and event direction agrees within 30 days.

- **OS** — 98% agreement; pick whichever source fits your time anchor (Liu's freeze vs current).
- **DSS** — 93% agreement; Liu flagged this as approximate, re-derived value is no more accurate. Use OS instead unless cancer-specific death matters.
- **PFI** — 96% agreement; re-derived is past Liu's reliability bar; good substitute for Liu's `PFI` with extended cohort.
- **DFI** — 90% agreement (up from 77% before the supplement integration). Where both Liu and we populated, event-direction agreement is **99.7%**. The bulk of the remaining 10% is patients where we have *extra coverage* Liu didn't have, not contradictions. For clean Liu reproduction use Liu's `DFI`; for broader coverage including post-2018 use `dfi_event`.

### Where the BCR biotab data lives

The Clinical Supplement biotab data that powers our DFI re-derivation is also surfaced as 7 configs per project (dashes in the project id become underscores, e.g. `TCGA_LUAD_clinical_supplement_patient`) in [`gabrielaltay/tcga-tabular-open`][tabular]:

- `<project>_clinical_supplement_patient` — initial BCR patient form
- `<project>_clinical_supplement_follow_up` — BCR follow-up encounters (one or more form versions per project)
- `<project>_clinical_supplement_nte` — new tumor events
- `<project>_clinical_supplement_drug` — drug records (drug name, dosage, response)
- `<project>_clinical_supplement_radiation` — radiation records
- `<project>_clinical_supplement_ablation` — ablation procedures (LIHC only)
- `<project>_clinical_supplement_omf` — Other Mutation Files (germline)

Per-project schemas (each project ships only the columns its biotab forms contain — BLCA has BCG-related fields, CHOL/LIHC have hepatic markers, etc.). Cross-project queries union with NULL padding via `concatenate_datasets`.

[cdr]: https://huggingface.co/datasets/gabrielaltay/tcga-pancanatlas-cdr
[patients]: https://huggingface.co/datasets/gabrielaltay/tcga-patients-open
[tabular]: https://huggingface.co/datasets/gabrielaltay/tcga-tabular-open

## Conclusions

1. **OS, DSS, PFI reproduce strongly** (98% / 93% / 96% match against Liu's curated values), with most disagreement explainable as data drift since the 2018 freeze.
2. **DFI is now usable** at **90.1% agreement** (up from 77.2% pre-supplement integration). Where both Liu and we populated, event-direction agreement is **99.7%**. The under-population gap shrank from 1,625 patients to 52 patients after we started fetching BCR biotab Clinical Supplements.
3. **Both sources earn their keep**. For users who want Liu's frozen values for direct reproducibility, `tcga-pancanatlas-cdr` is verbatim. For users who want a survival cohort that includes 2018+ patients and reflects current vital status, the re-derived values extend coverage.
4. **The BCR biotab integration is reusable**. Other Pan-Cancer Atlas papers that read BCR-original fields (rather than the harmonized API) can now reproduce against current data without hitting the same wall.
5. **Validation as a standing practice**. This report documents reproducing one of TCGA Pan-Cancer Atlas issue's headline papers from raw GDC data. The same template applies to Hoadley et al. 2018 (iClusters) and other Pan-Cancer Atlas reproductions — see [`dev_todo/reproduce_validate_program.md`](../../dev_todo/reproduce_validate_program.md).
"""


def main() -> None:
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(REPORT)
    print(f"Wrote {OUT.relative_to(HERE.parent.parent)}")


if __name__ == "__main__":
    main()
