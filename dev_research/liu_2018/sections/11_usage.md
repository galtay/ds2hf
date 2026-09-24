## Section 11 — How to use the two streams in practice

Each patient row in `gabrielaltay/tcga-patients-open` (and each `cases` row in `gabrielaltay/tcga-tabular-open`) carries **two parallel streams** of survival annotation:

- **`cdr_*` (curated, frozen)** — Liu's 2018 values lifted verbatim from `TCGA-CDR-SupplementalTableS1.xlsx`. Direct reproducibility against the paper. Filter on `cdr_matched == True` to restrict to the 11,160 patients Liu covered.
- **`{os,dss,pfi,dfi}_event` / `_time` (re-derived, live)** — the same four endpoints recomputed from the current GDC data using Liu's documented algorithm (`ds2hf_pipeline.tcga.survival`), augmented with `treatment_outcome_first_course` from the BCR biotab Clinical Supplements. Full coverage including 268 post-freeze patients.

### Pick based on what you need

**For exact reproducibility against Liu et al.** Use `cdr_*` columns. Filter to `cdr_matched == True`. You get exactly what's in `TCGA-CDR-SupplementalTableS1.xlsx`, and your results will line up bit-for-bit with the paper.

```python
import pyarrow.parquet as pq
patients = pq.read_table("TCGA-LUAD/data.parquet").to_pylist()
liu_cohort = [p for p in patients if p["cdr_matched"]]
# OS curve from cdr_OS / cdr_OS_time → reproduces Liu's published numbers.
```

**For maximum cohort size and current vital status.** Use `{os,dss,pfi,dfi}_event/_time`. Includes ~268 post-freeze patients and reflects the GDC's current data — patients who died after 2018 show as Dead.

```python
patients = pq.read_table("TCGA-LUAD/data.parquet").to_pylist()
modern_cohort = [p for p in patients if p["os_event"] is not None]
# OS curve uses current vital status; includes post-freeze patients.
```

**For audit / sanity-check.** Filter to `cdr_matched == True` and compare the two streams. Disagreement is a red flag worth investigating, especially for OS (where it's almost always data drift; see Section 6).

```python
matched = [p for p in patients if p["cdr_matched"]]
mismatches = [
    p for p in matched
    if p["cdr_OS"] is not None and p["os_event"] is not None
    and p["cdr_OS"] != p["os_event"]
]
# Each mismatch is a patient whose vital status changed since Liu's 2018 freeze.
```

### Rules of thumb per endpoint

Agreement rate = both correctly NA, or both populated and event direction agrees within 30 days.

- **OS** — 98% agreement; pick whichever stream fits your time anchor (Liu's freeze vs current).
- **DSS** — 93% agreement; Liu flagged this as approximate, re-derived value is no more accurate. Use OS instead unless cancer-specific death matters.
- **PFI** — 96% agreement; re-derived is past Liu's reliability bar; good substitute for `cdr_PFI` with extended cohort.
- **DFI** — 90% agreement (up from 77% before the supplement integration). Where both Liu and we populated, event-direction agreement is **99.7%**. The bulk of the remaining 10% is patients where we have *extra coverage* Liu didn't have, not contradictions. For clean Liu reproduction use `cdr_DFI`; for broader coverage including post-2018 use `dfi_event`.

### Where the BCR biotab data lives

The Clinical Supplement biotab data that powers our DFI re-derivation is also surfaced as 7 tables per project in [`gabrielaltay/tcga-tabular-open`][tabular]:

- `<project>/clinical_supplement_patient` — initial BCR patient form
- `<project>/clinical_supplement_follow_up` — BCR follow-up encounters (one or more form versions per project)
- `<project>/clinical_supplement_nte` — new tumor events
- `<project>/clinical_supplement_drug` — drug records (drug name, dosage, response)
- `<project>/clinical_supplement_radiation` — radiation records
- `<project>/clinical_supplement_ablation` — ablation procedures (LIHC only)
- `<project>/clinical_supplement_omf` — Other Mutation Files (germline)

Per-project schemas (each project ships only the columns its biotab forms contain — BLCA has BCG-related fields, CHOL/LIHC have hepatic markers, etc.). Cross-project queries union with NULL padding via `concatenate_datasets`.

[tabular]: https://huggingface.co/datasets/gabrielaltay/tcga-tabular-open

## Conclusions

1. **OS, DSS, PFI reproduce strongly** (98% / 93% / 96% match against Liu's curated values), with most disagreement explainable as data drift since the 2018 freeze.
2. **DFI is now usable** at **90.1% agreement** (up from 77.2% pre-supplement integration). Where both Liu and we populated, event-direction agreement is **99.7%**. The under-population gap shrank from 1,625 patients to 52 patients after we started fetching BCR biotab Clinical Supplements.
3. **Two-stream design earns its keep**. For users who want Liu's frozen values for direct reproducibility, `cdr_*` is verbatim. For users who want a survival cohort that includes 2018+ patients and reflects current vital status, the re-derived columns extend coverage.
4. **The BCR biotab integration is reusable**. Other Pan-Cancer Atlas papers that read BCR-original fields (rather than the harmonized API) can now reproduce against current data without hitting the same wall.
5. **Validation as a standing practice**. This report documents reproducing one of TCGA Pan-Cancer Atlas issue's headline papers from raw GDC data. The same template applies to Hoadley et al. 2018 (iClusters) and other Pan-Cancer Atlas reproductions — see [`dev_todo/reproduce_validate_program.md`](../../dev_todo/reproduce_validate_program.md).
