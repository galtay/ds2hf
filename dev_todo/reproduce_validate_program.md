# Reproduce-and-validate program for prior TCGA pan-cancer results

A standing engineering practice for this project: every time we add a
derived column or computed feature to the pipeline, we should be able to
reproduce the corresponding result from a published TCGA paper from the
same raw data, and validate that we get the same answer for the cases
they covered.

This emerged from the CDR (Liu et al. 2018) work in
`ds2hf_pipeline.tcga.survival`. That pattern — curated stream alongside
re-derived stream, with a per-project per-endpoint match-rate report
(`scripts/validate_survival.py`) — is the template.

## Why

We don't gain confidence from passing tests on synthetic data alone. The
GDC schema, field semantics, and modality conventions have all drifted
since 2018; reproducing a published result on the modern data is the
strongest available check that our transformations preserve the
underlying clinical/molecular truth.

## Pattern (the template)

For each paper:

1. **Curated stream**: download the supplementary table verbatim and
   surface its values as `<paper>_*` columns on the relevant entity
   (patient rows for clinical, sample/aliquot rows for molecular). md5
   pin the file. Mark unmatched cases with an audit flag (e.g.
   `cdr_matched`) rather than silently dropping.
2. **Re-derived stream**: implement the paper's documented algorithm
   against current GDC data, returning the same column shape. Aim for
   full coverage including post-freeze cases the curated stream misses.
3. **Validation script**: `scripts/validate_<paper>.py` reports per-project
   per-endpoint match rates (match / mismatch / both-NA / curated-only-NA
   / derived-only-NA). Iterate the algorithm until match rates are
   acceptable or the divergence is documented.
4. **Canary test**: at least one named patient in `tests/test_<paper>.py`
   with the published expected values, locked in as a regression guard.
5. **Memory entry**: per-endpoint match rates and the structural reasons
   for any persistent divergence, so future sessions don't relitigate
   work that's already settled.

## Open targets

- **Liu et al. 2018 (CDR)** — *in progress*, wired into the pipeline.
  Match rates: OS 98% / DSS 93% / PFI 96% / DFI 47%. DFI gap is
  structural — Liu's per-patient inclusion criteria for DFI go beyond
  the paper's documented spec; likely live in Table S2 we'd need to
  reverse-engineer or get the SAS source for.
  See [`dev_research/liu_2018/notes.md`](../dev_research/liu_2018/notes.md).

- **Hoadley et al. 2018 (Cell-of-Origin / iClusters)** — *next target*.
  28 integrative clusters across 10,000 tumors / 33 cancer types.
  Curated stream (the published Table S2 iCluster assignments) is
  cheap; re-derivation requires modalities we don't yet ship (DNA
  methylation, copy number, miRNA-Seq, RPPA) so it's gated on a
  modalities expansion. Adding the curated stream alone is still
  valuable — pan-X cohort labels (Pan-Squamous, Pan-GYN, Pan-Kidney,
  Pan-GI) are immediately useful for downstream filtering.
  See [`dev_research/hoadley_2018/notes.md`](../dev_research/hoadley_2018/notes.md).

## Future candidates

The Cell 2018 Pan-Cancer Atlas issue ships ~25 papers; many define
derived per-patient or per-sample features (immune subtypes, oncogenic
pathway alterations, DNA damage repair signatures, etc.). Anything with
a supplementary table listing per-case values is a viable target for
this pattern. We'll add a `dev_research/<paper>/` directory per paper
(holding the PDF, notes, and any analysis notebooks) and link them here
as we pick them up.
