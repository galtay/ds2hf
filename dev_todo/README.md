# dev_todo

Open engineering tasks that don't fit a current iteration. Plain markdown
files, no formal tracking. Items here describe *what* and *why*; pick one up
when there's room.

## Contents

- [`reproduce_validate_program.md`](reproduce_validate_program.md) — standing
  practice of reproducing published TCGA pan-cancer results from raw data and
  validating per-case agreement. Liu 2018 CDR wired in as the template;
  Hoadley 2018 iClusters queued as the next target.
- [`ssGSEA.md`](ssGSEA.md) — pathway activity as a derived molecular
  modality. Scoping plus a 100-sample proof-of-concept; the implementation
  and its GSVA validation are now in the repo
  (`ds2hf_pipeline.tcga.ssgsea`, `dev_research/ssgsea/`). What remains is the
  cohort-wide run and the published table.
