# ssGSEA — validation against Bioconductor GSVA

`ds2hf_pipeline.tcga.ssgsea` is a Python transcription of GSVA's `ssgsea.R`.
This directory holds what makes that transcription trustworthy: a small
deterministic fixture and the scores GSVA itself produced for it.

This follows the [reproduce-and-validate
pattern](../../dev_todo/reproduce_validate_program.md) already used for the
Liu 2018 CDR — an independent reference stream, a re-derived stream, and a
committed comparison — with GSVA playing the role Liu's workbook plays.

## Why not just depend on GSVA

The pipeline is a Python/uv workspace and adding R + Bioconductor to it
would end the "one `uv sync`" story for every consumer, to run an algorithm
that is about forty lines of numpy. So R is used **once**, offline, to
produce ground truth, and never at build time.

The transcription earns that: on a 100-sample × 19,938-gene TCGA matrix
scored against MSigDB Hallmark, agreement with GSVA 2.6.6 was

| | raw | normalized |
|---|---|---|
| Pearson / Spearman | 1.0000000000 | 1.0000000000 |
| max relative difference | 4.8e-13 | 6.0e-13 |
| per-sample pathway ordering | 100/100 identical | — |

Residuals are float64 summation-order noise.

## Files

| file | what it is |
|---|---|
| `fixture_expr.tsv` | 600 genes × 8 samples, seeded lognormal with ~31% exact zeros so the tie-handling path is exercised |
| `fixture_sets.gmt` | 5 synthetic sets sized 5 / 12 / 30 / 80 / 200 — the 5-gene one exists to prove the post-mapping size filter drops it |
| `fixture_expected_raw.csv` | raw ssGSEA scores from GSVA 2.6.6, `minSize=10 maxSize=500 alpha=0.25 normalize=FALSE` |
| `fixture_oracle.R` | the script that produced them |

`packages/ds2hf-pipeline/tests/test_ssgsea.py::test_matches_gsva_reference`
asserts we reproduce `fixture_expected_raw.csv` to within 1e-9 relative.

## Regenerating the reference

Needs Docker running; nothing else. The Bioconductor install dominates the
runtime (10–20 min); the scoring itself is instant.

```sh
cd dev_research/ssgsea
docker run --rm -v "$PWD":/work -w /work \
  bioconductor/bioconductor_docker:RELEASE_3_23 \
  bash -c "R -e 'BiocManager::install(c(\"GSVA\",\"GSEABase\"), ask=FALSE, update=FALSE)' \
           && Rscript fixture_oracle.R"
```

Only regenerate deliberately — if the numbers move, that is either a GSVA
behaviour change worth understanding or a bug in the transcription. Do not
refresh the fixture to make a failing test pass.

## Things worth not rediscovering

- **alpha weights ranks, not expression** (`Ra <- R^alpha`). Any strictly
  monotonic transform of expression is therefore a no-op: ssGSEA(TPM) and
  ssGSEA(log1p(TPM)) are bit-identical. Don't log-transform before scoring.
- **GSVA does not drop constant rows for ssGSEA** despite warning about
  them. Our agreement was measured while retaining all 19,938 genes,
  including 504 all-zero ones. Dropping them would shift every rank.
- **Normalization is a single global scalar** (`score / (max - min)` over
  the whole matrix), so chunked scoring is exactly decomposable — but the
  divisor depends on the cohort *and* the gene-set collection scored
  together. Adding MSigDB Reactome to a Hallmark run widened it ~49% on
  TCGA data. Publish raw scores plus the divisor; see
  `ssgsea.normalize_global`.
- **"Unrelated pathway" is not a valid negative control.** ssGSEA ranks
  within a sample, so pathways are positively correlated by construction
  (mean pairwise r ≈ +0.19 across the 50 Hallmarks). Use size-matched
  random gene sets instead.
