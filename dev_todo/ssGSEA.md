Yes. I think this is a very tractable extension of the dataset. I would treat it as a **derived molecular modality**: RNA expression in, pathway activity out.

The most important implementation detail is that **ssGSEA is fundamentally sample-wise, but the standard GSVA implementation applies a final normalization using the range of scores over the entire input dataset**. So if we want pan-TCGA scores that are comparable across cancer types, we should not independently normalize each TCGA project. The GSVA source explicitly normalizes using the minimum and maximum over the entire calculated dataset. ([GitHub][1])

## Proposed pipeline

```text
HF TCGA RNA expression
        │
        │ select TPM
        ▼
sample × gene expression
        │
        │ map gene_name → MSigDB symbols
        ▼
expression matrix
genes × samples
        │
        ├── MSigDB Hallmark
        ├── Reactome
        └── WikiPathways
        │
        ▼
ssGSEA
        │
        ▼
pathway × sample scores
        │
        ▼
transpose / join identifiers
        │
        ▼
HF derived pathway dataset
```

### 1. Define the unit of observation

I would calculate a score for **each RNA-seq sample/aliquot**, not each patient.

So internally we'd want something like:

```python
{
    "case_id": "...",
    "sample_id": "...",
    "aliquot_id": "...",
    "project_id": "TCGA-LUAD",
    "sample_type": "Primary Tumor",
    "genes": [...],
    "tpm": [...]
}
```

If a patient has a tumor RNA sample and a normal RNA sample, they get distinct pathway-score records. That's biologically the right level.

### 2. Use `tpm_unstranded`

For ssGSEA I would use the GDC STAR **TPM** column rather than raw counts.

ssGSEA works by ranking genes **within each sample** and computing a weighted random walk over those ranks. The current GSVA implementation uses `alpha=0.25` by default. ([Bioconductor][2])

That gives us a nice property:

```text
TPM               log1p(TPM)
1000      #1      6.91       #1
500       #2      6.22       #2
12        #3      2.56       #3
...
```

Because `log1p()` is monotonic,

```python
ssGSEA(TPM) ≈ ssGSEA(log1p(TPM))
```

apart from numerical/tie handling.

So I'd just use TPM directly and avoid adding an unnecessary transformation.

I would **not** use raw counts, because TPM normalization changes rankings between genes of different lengths; raw counts therefore do not necessarily give the same within-sample gene ordering.

### 3. Use gene symbols initially

MSigDB and your expression data need to use identical identifiers. GSVA explicitly warns that if ordinary matrices/lists are supplied, identifier matching is the user's responsibility. ([Bioconductor][3])

For your data, the easiest route is:

```text
GDC gene_name
      ↓
MSigDB gene symbols
```

rather than:

```text
ENSG00000141510.18
      ↓ annotation mapping
TP53
      ↓
MSigDB
```

I would keep the Ensembl ID around for provenance, but score using `gene_name`.

There are two things we'd explicitly QC here:

```python
genes_with_symbol
genes_with_duplicate_symbol
genes_matching_msigdb
```

For duplicate gene symbols, I'd inspect how frequent they actually are before deciding on a collapse strategy. Ideally they are rare enough that we can define a simple deterministic rule.

### 4. Start with Hallmark

I would **not** start with 13,000 gene sets.

The first release I'd make is the **50 MSigDB Hallmark gene sets**. Hallmarks were specifically constructed as coherent signatures representing well-defined biological states/processes, and they eliminate a lot of redundancy present in larger pathway collections. ([GSEA MSigDB][4])

That gets us:

```text
~11,000 TCGA RNA samples × 50 pathways
```

Examples:

```text
HALLMARK_APOPTOSIS
HALLMARK_E2F_TARGETS
HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
HALLMARK_G2M_CHECKPOINT
HALLMARK_HYPOXIA
HALLMARK_INTERFERON_GAMMA_RESPONSE
HALLMARK_KRAS_SIGNALING_UP
HALLMARK_MTORC1_SIGNALING
HALLMARK_P53_PATHWAY
HALLMARK_TNFA_SIGNALING_VIA_NFKB
...
```

That's immediately useful as a compact representation of TCGA.

Then I'd do:

```text
v1:
    H — Hallmark

v2:
    C2:CP:REACTOME
    C2:CP:WIKIPATHWAYS
```

and perhaps eventually other collections.

Current MSigDB is **2026.1.Hs**, released in January 2026, so I would pin exactly that version in the dataset metadata rather than silently tracking "latest." ([GSEA MSigDB][4])

### 5. Run the canonical Bioconductor implementation

Even though most of your pipeline is Python, for a public reference dataset I lean toward using the **Bioconductor `GSVA` implementation** as the ground truth rather than reimplementing ssGSEA ourselves.

The current API is essentially:

```r
library(GSVA)

par <- ssgseaParam(
    expr,
    gene_sets,
    minSize = 10,
    maxSize = 500,
    alpha = 0.25,
    normalize = FALSE
)

scores <- gsva(par)
```

`expr` is:

```text
genes × samples
```

and the result is:

```text
gene sets × samples
```

GSVA 2.6.6 is the current Bioconductor release, and `ssgseaParam()` currently defaults to `alpha=0.25` and `normalize=TRUE`. ([Bioconductor][5])

Notice that I deliberately put:

```r
normalize = FALSE
```

there.

### 6. Handle normalization ourselves

This is the subtle part.

GSVA's standard ssGSEA does:

```text
raw ssGSEA scores
        ↓
divide by (global maximum - global minimum)
```

and its implementation describes this as normalization using the **entire dataset**. ([GitHub][1])

That means this would be bad:

```text
TCGA-BRCA → ssGSEA normalize
TCGA-LUAD → ssGSEA normalize
TCGA-LGG  → ssGSEA normalize
...
```

because each cancer would get a different scale.

Instead I'd calculate:

```text
TCGA-BRCA ─┐
TCGA-LUAD ─┤
TCGA-LGG  ─┤ → raw ssGSEA → global min/max → normalize once
...        ─┘
```

There are two implementation strategies.

**A. Simplest:** construct the full ~20k gene × ~11k sample matrix and calculate everything in one GSVA call.

That's roughly:

```text
20,000 × 11,000 = 220 million values
```

At float64 that's only ~1.76 GB for the raw matrix, though working memory will be higher. This isn't a particularly large computation on a decent workstation.

**B. More scalable:** calculate `normalize=FALSE` in chunks, retain the raw scores, then reproduce GSVA's final global normalization:

```python
score_range = raw_scores.max() - raw_scores.min()
scores = raw_scores / score_range
```

The GSVA source shows exactly this range-based division. ([GitHub][1])

I actually prefer **B**, because it makes the pipeline easily resumable and lets you process directly from Hugging Face without materializing one monster input object.

### 7. Filtering gene sets

I would impose:

```text
minSize = 10
maxSize = 500
```

**after intersecting the pathway genes with TCGA genes**.

GSVA itself filters unmapped genes and then filters gene sets according to their post-mapping sizes. ([Bioconductor][3])

For every pathway I'd retain metadata such as:

```python
{
    "pathway": "HALLMARK_P53_PATHWAY",
    "msigdb_version": "2026.1.Hs",
    "original_gene_count": 200,
    "matched_gene_count": 198,
    "match_fraction": 0.99
}
```

That is valuable provenance and would make your version better than many precomputed resources.

## Output format

I would probably **not append 50 or 1,500 columns to the existing RNA table**.

Instead make a derived config/table:

```text
ssgsea_hallmark
```

with records like:

```python
{
    "case_id": "....",
    "sample_id": "....",
    "project_id": "TCGA-LUAD",
    "sample_type": "Primary Tumor",

    "pathway_names": [
        "HALLMARK_ADIPOGENESIS",
        "HALLMARK_ALLOGRAFT_REJECTION",
        ...
    ],

    "scores": [
        0.2318,
        -0.0934,
        ...
    ]
}
```

Or, because Hallmark is only 50 dimensions, wide format is also perfectly reasonable:

```text
sample_id
project_id
HALLMARK_APOPTOSIS
HALLMARK_E2F_TARGETS
HALLMARK_G2M_CHECKPOINT
...
```

For Reactome with ~1,000+ features, arrays start looking considerably better.

## Validation

Before publishing it, I'd build a fairly strong validation suite.

For example:

```text
1. TPM vs log1p(TPM)
   → scores should be essentially identical

2. Chunked raw ssGSEA vs single-call raw ssGSEA
   → exactly/nearly identical

3. Global normalization
   → chunked + post-normalization equals normal GSVA call

4. Gene-set coverage
   → inspect matched/original gene counts

5. Biological positive controls
```

For the biological controls, some easy ones are:

```text
MKI67 expression
       ↕ positive correlation
HALLMARK_E2F_TARGETS
HALLMARK_G2M_CHECKPOINT
```

and

```text
CD8A / PTPRC
       ↕
HALLMARK_INTERFERON_GAMMA_RESPONSE
```

We can also look for expected cancer-specific distributions rather than assuming the method is correct simply because the software ran.

## One additional issue: licensing/provenance

Current MSigDB releases are generally under **CC BY 4.0**, but some constituent gene-set collections have additional restrictions, notably some KEGG-derived sets. ([GSEA MSigDB][6])

For a public Hugging Face dataset, I would therefore start with:

```text
Hallmark
Reactome
WikiPathways
```

and explicitly record:

```yaml
gene_set_source: MSigDB
msigdb_version: 2026.1.Hs
method: ssGSEA
implementation: GSVA
gsva_version: 2.6.6
alpha: 0.25
min_size: 10
max_size: 500
normalization: global_range
expression_unit: TPM
gene_identifier: HGNC_symbol
```

That gives us a highly reproducible derived dataset. ([GSEA MSigDB][6])

### What I would build first

The first experiment can be extremely small:

```text
your HF TCGA data
       ↓
100 randomly selected RNA samples
       ↓
TPM × gene symbol
       ↓
MSigDB 2026.1 Hallmark
       ↓
Bioconductor GSVA ssGSEA
       ↓
100 × 50 dataframe
```

Then inspect those 50 scores and the gene-set mapping statistics. Once that looks right, scaling it to all of TCGA is mostly an engineering problem.

I think the next useful step is to look directly at **how your gene-expression records are represented in `tcga-tabular-open` and write the actual extraction → matrix → GSVA pipeline around that schema**. That should only require a fairly small amount of code.

---

# Findings from the 100-sample proof-of-concept (2026-08-12)

Everything below was measured, not assumed. Working code lives outside the
repo for now; promote to `dev_research/ssgsea/` when we build for real.

## Verified against the GSVA source

**The α-weight is applied to ranks, not expression** (`Ra <- R^alpha`, where
`R` is the rank matrix). This makes the log1p argument above *exact* rather
than approximate: any strictly monotonic transform gives identical ranks and
therefore identical scores. Confirmed empirically — `ssGSEA(TPM)` and
`ssGSEA(log1p(TPM))` were bit-identical (max abs diff 0.0).

**Normalization is `scores / (max - min)` accumulated over the entire
matrix** — a pure scaling with no shift, so sign and zero are preserved.
Because it is a single global scalar, strategy B is not merely "more
scalable", it is *provably equivalent* to strategy A. Confirmed: raw scores
computed in 3 chunks were bit-identical to the single call, and normalizing
with an accumulated divisor reproduced the single-call result exactly.

## Decisions now pinned

- **Gene universe: protein-coding only.** This is a real methodological
  choice, not a detail — ssGSEA ranks *within the supplied matrix*, so
  including the ~40k lncRNA/pseudogene rows shifts every rank and changes
  every score. Pin it in the metadata block alongside `alpha`.
- **Duplicate symbols are close to a non-issue.** Measured on the actual
  GENCODE v36 universe: 60,660 genes → 59,427 unique symbols, but the
  duplicates are overwhelmingly non-coding repeats (`Y_RNA` ×756,
  `Metazoa_SRP` ×170). Restricted to protein-coding: 19,962 genes, only 24
  duplicated symbols — and **18 of those are `_PAR_Y` copies that are all
  exactly 0.0 TPM** (GDC's STAR pipeline assigns pseudoautosomal reads to
  the X copy). Dropping `_PAR_Y` is lossless. That leaves **6** genuinely
  ambiguous symbols; collapsing by max TPM is sufficient.
- **Final gene universe: 19,938 symbols.**

## The normalization is a versioning hazard

The divisor is the range over *the entire computed matrix*, so it depends on
dataset composition. Measured: the divisor over 100 samples vs the first 50
differs enough to shift a given sample's normalized score by **2.8%** —
from nothing but cohort membership. Adding Reactome later, or new TCGA
samples, would silently move every published Hallmark score.

**Therefore: ship `score_raw` plus the divisor**, and optionally
`score_normalized` alongside. Consumers can then renormalize against any
cohort, and v1→v2 values stay reproducible. Same instinct as re-deriving
Liu's endpoints instead of baking in his frozen values.

## Measurements

| quantity | value |
|---|---|
| RNA-seq aliquots pan-TCGA | 11,505 (10,517 cases) |
| Gene universe after filtering | 19,938 protein-coding symbols |
| Hallmark sets surviving mapping + size filter | **50 / 50** |
| Gene-set match fraction | median 1.000, min 0.983 |
| TPM matrix, full cohort | 0.92 GB float32 |
| Matrix extraction | ~5 min serial, trivially parallel |
| Scoring throughput | 7 ms/sample → **~1.4 min** for all 11,505, single core |

The whole computation is small. Chunking is a convenience, not a necessity.

## Biological validation (n=100, 28 projects)

Spearman ρ between marker gene TPM and pathway score:

| marker | pathway | ρ |
|---|---|---|
| MKI67 | HALLMARK_G2M_CHECKPOINT | **+0.917** |
| MKI67 | HALLMARK_E2F_TARGETS | **+0.884** |
| COL1A1 | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | **+0.818** |
| CD8A | HALLMARK_ALLOGRAFT_REJECTION | **+0.826** |
| PTPRC | HALLMARK_INFLAMMATORY_RESPONSE | **+0.759** |
| HIF1A | HALLMARK_HYPOXIA | +0.234 |

Tumor vs solid-tissue-normal: E2F / G2M / MYC targets highest in tumour,
TNFA / inflammatory / IL6-JAK-STAT3 highest in normal. Both expected.

Two caveats worth carrying forward:

- **HIF1A is a poor hypoxia marker and the weak ρ is correct biology.** HIF1A
  is regulated by protein stabilisation, not transcription, so its mRNA
  should *not* track the hypoxia program. Don't "fix" this.
- **"Unrelated pathway" is not a valid negative control.** ssGSEA ranks
  within a sample, so pathways are correlated by construction (mean pairwise
  r = +0.19 across the 50 Hallmarks; 21% of pairs exceed |r| 0.5). MKI67 vs
  BILE_ACID_METABOLISM gives ρ = −0.585 — compositionality, not a bug. The
  correct control is **size-matched random gene sets**: |ρ| median 0.090 vs
  the real G2M set's 0.917.

## Implementation: Python, with GSVA as an offline oracle

Pipeline stays pure Python — R never enters `uv sync`. ssGSEA's core is ~40
lines of numpy (rank → descending order → weighted random walk → sum), and
it is transcribed directly from `ssgsea.R`.

GSVA runs **once** in `bioconductor/bioconductor_docker` on these 100
samples; its output is committed as a validation fixture. That is exactly
the `reproduce_validate_program.md` template, with GSVA playing the role
Liu's CDR workbook plays.

### Oracle result: agreement is exact

Ran GSVA **2.6.6** (Bioconductor `RELEASE_3_23`) over the same 19,938 x 100
matrix and the same GMT, `minSize=10 maxSize=500 alpha=0.25`, both
`normalize=FALSE` and `TRUE`:

| | raw | normalized |
|---|---|---|
| Pearson / Spearman vs ours | 1.0000000000 | 1.0000000000 |
| max abs difference | 2.0e-11 (on values ~1e3) | 2.0e-15 |
| max **relative** difference | 4.8e-13 | 6.0e-13 |
| per-sample pathway ordering | 100/100 identical | — |

Residuals are float64 summation-order noise, not algorithmic divergence.
**The Python implementation is the reference-equivalent; R is not needed in
the pipeline.**

One near-miss worth recording: GSVA warns `504 rows with constant values`
(the all-zero genes in this 100-sample subset). Because we agree to 1e-13
while retaining all 19,938 genes, GSVA demonstrably does **not** drop
constant rows for ssGSEA — the warning is informational. Had it dropped
them the gene universe would differ, ranks would shift, and scores would
diverge materially. Re-check this if the GSVA version changes.

`gseapy.ssgsea` is a viable alternative to our own implementation — with
`sample_norm_method='rank'`, `correl_norm_type='rank'` it reduces to
`rank**alpha`, matching GSVA — but its defaults differ (`min_size` 15 vs 10,
`max_size` 2000 vs 500) and it does no global range normalization. Decide
after the oracle comparison shows how each tracks GSVA.

## Suggested output schema

Long format, one row per (aliquot, pathway) — matches how
`gene_expression_quantification` already works in the tabular layout, keeps
Hallmark at 11,505 × 50 ≈ 575k rows, and extends to Reactome's ~1,700 sets
with no schema change. Wide format would need a new schema per collection.

```
case_id, case_submitter_id, sample_id, sample_submitter_id, aliquot_id,
sample_type, pathway, score_raw, score_normalized
```

with collection-level provenance (msigdb version + GMT md5, alpha, min/max
size, gene universe, divisor) on the dataset card.

MSigDB **2026.1.Hs** confirmed current. Hallmark GMT pinned:
`h.all.v2026.1.Hs.symbols.gmt`, md5 `367eec875967c2cfbf664a1a065b7b8d`,
50 sets, sizes 32–200 (so `maxSize=500` never binds for Hallmark).

## Remaining open questions

1. ~1,000 cases have multiple RNA aliquots — confirm per-aliquot is the
   published unit, with sample_type carried for tumor/normal filtering.
2. Whether to publish `score_normalized` at all, given the composition
   dependence, or only `score_raw` + divisor.
3. Whether v1 ships Hallmark alone or Hallmark + Reactome + WikiPathways
   together. Note this is not purely additive: if scores are normalized per
   run, the collection composition changes the divisor and therefore every
   Hallmark value. Shipping `score_raw` + divisor makes the question moot.

## Promoting the proof-of-concept

Working code (matrix build, ssGSEA, validation, oracle) currently lives in
a scratch dir. To make it real:

- `dev_research/ssgsea/` — the oracle fixture (`gsva_raw.csv`), `oracle.R`,
  `run_oracle.sh`, and the comparison script, following the
  `dev_research/liu_2018/` pattern.
- `ds2hf_pipeline/tcga/ssgsea.py` — matrix build + scoring.
- `tests/test_ssgsea.py` — canary against a handful of committed oracle
  values, plus the log1p-invariance and chunk-equivalence properties (all
  cheap, all deterministic).
- New tabular table `ssgsea_hallmark`; the `--table` append path means it
  costs one scoped build rather than a full re-derive.

[1]: https://github.com/rcastelo/GSVA/blob/devel/R/ssgsea.R?utm_source=chatgpt.com "GSVA/R/ssgsea.R at devel · rcastelo/GSVA · GitHub"
[2]: https://bioconductor.org/packages//release/bioc/vignettes/GSVA/inst/doc/GSVA.html?utm_source=chatgpt.com "GSVA: gene set variation analysis"
[3]: https://bioconductor.org/packages//release/bioc/vignettes/GSVA/inst/doc/GSVA.html "GSVA: gene set variation analysis"
[4]: https://www.gsea-msigdb.org/gsea/msigdb?utm_source=chatgpt.com "GSEA | MSigDB"
[5]: https://bioconductor.org/packages//release/bioc/html/GSVA.html?utm_source=chatgpt.com "Bioconductor - GSVA"
[6]: https://www.gsea-msigdb.org/gsea/msigdb_license_terms.jsp?utm_source=chatgpt.com "GSEA | License Terms for MSigDB released after April 2017"

---

## Gene universe: protein-coding + functional IG/TR segments (2026-08-12)

Revised after measuring MSigDB coverage. The universe is **protein-coding
plus the 409 functional immunoglobulin / T-cell-receptor segments**
(`IG_V/C/D/J_gene`, `TR_V/C/D/J_gene`), `_PAR_Y` dropped — 20,347 symbols.

**Why.** GENCODE gives Ig/TCR segments their own biotypes, so a plain
protein-coding filter drops all of them, and MSigDB's B-cell-receptor,
Fc-receptor and complement pathways are largely built from them. The tumour
immune microenvironment is central to oncology, so under-representing it is
a real limitation rather than a cosmetic one.

| pathway | before | after |
|---|---|---|
| REACTOME_CD22_MEDIATED_BCR_REGULATION | 5/61 (8%) | 61/61 (100%) |
| REACTOME_FCGR_ACTIVATION | 12/69 (17%) | 69/69 (100%) |
| REACTOME_SCAVENGING_OF_HEME_FROM_PLASMA | 13/69 (19%) | 69/69 (100%) |
| REACTOME_INITIAL_TRIGGERING_OF_COMPLEMENT | 22/80 (28%) | 79/80 (99%) |

Cost: existing scores move a median of 0.21 of a pathway's cross-sample SD
(p95 0.54, per-pathway Spearman >= 0.986). Pseudogene biotypes (~237) are
excluded — near-uniformly zero, they would only pad the low ranks.

These genes are genuinely measured by this assay: IG_C_gene is 71% above
1 TPM versus 52% for protein-coding, because B cells transcribe
polyadenylated Ig mRNA. IGKC ranks ~400 of 60,660 in an infiltrated tumour.

**Caveat to document on the card:** V/D/J segments are somatically
rearranged and clonally expanded, so their bulk expression reports
lymphocyte infiltration and clonality rather than regulation of a fixed
locus. That is what makes them informative about the microenvironment, but
it is a different kind of measurement from the rest of the matrix.

## Coverage by collection, against that universe

| collection | sets | mean match fraction | sets <90% covered |
|---|---|---|---|
| Hallmark | 50 | 0.9964 | 0 |
| Reactome | 1,839 | 0.9935 | 24 |
| WikiPathways | 925 | 0.9776 | 49 |

**WikiPathways is the weakest, and its gap is not fixable from this assay.**
Its poorly-covered sets are miRNA-centric (`WP_MIRNAS_INVOLVED_IN_DNA_DAMAGE_RESPONSE`
15/68, `WP_CELL_DIFFERENTIATION_INDEX` 14/50), and TCGA RNA-Seq is poly-A
selected — mature miRNAs are not polyadenylated, so STAR gene counts leave
them **94.8% exactly zero**. Adding the 1,881 miRNA rows would pad the low
ranks and inform nothing; adding `Mt_rRNA` would be worse still, dominating
the top ranks at ~14,600 TPM.

The right instrument is TCGA's separate miRNA-Seq assay (miRNA Expression
Quantification, 11,441 files) — a future modality, not a gene-universe fix.
Until then `matched_gene_count` / `original_gene_count` flag the affected
sets on every row.

## No maximum gene-set size

`maxSize=500` was inherited convention (Broad GSEA desktop uses 15/500;
gseapy 15/2000), **not a GSVA default — GSVA ships minSize=1, maxSize=Inf**.
The permutation-testing rationale doesn't transfer: ssGSEA runs no
permutation test and we publish descriptive scores. Measured leave-one-gene-out
sensitivity declines smoothly with size (97.9% of a score's cross-sample SD
at n=5, 62% at n=10, 17% at n=100) with no natural breakpoint, and large
sets are the *most* stable. So: no maximum. A minimum is still defensible
on stability grounds; the exact floor is a judgement call, and the
per-pathway gene counts we ship let consumers set their own.

## Which collections to ship, measured (2026-08-12)

Coverage of every plausible MSigDB 2026.1.Hs collection against our
20,347-symbol universe. `mean match` is the average fraction of a set's
genes present; `<90%` is the share of sets below 90% covered.

| collection | sets | mean match | <90% | verdict |
|---|---|---|---|---|
| C2:CP:PID | 196 | 0.9969 | 0.5% | **include** — NCI-Nature cancer signalling, best coverage measured |
| C2:CP:BIOCARTA | 292 | 0.9922 | 3.1% | optional — clean but older curation, small sets |
| C4:CM cancer modules | 431 | 0.9904 | 0.0% | optional — oncology-specific, 2005-era method |
| C3:MIR miRNA targets | 2,377 | 0.9857 | 0.6% | optional — these are *targets* (mRNAs), so unaffected by the miRNA problem |
| C4:CGN neighbourhoods | 427 | 0.9852 | 0.2% | optional |
| C4:3CA cancer cell atlas | 148 | 0.9780 | 4.7% | **include** — scRNA-derived tumour-cell meta-programs |
| C7 ImmuneSigDB | 4,872 | 0.9744 | 3.8% | defer — largely mouse vaccine/perturbation contrasts |
| C6 oncogenic signatures | 189 | 0.9723 | 1.6% | **include** — the one collection built for oncogene/TSG readout |
| C5 GO:BP | 7,538 | 0.9637 | 10.2% | defer — generic and large |
| C8 cell-type signatures | 866 | 0.9110 | 28.3% | maybe — relevant to TIME but lncRNA-heavy markers |
| C3:TFT TF targets | 506 | 0.7672 | 97.2% | **no** — needs ~16.9k lncRNA in the universe |
| C1 positional/cytoband | 302 | 0.3926 | 100.0% | **no** — cytobands enumerate pseudogenes; CNV is the right instrument |
| C2:CP:KEGG legacy | 186 | — | — | **no** — restricted licence |

**Proposed v2: Reactome + C6 + C4:3CA + C2:CP:PID.** A broad canonical
backbone plus the three oncology-specific collections, all above 0.97
coverage. 533 sets beyond Reactome, ~6.1M rows, ~75 MB.

**WikiPathways is deferred until TCGA miRNA-Seq is in the dataset.** Its
weak sets are miRNA-centric and this assay cannot measure them, so shipping
it now would publish 49 sets scoring something other than what their names
claim. That is a modality gap, not a curation gap: the fix is the separate
miRNA Expression Quantification assay (11,441 files), after which
WikiPathways coverage should be re-measured rather than assumed.

Note the general lesson: a collection's usability here is decided by whether
its genes live in the assay's universe, not by how good the curation is.
C3:TFT is well curated and unusable; C1 is fine and unusable. Measure
coverage before adding any collection.
