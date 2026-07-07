# rsomics-voom

voom log2-CPM transform with mean-variance precision weights for RNA-seq count
matrices. A single-binary Rust reimplementation of limma's `voom()`.

Given a gene x sample count matrix it produces the two matrices a downstream
weighted linear model needs:

- **E** — log2 counts-per-million, `log2((count + 0.5) / (lib.size + 1) * 1e6)`.
- **weights** — per-observation precision weights from the mean-variance trend,
  `1 / predicted_sd(fitted_log2_count)^4`.

```
rsomics-voom counts.tsv -o E.tsv -w weights.tsv
rsomics-voom counts.tsv -w weights.tsv > E.tsv   # E to stdout
```

Input is a TSV with a header row of sample IDs and gene IDs in the first
column. The design is intercept-only (one mean per gene), matching
`voom(counts)` with no explicit design.

## Method

1. Library sizes are column sums; `E = log2((count + 0.5)/(lib.size + 1) * 1e6)`.
2. A gene-wise intercept-only fit gives the mean log2-CPM (`Amean`) and residual
   standard deviation (`sigma`, residual df `= n - 1`).
3. The mean-variance trend is `sqrt(sigma)` against
   `Amean + mean(log2(lib.size + 1)) - log2(1e6)`, with all-zero-count genes
   dropped, smoothed by Cleveland LOWESS. The span defaults to `0.5`
   (`voom()`'s default); `--span <f>` sets a fixed span, and `--adaptive-span`
   selects it from the gene count as `min(0.3 + 0.7 * (50/ngenes)^(1/3), 1)`
   (`voom(adaptive.span = TRUE)`), overriding `--span`.
4. Each observation's fitted log2-count is read off the smoothed trend (linear
   interpolation, endpoint-clamped) to give its predicted standard deviation;
   the weight is the inverse fourth power.

## Origin

This crate is an independent Rust reimplementation of limma's `voom` based on:

- The published method (Law, Chen, Shi & Smyth, "voom: precision weights unlock
  linear model analysis tools for RNA-seq read counts", Genome Biology 2014,
  15:R29, DOI 10.1186/gb-2014-15-2-r29).
- The voom mean-variance modelling description and the public LOWESS algorithm
  (W. S. Cleveland, "Robust Locally Weighted Regression and Smoothing
  Scatterplots", JASA 1979, 74:829-836).
- Black-box behaviour testing against the `voom()` binary output.

No source code from the GPL upstream was used as reference during
implementation. Test fixtures are independently generated.

Output is verified value-exact (relative deviation < 1e-6) against limma 3.62.1
voom across several matrix shapes: the default matches `voom(counts, design)`
(fixed span 0.5) and `--adaptive-span` matches `voom(counts, design,
adaptive.span = TRUE)`. The committed goldens in `tests/golden/` — including a
constant-gene degenerate matrix where the two spans diverge 73% — let CI
validate without an R install.

License: MIT OR Apache-2.0.
Upstream credit: limma <https://bioconductor.org/packages/limma/> (GPL >= 2).
