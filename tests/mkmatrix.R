#!/usr/bin/env Rscript
# Generate a reproducible negative-binomial-ish gene x sample count matrix.
# Usage: mkmatrix.R <ngenes> <nsamples> <seed> <out.tsv>
args <- commandArgs(trailingOnly = TRUE)
ng <- as.integer(args[1]); ns <- as.integer(args[2]); seed <- as.integer(args[3]); out <- args[4]
set.seed(seed)
base <- rgamma(ng, shape = 0.6, scale = 60)
libfac <- runif(ns, 0.7, 1.4)
m <- matrix(0L, ng, ns)
for (j in seq_len(ns)) m[, j] <- rpois(ng, base * libfac[j])
rownames(m) <- sprintf("ENSG%08d", seq_len(ng))
colnames(m) <- sprintf("S%03d", seq_len(ns))
con <- file(out, "w")
writeLines(paste(c("gene", colnames(m)), collapse = "\t"), con)
for (i in seq_len(ng)) writeLines(paste(c(rownames(m)[i], m[i, ]), collapse = "\t"), con)
close(con)
