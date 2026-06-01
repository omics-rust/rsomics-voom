#!/usr/bin/env Rscript
# voom oracle: read a gene x sample count matrix TSV (header = sample ids,
# first column = gene ids), run limma::voom with the default intercept-only
# design, and write the E matrix and the weights matrix as TSVs.
#
# Usage: voom_oracle.R <counts.tsv> <E_out.tsv> <weights_out.tsv>
suppressMessages(library(limma))

args <- commandArgs(trailingOnly = TRUE)
counts_path <- args[1]
e_out <- args[2]
w_out <- args[3]

m <- as.matrix(read.delim(counts_path, row.names = 1, check.names = FALSE))
v <- voom(m)
dimnames(v$E) <- dimnames(m)
dimnames(v$weights) <- dimnames(m)

write_tsv <- function(mat, path) {
  con <- file(path, "w")
  writeLines(paste(c("gene", colnames(mat)), collapse = "\t"), con)
  for (i in seq_len(nrow(mat))) {
    vals <- formatC(mat[i, ], digits = 10, format = "g", flag = "")
    writeLines(paste(c(rownames(mat)[i], trimws(vals)), collapse = "\t"), con)
  }
  close(con)
}

write_tsv(v$E, e_out)
write_tsv(v$weights, w_out)
