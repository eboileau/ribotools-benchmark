#!/usr/bin/env Rscript --vanilla

source(file.path(Sys.getenv("SCRIPT_DIR"), "utils.R"))
opts <- parse_args()

library(xtail)

ribo_counts <- as.matrix(read.csv(opts$ribo, row.names = 1, check.names = FALSE))
rna_counts  <- as.matrix(read.csv(opts$rna, row.names = 1, check.names = FALSE))
samples <- read.csv(opts$samples)

# assign condition labels to the columns of the mRNA and RPF data
# pick one assay (both are ordered identically) or use 1:8
condition <- factor(samples$condition[samples$assay=="ribo"])

# fails here...
res <- xtail(ribo_counts, rna_counts, condition)
