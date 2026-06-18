#!/usr/bin/env Rscript --vanilla

source(file.path(Sys.getenv("SCRIPT_DIR"), "utils.R"))
opts <- parse_args()

library(riborex)

ribo <- as.matrix(read.csv(opts$ribo, row.names = 1, check.names = FALSE))
rna  <- as.matrix(read.csv(opts$rna, row.names = 1, check.names = FALSE))
samples <- read.csv(opts$samples)

# We need to prepare two vectors to indicate the treatments of samples in RNA- and Ribo-seq data.
condition <- factor(samples$condition[samples$assay=="ribo"])

# use default engine
res.deseq2 <- riborex(rna, ribo, condition, condition)
res <- data.frame(res.deseq2)

# unify output
res$effect <- res$log2FoldChange
res$score <- -log10(res$padj + 1e-300)
res$class <- "background"
res$class[res$padj<opts$alpha] <- "translation"
res$identifier <- rownames(res)
res <- res[,c("identifier", "class", "effect", "score")]

write.csv(res, opts$out, row.names = FALSE, quote = FALSE)
