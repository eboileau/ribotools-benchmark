#!/usr/bin/env Rscript --vanilla

source(file.path(Sys.getenv("SCRIPT_DIR"), "utils.R"))
opts <- parse_args()

library(anota2seq)

ribo_counts <- as.matrix(read.csv(opts$ribo, row.names = 1, check.names = FALSE))
rna_counts  <- as.matrix(read.csv(opts$rna, row.names = 1, check.names = FALSE))
samples <- read.csv(opts$samples)

# phenoVec must describe the sample class for corresponding columns in dataT and dataP
# pick one assay (both are ordered identically) or use 1:8
condition <- factor(samples$condition[samples$assay=="ribo"])

# defaults - filterZeroGenes should have no effect on the simulated data
ads <- anota2seqDataSetFromMatrix(
  dataP = ribo_counts,
  dataT = rna_counts,
  phenoVec = condition,
  dataType = "RNAseq",
  normalize = TRUE,
  filterZeroGenes = FALSE,
)

# complete analysis using the one-step procedure function
# no "effect size filter"
ads <- anota2seqRun(
  Anota2seqDataSet = ads,
  thresholds = list(
    maxPAdj = 0.05,
    minEff = NULL
  ),
)

# regulatory modes: translatedmRNA, totalmRNA, translation, buffering
# p-values and group effect size (log2FC)
# for buffering: a positive `apvEff` means buffered down and a negative `apvEff` means buffered up
res <- anota2seqGetOutput(
  ads,
  output = "singleDf",
  selContrast = 1,
  getRVM = TRUE
)

# unify output
res$effect <- res$totalmRNA.apvEff
res$score <- -log10(res$totalmRNA.apvRvmPAdj + 1e-300)
res$class <- res$singleRegMode
res$effect[res$singleRegMode == "translation"] <-
    res$translation.apvEff[res$singleRegMode == "translation"]
res$score[res$singleRegMode == "translation"] <-
    -log10(res$translation.apvRvmPAdj[res$singleRegMode == "translation"] + 1e-300)
res$effect[res$singleRegMode == "buffering"] <-
    -res$buffering.apvEff[res$singleRegMode == "buffering"]
res$score[res$singleRegMode == "buffering"] <-
    -log10(res$buffering.apvRvmPAdj[res$singleRegMode == "buffering"] + 1e-300)
res <- res[,c("identifier", "class", "effect", "score")]

write.csv(res, opts$out, row.names = FALSE, quote = FALSE)
