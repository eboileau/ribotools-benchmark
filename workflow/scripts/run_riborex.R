#!/usr/bin/env Rscript --vanilla

source(file.path(Sys.getenv("SCRIPT_DIR"), "utils.R"))
opts <- parse_args()

library(riborex)
cat("Riborex", file=opts$out)

# ribo_counts <- as.matrix(read.csv(opt$ribo,   row.names = 1, check.names = FALSE))
# rna_counts  <- as.matrix(read.csv(opt$rna,    row.names = 1, check.names = FALSE))
# labels      <- read.csv(opt$labels)
# 
# condition <- labels$condition   # expects a 'condition' column
# 
# res       <- xtail(ribo_counts, rna_counts, condition, bins = 1000)
# res_table <- resultsTable(res)
# 
# out <- data.frame(
#   gene_id    = rownames(res_table),
#   log2FC_TE  = res_table$log2FC_TE_v1,
#   log2FC_RPF = NA_real_,
#   log2FC_RNA = NA_real_,
#   padj_TE    = res_table$pvalue.adjust,
#   score      = -log10(res_table$pvalue.adjust + 1e-300),
#   tool       = "xtail",
#   stringsAsFactors = FALSE
# )
# 
# write.csv(out, opt$out, row.names = FALSE, quote = FALSE)
