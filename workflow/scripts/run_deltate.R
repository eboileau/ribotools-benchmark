#!/usr/bin/env Rscript --vanilla

source(file.path(Sys.getenv("SCRIPT_DIR"), "utils.R"))
opts <- parse_args()

library(DESeq2)
cat("deltaTE", file=opts$out)

# ribo_counts <- as.matrix(read.csv(opt$ribo,   row.names = 1, check.names = FALSE))
# rna_counts  <- as.matrix(read.csv(opt$rna,    row.names = 1, check.names = FALSE))
# labels      <- read.csv(opt$samples)
# 
# # Build combined count matrix and coldata for interaction model
# counts_combined <- cbind(ribo_counts, rna_counts)
# 
# col_data <- data.frame(
#   condition   = factor(rep(labels$condition, 2)),
#   libtype     = factor(c(rep("RPF", ncol(ribo_counts)),
#                          rep("RNA", ncol(rna_counts)))),
#   row.names   = colnames(counts_combined)
# )
# col_data$condition <- relevel(col_data$condition, ref = levels(col_data$condition)[1])
# 
# # Three DESeq2 models
# dds_te <- DESeqDataSetFromMatrix(counts_combined, col_data,
#                                   design = ~ condition + libtype + condition:libtype)
# dds_ribo <- DESeqDataSetFromMatrix(ribo_counts,
#                                     data.frame(condition = factor(labels$condition),
#                                                row.names = colnames(ribo_counts)),
#                                     design = ~ condition)
# dds_rna  <- DESeqDataSetFromMatrix(rna_counts,
#                                     data.frame(condition = factor(labels$condition),
#                                                row.names = colnames(rna_counts)),
#                                     design = ~ condition)
# 
# dds_te   <- DESeq(dds_te)
# dds_ribo <- DESeq(dds_ribo)
# dds_rna  <- DESeq(dds_rna)
# 
# lvls    <- levels(col_data$condition)
# ref     <- lvls[1]
# trt     <- lvls[2]
# 
# res_te   <- results(dds_te,   name = paste0("condition", trt, ".libtypeRPF"))
# res_ribo <- results(dds_ribo, contrast = c("condition", trt, ref))
# res_rna  <- results(dds_rna,  contrast = c("condition", trt, ref))
# 
# out <- data.frame(
#   gene_id    = rownames(res_te),
#   log2FC_TE  = res_te$log2FoldChange,
#   log2FC_RPF = res_ribo$log2FoldChange,
#   log2FC_RNA = res_rna$log2FoldChange,
#   padj_TE    = res_te$padj,
#   score      = -log10(res_te$padj + 1e-300),
#   tool       = "deltate",
#   stringsAsFactors = FALSE
# )
# 
# # Replace NA padj (filtered genes) with worst score
# out$padj_TE[is.na(out$padj_TE)] <- 1
# out$score  <- -log10(out$padj_TE + 1e-300)
# 
# write.csv(out, opt$out, row.names = FALSE, quote = FALSE)
