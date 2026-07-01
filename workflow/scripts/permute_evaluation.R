#!/usr/bin/env Rscript --vanilla

# Benchmarking:
# permutation calibration test or "power characterisation"

source(file.path(Sys.getenv("SCRIPT_DIR"), "utils.R"))
opts <- parse_args_permute()

df <- do.call(rbind, lapply(opts$results, read.csv, stringsAsFactors = FALSE))

df$fpr_te <- df$n_sig_te / df$n_tested_te
df$fpr_ribo <- df$n_sig_ribo / df$n_tested_ribo
df$fpr_rna <- df$n_sig_rna / df$n_tested_rna

summarize <- function(df, col) {
  agg_mean <- aggregate(df[[col]], by = list(tool = df$tool), FUN = mean, na.rm = TRUE)
  agg_sd <- aggregate(df[[col]], by = list(tool = df$tool), FUN = sd,   na.rm = TRUE)
  agg_n <- aggregate(df[[col]], by = list(tool = df$tool), FUN = function(x) sum(!is.na(x)))
  data.frame(
    tool = agg_mean$tool,
    mean = agg_mean$x,
    sd = agg_sd$x,
    n = agg_n$x
  )
}

te_summary <- summarize(df, "fpr_te") 
ribo_summary <- summarize(df, "fpr_ribo")
rna_summary <- summarize(df, "fpr_rna")
names(te_summary)[2:4] <- paste0("te_",   names(te_summary)[2:4])
names(ribo_summary)[2:4] <- paste0("ribo_", names(ribo_summary)[2:4])
names(rna_summary)[2:4] <- paste0("rna_",  names(rna_summary)[2:4])

summary_table <- Reduce(function(x, y) merge(x, y, by = "tool", all = TRUE),
                         list(te_summary, ribo_summary, rna_summary))

write.csv(summary_table, opts$out, row.names = FALSE, quote = FALSE)
