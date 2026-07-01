#!/usr/bin/env Rscript --vanilla

source(file.path(Sys.getenv("SCRIPT_DIR"), "utils.R"))
opts <- parse_args()

library(yaml)
library(openxlsx)

method <- "LRT"

# read count and sample tables
ribo <- as.matrix(read.csv(opts$ribo, row.names = 1, check.names = FALSE))
rna  <- as.matrix(read.csv(opts$rna, row.names = 1, check.names = FALSE))

# default lexicographic ordering
samples <- read.csv(opts$samples)
conditions <- as.character(levels(as.factor((samples$condition))))
contrast_key <- paste0(conditions[2], "_vs_", conditions[1])
contrast_list <- list()
contrast_list[[contrast_key]] <- c(conditions[2], conditions[1])

alpha <- as.numeric(opts$alpha)

# define allowed permutations 
enumerate_col_perms <- function(samples) {
  # one assay to define the valid permutations
  condition <- as.character(samples$condition[samples$assay=="ribo"])
  lvls <- levels(as.factor(condition))
  n <- length(condition)
  n_ctrl <- sum(condition == lvls[1])
  # which column indices get "control"
  all_combn <- combn(n, n_ctrl, simplify = FALSE)
  # remove original assignment
  valid <- Filter(function(x) !setequal(x, which(condition == lvls[1])), all_combn)
  valid <- Filter(function(x) !setequal(x, which(condition == lvls[2])), valid)  
  # full column permutation orders
  lapply(valid, function(idx) {
    c(idx, setdiff(seq_len(n), idx))
  })
}

# calling 
run_ribotools <- function(perm_cols, iter) {
  ribo_perm <- ribo[, perm_cols]
  rna_perm  <- rna[, perm_cols]
  colnames(ribo_perm) <- colnames(ribo)
  colnames(rna_perm)  <- colnames(rna)
  # write data...
  merge <- cbind(ribo_perm, rna_perm)
  output_dir <- dirname(opts$out)
  count_file <-  file.path(output_dir, paste0("counts", iter, ".csv"))
  write.csv(merge, count_file, row.names = TRUE, quote = FALSE)
  # ... and config
  cfg <- list(
    tea_data = output_dir,
    contrasts = contrast_list,
    sample_table = opts$samples,
    count_table = count_file
  )
  yaml_file <- file.path(output_dir, paste0("config", iter, ".yaml"))
  write_yaml(cfg, yaml_file)
  # call Ribotools using default lfcThreshold
  status <- system2(
    "run-tea",
    args = c("--method", method, "--alpha", alpha, "--lfcThreshold", 0, "--delim", "CSV", yaml_file),
    stdout = FALSE,
    stderr = FALSE
  )
  if (status != 0) {
    stop("run_ribotools.R failed with exit status ", status)
  }
  # evaluate 
  res <- read.xlsx(file.path(output_dir, method, contrast_key, paste0(contrast_key, ".xlsx")))
  df <- data.frame(
    tool = "ribotools_lrt",
    iteration = iter,
    perm_cols = paste(perm_cols, collapse=";"),
    n_sig_te = sum(res$padj.dte < alpha, na.rm = TRUE),
    n_sig_ribo = sum(res$padj.ribo < alpha, na.rm = TRUE),
    n_sig_rna = sum(res$padj.rna < alpha, na.rm = TRUE),
    n_tested_te = sum(!is.na(res$pvalue.dte)),
    n_tested_ribo = sum(!is.na(res$pvalue.ribo)),
    n_tested_rna = sum(!is.na(res$pvalue.rna)),
    stringsAsFactors = FALSE
  )
  # force clean...
  unlink(c(count_file, yaml_file), recursive = FALSE, force = FALSE, expand = FALSE)
  unlink(c(file.path(output_dir, method)), recursive = TRUE, force = FALSE, expand = FALSE)
  df
}

valid_perms <- enumerate_col_perms(samples)
res <- lapply(seq_along(valid_perms), function(idx) {
  run_ribotools(perm_cols=valid_perms[[idx]], iter=idx)
})
out <- do.call(rbind, res)
write.csv(out, opts$out, row.names = FALSE, quote = FALSE)
