#!/usr/bin/env Rscript --vanilla

source(file.path(Sys.getenv("SCRIPT_DIR"), "utils.R"))
opts <- parse_args_permute()

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

# calling 
run_ribotools <- function(seed) {
  set.seed(seed)
  perm_cols <- sample(ncol(ribo))
  ribo_perm <- ribo[, perm_cols]
  rna_perm  <- rna[, perm_cols]
  colnames(ribo_perm) <- colnames(ribo)
  colnames(rna_perm)  <- colnames(rna)
  # write data...
  merge <- cbind(ribo_perm, rna_perm)
  output_dir <- dirname(opts$out)
  count_file <-  file.path(output_dir, paste0("counts", seed, ".csv"))
  write.csv(merge, count_file, row.names = TRUE, quote = FALSE)
  # ... and config
  cfg <- list(
    tea_data = output_dir,
    contrasts = contrast_list,
    sample_table = opts$samples,
    count_table = count_file
  )
  yaml_file <- file.path(output_dir, paste0("config", seed, ".yaml"))
  write_yaml(cfg, yaml_file)
  # call Ribotools using default lfcThreshold
  status <- system2(
    "run-tea",
    args = c("--method", method, "--alpha", alpha, "--delim", "CSV", yaml_file),
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
    iteration = seed,
    n_sig_te = sum(res$padj.dte < alpha, na.rm = TRUE),
    n_sig_ribo = sum(res$padj.ribo < alpha, na.rm = TRUE),
    n_sig_rna = sum(res$padj.rna < alpha, na.rm = TRUE),
    n_tested = sum(!is.na(res)),
    stringsAsFactors = FALSE
  )
  # force clean...
  unlink(c(count_file, yaml_file), recursive = FALSE, force = FALSE, expand = FALSE)
  unlink(c(file.path(output_dir, method)), recursive = TRUE, force = FALSE, expand = FALSE)
  df
}

res <- lapply(seq_len(opts$permutations), run_ribotools)
out <- do.call(rbind, res)
write.csv(out, opts$out, row.names = FALSE, quote = FALSE)
