#!/usr/bin/env Rscript --vanilla

source(file.path(Sys.getenv("SCRIPT_DIR"), "utils.R"))
opts <- parse_args()

library(yaml)
library(openxlsx)

# merge count tables and write the config
ribo <- as.matrix(read.csv(opts$ribo, row.names = 1, check.names = FALSE))
rna  <- as.matrix(read.csv(opts$rna, row.names = 1, check.names = FALSE))
merge <- cbind(ribo, rna)
output_dir <- dirname(opts$out)
count_file <-  file.path(output_dir, "counts.csv")
write.csv(merge, count_file, row.names = TRUE, quote = FALSE)

# default lexicographic ordering
samples <- read.csv(opts$samples)
conditions <- as.character(levels(as.factor((samples$condition))))
contrast_key <- paste0(conditions[2], "_vs_", conditions[1])
contrast_list <- list()
contrast_list[[contrast_key]] <- c(conditions[2], conditions[1])

cfg <- list(
  tea_data = output_dir,
  contrasts = contrast_list,
  sample_table = opts$samples,
  count_table = count_file
)
yaml_file <- file.path(output_dir, "config.yaml")
write_yaml(cfg, yaml_file)

# call Ribotools
method <- "deltaTE"
status <- system2(
  "run-tea",
  args = c("--method", method, "--alpha", as.numeric(opts$alpha), "--lfcThreshold", 0, "--delim", "CSV", yaml_file),
  stdout = FALSE,
  stderr = FALSE
)
if (status != 0) {
  stop("run_ribotools.R failed with exit status ", status)
}

# unify output
res <- read.xlsx(file.path(output_dir, method, contrast_key, paste0(contrast_key, ".xlsx")))
tmp <- read.table(file.path(output_dir, method, contrast_key, "buffered.txt"), stringsAsFactors = F)
buffered <- as.vector(t(as.matrix(tmp[-1])))
tmp <- read.table(file.path(output_dir, method, contrast_key, "exclusive.txt"), stringsAsFactors = F)
exclusive <- as.vector(t(as.matrix(tmp[-1])))
tmp <- read.table(file.path(output_dir, method, contrast_key, "forwarded.txt"), stringsAsFactors = F)
forwarded <- as.vector(t(as.matrix(tmp[-1])))
tmp <- read.table(file.path(output_dir, method, contrast_key, "intensified.txt"), stringsAsFactors = F)
intensified <- as.vector(t(as.matrix(tmp[-1])))

res$effect <- res$log2FC.dte
res$score <- -log10(res$padj.dte + 1e-300)

res$effect[res$id %in% buffered] <- res$log2FC.rna[res$id %in% buffered] 
res$score[res$id %in% buffered] <- -log10(res$padj.rna[res$id %in% buffered] + 1e-300)
res$effect[res$id %in% forwarded] <- res$log2FC.rna[res$id %in% forwarded] 
res$score[res$id %in% forwarded] <- -log10(res$padj.rna[res$id %in% forwarded] + 1e-300)

res$class <- "background"
res$class[res$id %in% forwarded] <- "abundance"
res$class[res$id %in% intensified] <- "abundance"
res$class[res$id %in% exclusive] <- "translation"
res$class[res$id %in% buffered] <- "buffering"

res$identifier <- res$id
res <- res[,c("identifier", "class", "effect", "score")]

write.csv(res, opts$out, row.names = FALSE, quote = FALSE)
