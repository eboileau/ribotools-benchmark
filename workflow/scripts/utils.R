
msg <- "run_tool.R --ribo RIBO --rna RNA --samples SAMPLES --out OUT --alpha ALPHA\n"
keys.expected <- c("alpha", "out", "ribo", "rna", "samples")

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args)!=10) { stop(msg, call.=FALSE) }
  keys <- sub("^--", "", args[seq(1, length(args), 2)])
  if (!(identical(sort(keys), keys.expected))) { stop(msg, call.=FALSE) }
  vals <- args[seq(2, length(args), 2)]
  setNames(as.list(vals), keys)
}
