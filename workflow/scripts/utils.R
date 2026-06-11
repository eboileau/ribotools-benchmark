
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args)!=8) { stop("run_tool.R --ribo RIBO --rna RNA --samples SAMPLES --out OUT\n", call.=FALSE) }
  keys <- sub("^--", "", args[seq(1, length(args), 2)])
  vals <- args[seq(2, length(args), 2)]
  setNames(as.list(vals), keys)
}
