
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

msg_evaluation <- "run_evaluation.R --rna RNA --translation TRAN --abundance ABUN --buffering BUFF --results FILE1 FILE2 [FILE3, ...] --out OUT"
keys_evaluation.expected <- c("abundance", "buffering", "out", "results", "rna", "translation")

parse_args_evaluation <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  res <- list()
  i <- 1
  while (i <= length(args)) {
    arg <- args[i]
    if (!startsWith(arg, "--")) {
      stop("Unexpected argument: ", arg)
    }
    key <- sub("^--", "", arg)
    i <- i + 1
    values <- character()
    while (i <= length(args) && !startsWith(args[i], "--")) {
      values <- c(values, args[i])
      i <- i + 1
    }
    res[[key]] <- values
  }
  if (!(identical(sort(names(res)), keys_evaluation.expected))) { stop(msg_evaluation, call.=FALSE) }
  for (key in keys_evaluation.expected) {
    if (key!="results") {if (length(res[[key]])!=1) { stop(msg_evaluation, call.=FALSE) }}
  }
  if (length(res[["results"]])<2) { stop(msg_evaluation, call.=FALSE) }
  res
}
