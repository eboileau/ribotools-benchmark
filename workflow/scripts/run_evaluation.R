#!/usr/bin/env Rscript --vanilla

# Benchmarking: evaluate performance across 4 test cases (see below for description)

# Input - output from run_tool.R
# ------------------------------
# identifier: gene ID
# class: predicted class - terminology is fixed
# effect: log2FC or equivalent effect size
# score: -log10(adjpval)

source(file.path(Sys.getenv("SCRIPT_DIR"), "utils.R"))
opts <- parse_args_evaluation()

suppressPackageStartupMessages({
  library(pROC)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(RColorBrewer)
})

# extract tool names from files
tool_files <- opts$results |>
  as.list() |>
  set_names(~ basename(dirname(.x)))

# exclude tools that do not output classes for classification accuracy
tool_to_exclude <- c("riborex")
class_tools <- setdiff(names(tool_files), tool_to_exclude)
class_levels <- c("translation", "abundance", "buffering", "background")

# test cases - self-explanatory
test_cases <- list(
  translation_strict = list(
    label = "Translation (strict)",
    tp_cls = "translation",
    tn_cls = c("background", "abundance", "buffering"),
  ),
  translation_only = list(
    label = "Translation (only)",
    tp_cls = "translation",
    tn_cls = "background",
  ),
  any_regulated = list(
    label = "Any regulated",
    tp_cls = c("translation", "abundance", "buffering"),
    tn_cls = "background",
  ),
  buffering = list(
    label = "Buffering",
    tp_cls = "buffering",
    tn_cls = c("background", "translation", "abundance"),
  )
)

# reference gene list - ribo and rna must have the same gene set
reference_genes  <- rownames(as.matrix(read.csv(opts$rna, row.names = 1, check.names = FALSE)))
n_ref <- length(reference_genes)

# ground truth
read_gt <- function(path) {
  g <- trimws(readLines(path, warn = FALSE))
  g[nchar(g) > 0]
}
gt_translation <- read_gt(opts$translation)
gt_abundance   <- read_gt(opts$abundance)
gt_buffering   <- read_gt(opts$buffering)

# read results
assign_gt_class <- function(id) {
  case_when(
    id %in% gt_translation ~ "translation",
    id %in% gt_abundance   ~ "abundance",
    id %in% gt_buffering   ~ "buffering",
    TRUE                   ~ "background"
  )
}

read_results <- function(path, tool) {
  df <- read.csv(path)
  n_reported <- nrow(df)
  n_in_ref <- sum(df$identifier %in% reference_genes)
  missing_genes <- setdiff(reference_genes, df$identifier)
  if (length(missing_genes) > 0) {
    pad_df <- tibble(
      identifier = missing_genes,
      class = "background",
      effect = 0,
      score = 0  # -log10(1) = 0
    )
    df <- bind_rows(df, pad_df)
  }
  df <- df %>%
    mutate(
      tool = tool,
      gt_class = assign_gt_class(identifier),
      patched = identifier %in% missing_genes
    )
  list(
    data = df,
    stats = tibble(
      tool = tool,
      n_ref = n_ref,
      n_reported = n_reported,
      n_in_ref = n_in_ref,
      n_padded = length(missing_genes)
    )
  )
}

results <- imap(tool_files, read_results)
data <- map(results, "data")
stats <- map_dfr(results, "stats")

write.csv(stats, file.path(opts$out, "stats.csv"), row.names = FALSE, quote = FALSE)

stopifnot(length(unique(map_int(data, nrow))) == 1)

## A. ROC + AUC - all tools
# does a tool's score rank true regulated genes above background? (no class)

compute_roc <- function(df, tc) {
  eval_df <- df %>%
    filter(gt_class %in% c(tc$tp_cls, tc$tn_cls)) %>%
    mutate(binary_label = as.integer(gt_class %in% tc$tp_cls))
  n_pos <- sum(eval_df$binary_label)
  n_neg <- sum(!eval_df$binary_label)
  roc_obj <- roc(
    response = eval_df$binary_label,
    predictor = eval_df$score,
    direction = "<",
    quiet = TRUE
  )
  auc_val <- as.numeric(auc(roc_obj))
  ci_obj  <- ci.auc(roc_obj,
    method = "delong",
    boot.n = 2000,
    conf.level = 1 - 0.05
  )
  list(
    roc_obj = roc_obj,
    auc = auc_val,
    ci_lower = ci_obj[1],
    ci_upper = ci_obj[3],
    n_pos = n_pos,
    n_neg = n_neg,
    tc_label = tc$label
  )
}

evaluation <- map(data, function(df) {
  map(test_cases, function(tc) compute_roc(df, tc))
})

final_results <- imap_dfr(evaluation, function(tc_list, tool_nm) {
  imap_dfr(tc_list, function(res, tc_nm) {
    tibble(
      tool = tool_nm,
      test_case = tc_nm,
      test_label = res$tc_label,
      auc = round(res$auc, 4),
      ci_lower = round(res$ci_lower, 4),
      ci_upper = round(res$ci_upper, 4),
      n_pos = res$n_pos,
      n_neg = res$n_neg
    )
  })
})

write.csv(final_results, file.path(opts$out, "final_results.csv"), row.names = FALSE, quote = FALSE)

roc_df <- imap_dfr(evaluation, function(tc_list, tool_nm) {
  imap_dfr(tc_list, function(res, tc_nm) {
    if (is.null(res)) return(NULL)
    ro <- res$roc_obj
    tibble(
      tool = tool_nm,
      test_case  = tc_nm,
      test_label = res$tc_label,
      fpr = 1 - ro$specificities,
      tpr = ro$sensitivities
    )
  })
})

auc_labels <- final_results %>%
  filter(!is.na(auc)) %>%
  mutate(legend_label = sprintf(
    "%s  (AUC\u202f=\u202f%.3f\u202f[%.3f\u2013%.3f])",
    tool, auc, ci_lower, ci_upper
))

tc_order <- map_chr(test_cases, "label")
roc_df <- roc_df %>%
  left_join(auc_labels %>% select(tool, test_case, legend_label),
            by = c("tool", "test_case")) %>%
  mutate(test_label = factor(test_label, levels = tc_order))
  
# colors... adjust for more tools.
n_tools <- length(tool_files)
  tool_colours <- setNames(
  RColorBrewer::brewer.pal(max(3, n_tools), "Dark2")[seq_len(n_tools)],
  names(tool_files)
)
legend_colour_map <- auc_labels %>%
  select(tool, legend_label) %>%
  distinct() %>%
  mutate(colour = tool_colours[tool]) %>%
  { setNames(.$colour, .$legend_label) }

# plot ROC
# Error bars: 95% CI (delong), gene universe: full set genes (missing padded with score=0)
# score: -log10(adjpval)
p_roc <- ggplot(roc_df, aes(x = fpr, y = tpr,
                            colour = legend_label,
                            group = tool)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_line(linewidth = 0.85, alpha = 0.9) +
  facet_wrap(~ test_label, nrow = 2, ncol = 2) +
  scale_colour_manual(values = legend_colour_map, name = NULL) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     expand = expansion(mult = 0.02)) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     expand = expansion(mult = 0.02)) +
  labs(
    x = "1 - Specificity",
    y = "Sensitivity"
  ) +
  coord_equal() +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    legend.key.width = unit(1.5, "cm"),
    strip.background = element_rect(fill = "grey92"),
    strip.text = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  ) +
  guides(colour = guide_legend(ncol = 1))
  
# ggsave("roc.pdf", p_roc, width = 12, height = 9, device = "pdf")
  
# plot AUC
# Error bars: 95% CI (delong), gene universe: full set genes (missing padded with score=0)
bar_df <- final_results %>%
  filter(!is.na(auc)) %>%
  mutate(
    test_label = factor(test_label, levels = tc_order),
    tool = factor(tool, levels = names(tool_files))
  )

p_bar <- ggplot(bar_df,
                aes(x = tool, y = auc,
                    fill = tool, ymin = ci_lower, ymax = ci_upper)) +
  geom_col(width = 0.65, alpha = 0.85) +
  geom_errorbar(width = 0.25, linewidth = 0.6, colour = "grey30") +
  geom_hline(yintercept = 0.5, linetype = "dashed",
             colour = "grey40", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.3f", auc), y = ci_upper),
            vjust = -0.4, size = 2.8, colour = "grey20") +
  facet_wrap(~ test_label, nrow = 2, ncol = 2) +
  scale_fill_manual(values = tool_colours, guide = "none") +
  scale_y_continuous(
    limits = c(0, 1.10),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = NULL,
    y = "Area under ROC curve (AUC)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, size = 9),
    strip.background = element_rect(fill = "grey92"),
    strip.text = element_text(face = "bold", size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )
  
# ggsave("bars.pdf", p_bar, width = 12, height = 9, device = "pdf")
   
## B. CLASSIFICATION - only tools that output a class label
# per-class precision/recall/F1/MCC

# binary MCC for one class vs rest
binary_mcc <- function(tp, fp, tn, fn) {
  denom <- sqrt(as.double((tp + fp)) * as.double((tp + fn)) * as.double((tn + fp)) * as.double((tn + fn)))
  if (denom == 0) return(NA_real_)
    (tp * tn - fp * fn) / denom
}
  
class_metrics <- function(pred, truth, levels) {
  # coerce to factor with identical levels so table is always square
  pred <- factor(pred, levels = levels)
  truth <- factor(truth, levels = levels)
  cm <- table(Predicted = pred, Truth = truth)
  map_dfr(levels, function(cls) {
    tp <- cm[cls, cls]
    fp <- sum(cm[cls, ]) - tp # predicted cls, truth != cls
    fn <- sum(cm[, cls]) - tp # truth cls, predicted != cls
    tn <- sum(cm) - tp - fp - fn
    precision <- if ((tp + fp) == 0) NA_real_ else tp / (tp + fp)
    recall <- if ((tp + fn) == 0) NA_real_ else tp / (tp + fn)
    f1 <- if (is.na(precision) || is.na(recall) ||
              (precision + recall) == 0) NA_real_ else
              2 * precision * recall / (precision + recall)
    mcc <- binary_mcc(tp, fp, tn, fn)
    tibble(
      class = cls,
      n_truth = as.integer(tp + fn), # total true positives in ground truth
      n_pred = as.integer(tp + fp), # total predicted as this class
      TP = as.integer(tp), FP = as.integer(fp),
      TN = as.integer(tn), FN = as.integer(fn),
      precision = round(precision, 4),
      recall = round(recall, 4),
      F1 = round(f1, 4),
      MCC = round(mcc, 4)
    )
  })
}
  
cls_metrics_df <- map_dfr(class_tools, function(tool_nm) {
  df <- data[[tool_nm]]
  # map class values to levels - names must match
  pred <- df$class
  truth <- df$gt_class
  metrics <- class_metrics(pred, truth, class_levels)
  metrics$tool <- tool_nm
  metrics
}) %>%
  select(tool, class, n_truth, n_pred, TP, FP, TN, FN,
        precision, recall, F1, MCC)
           
write.csv(stats, file.path(opts$out, "classification_metrics.csv"), row.names = FALSE, quote = FALSE)

# plot confution matrices
# Cell = count (% of column / truth class). Diagonal = correct calls
cm_long <- map_dfr(class_tools, function(tool_nm) {
    df    <- data[[tool_nm]]
    pred  <- factor(df$class, levels = class_levels)
    truth <- factor(df$gt_class,   levels = class_levels)
    cm    <- as.data.frame(table(Predicted = pred, Truth = truth))
    cm$tool <- tool_nm
    cm
  }) %>%
    group_by(tool, Truth) %>%
    mutate(pct = 100 * Freq / sum(Freq)) %>%   # % of truth class predicted as X
    ungroup() %>%
    mutate(
      tool      = factor(tool, levels = class_tools),
      Predicted = factor(Predicted, levels = rev(class_levels)),  # y-axis top-down
      Truth     = factor(Truth,     levels = class_levels)
    )

p_cm <- ggplot(cm_long, aes(x = Truth, y = Predicted, fill = pct)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%d\n(%.0f%%)", Freq, pct)),
              size = 2.8, lineheight = 1.1) +
  facet_wrap(~ tool, ncol = length(class_tools)) +
  scale_fill_gradientn(
    colours = c("white", "#FFF3B0", "#E09F3E", "#9E2A2B"),
    limits  = c(0, 100),
    name    = "% of\ntruth class"
  ) +
  labs(
    x = "Ground truth class",
    y = "Tool predicted class"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    strip.background = element_rect(fill = "grey92"),
    strip.text = element_text(face = "bold", size = 10),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold")
  )

n_ct <- length(class_tools)
# ggsave("confusion_matrices.pdf", p_cm,
#          width  = max(6, 4 * n_ct),
#          height = 5,
#          device = "pdf")
         
# plot bars F1 and MCC
# F1: harmonic mean of precision & recall
# MCC: Matthews Correlation Coefficient (robust to class imbalance) - e.g. background
metrics_long <- cls_metrics_df %>%
  select(tool, class, F1, MCC) %>%
  pivot_longer(c(F1, MCC), names_to = "metric", values_to = "value") %>%
  mutate(
    tool = factor(tool,   levels = class_tools),
    class = factor(class,  levels = class_levels),
    metric = factor(metric, levels = c("F1", "MCC"))
  )
  
classif_colours <- tool_colours[class_tools]
p_metrics <- ggplot(metrics_long,
                   aes(x = tool, y = value, fill = tool)) +
  geom_col(width = 0.65, alpha = 0.85, position = "dodge") +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
  facet_grid(metric ~ class, scales = "free_y") +
  scale_fill_manual(values = classif_colours, name = "Tool") +
  scale_y_continuous(breaks = c(-0.5, 0, 0.25, 0.5, 0.75, 1)) +
  labs(
    x = NULL,
    y = "Score"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, size = 8),
    strip.background = element_rect(fill = "grey92"),
    strip.text = element_text(face = "bold", size = 9),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 8)
  )

# ggsave("classification_metrics.pdf", p_metrics,
#          width  = max(8, 2.5 * length(class_levels)),
#          height = 6,
#          device = "pdf")
