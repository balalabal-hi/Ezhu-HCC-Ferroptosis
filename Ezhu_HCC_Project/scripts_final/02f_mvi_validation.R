#!/usr/bin/env Rscript

# 02f_mvi_validation.R
# MVI supplemental validation using PMC8692135 (ijbsv18p0261s2.csv)

options(device = pdf)
graphics.off()
Sys.setenv(DISPLAY = "")

suppressPackageStartupMessages({
  library(tidyverse)
  library(pROC)
  library(cowplot)
})

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Package 'data.table' is required for MVI validation.")
}

proc_dir <- "data/processed"
raw_dir <- "data/raw"
res_dir <- "results"
plot_dir <- "plots/supplementary"

if (!dir.exists(res_dir)) dir.create(res_dir, recursive = TRUE)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

mvi_path <- file.path(raw_dir, "ijbsv18p0261s2.csv")
alt_path <- file.path(Sys.getenv("HOME"), "Downloads", "ijbsv18p0261s2.csv")

if (!file.exists(mvi_path)) {
  if (file.exists(alt_path)) {
    file.copy(alt_path, mvi_path)
    message("[MVI] Copied MVI file from Downloads to data/raw/.")
  } else {
    stop("Missing ijbsv18p0261s2.csv in data/raw/ and Downloads.")
  }
}

coef_path <- file.path(res_dir, "prognostic_model_coef_ezhu.csv")
if (!file.exists(coef_path)) {
  stop("Missing prognostic_model_coef_ezhu.csv. Run prognostic model first.")
}

coef_df <- read.csv(coef_path)
coef_vec <- setNames(coef_df$Coefficient, coef_df$Gene)

header <- data.table::fread(mvi_path, nrows = 0, showProgress = FALSE)
available_cols <- colnames(header)

needed_cols <- c("Sample_ID", "MVI", names(coef_vec))
select_cols <- intersect(needed_cols, available_cols)
missing_cols <- setdiff(needed_cols, select_cols)

if (length(missing_cols) > 0) {
  message("[MVI] Missing columns in MVI file: ", paste(missing_cols, collapse = ", "))
}

if (!all(c("Sample_ID", "MVI") %in% select_cols)) {
  stop("MVI file must include Sample_ID and MVI columns.")
}

mvi_dt <- data.table::fread(mvi_path, select = select_cols, showProgress = FALSE)
mvi_dt <- as.data.frame(mvi_dt)

# Keep MVI vs non_MVI only
mvi_dt$MVI <- as.character(mvi_dt$MVI)
keep <- mvi_dt$MVI %in% c("MVI", "non_MVI")
mvi_dt <- mvi_dt[keep, , drop = FALSE]

expr_cols <- intersect(names(coef_vec), colnames(mvi_dt))
if (length(expr_cols) < 3) {
  stop("Too few model genes found in MVI file.")
}

expr_mat <- as.matrix(mvi_dt[, expr_cols, drop = FALSE])
mode(expr_mat) <- "numeric"
expr_mat <- t(expr_mat)

zscore_by_gene <- function(mat) {
  z <- t(scale(t(mat)))
  z[is.na(z)] <- 0
  z
}

expr_scaled <- zscore_by_gene(expr_mat)
common_genes <- intersect(rownames(expr_scaled), names(coef_vec))

risk_score <- colSums(expr_scaled[common_genes, , drop = FALSE] * coef_vec[common_genes])

mvi_res <- data.frame(
  sample_id = mvi_dt$Sample_ID,
  mvi_group = mvi_dt$MVI,
  risk_score = as.numeric(risk_score)
)

write.csv(mvi_res, file.path(res_dir, "mvi_risk_score.csv"), row.names = FALSE)

wilcox_p <- wilcox.test(risk_score ~ mvi_group, data = mvi_res)$p.value
roc_obj <- pROC::roc(mvi_res$mvi_group, mvi_res$risk_score, levels = c("non_MVI", "MVI"))
auc_val <- as.numeric(pROC::auc(roc_obj))

stats <- data.frame(
  Metric = c("MVI_n", "non_MVI_n", "Wilcoxon_P", "AUC"),
  Value = c(sum(mvi_res$mvi_group == "MVI"), sum(mvi_res$mvi_group == "non_MVI"), wilcox_p, auc_val)
)
write.csv(stats, file.path(res_dir, "mvi_validation_stats.csv"), row.names = FALSE)

p <- ggplot(mvi_res, aes(x = mvi_group, y = risk_score, fill = mvi_group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.2, alpha = 0.6) +
  scale_fill_manual(values = c("non_MVI" = "#2E9FDF", "MVI" = "#E7B800")) +
  labs(title = "Risk Score by MVI Status", x = "MVI Group", y = "Risk Score") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

p_roc <- ggplot(data.frame(tpr = roc_obj$sensitivities, fpr = 1 - roc_obj$specificities),
                aes(x = fpr, y = tpr)) +
  geom_path(linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  coord_equal() +
  labs(title = sprintf("MVI ROC (AUC = %.3f)", auc_val), x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 12)

fig <- cowplot::plot_grid(p, p_roc, ncol = 2, labels = c("A", "B"))

ggsave(file.path(plot_dir, "FigureS_mvi_validation.pdf"), fig, width = 9, height = 4.5)

message("[MVI] 完成")
message("  - results/mvi_risk_score.csv")
message("  - results/mvi_validation_stats.csv")
message("  - plots/supplementary/FigureS_mvi_validation.pdf")
