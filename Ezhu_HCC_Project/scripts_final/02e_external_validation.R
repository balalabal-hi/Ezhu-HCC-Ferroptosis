#!/usr/bin/env Rscript

# 02e_external_validation.R
# 外部验证：TCGA-LIHC + GSE76427 + GSE10143
# 输出：risk score、外部验证统计表、可选图版

# 禁用交互式图形设备
options(device = pdf)
graphics.off()
Sys.setenv(DISPLAY = "")

suppressPackageStartupMessages({
  library(tidyverse)
  library(GEOquery)
  library(survival)
  library(survminer)
  library(timeROC)
  library(cowplot)
})

if (!dir.exists("data/processed")) {
  stop("请在项目根目录运行此脚本")
}

proc_dir <- "data/processed"
res_dir <- "results"
plot_dir <- "plots/supplementary"

dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

coef_df <- read.csv(file.path(res_dir, "prognostic_model_coef_ezhu.csv"))
coef_vec <- setNames(coef_df$Coefficient, coef_df$Gene)
model_genes <- names(coef_vec)

calc_risk_score <- function(expr_mat, coef_vec) {
  common <- intersect(rownames(expr_mat), names(coef_vec))
  if (length(common) < 1) {
    stop("可用基因过少，无法计算风险评分")
  }
  score <- colSums(expr_mat[common, , drop = FALSE] * coef_vec[common])
  score
}

map_expression_to_symbols <- function(es) {
  expr_raw <- exprs(es)
  feat <- fData(es)
  gene_col <- grep("gene symbol|gene_symbol|symbol", colnames(feat), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(gene_col)) {
    stop("找不到基因符号列")
  }
  symbols <- feat[[gene_col]]
  symbols <- trimws(symbols)
  symbols <- sub(" ///.*$", "", symbols)
  symbols <- sub(" //.*$", "", symbols)
  symbols <- sub(";.*$", "", symbols)

  keep <- !is.na(symbols) & symbols != ""
  expr_raw <- expr_raw[keep, , drop = FALSE]
  symbols <- symbols[keep]

  expr_sum <- rowsum(expr_raw, symbols)
  expr_cnt <- rowsum(matrix(1, nrow(expr_raw), ncol(expr_raw)), symbols)
  expr_sum / expr_cnt
}

map_expression_with_gpl <- function(expr_raw, gpl_id) {
  gpl <- getGEO(gpl_id)
  anno <- Table(gpl)
  id_col <- grep("^ID$|^ID_REF$", colnames(anno), ignore.case = TRUE, value = TRUE)[1]
  gene_col <- grep("gene symbol|gene_symbol|symbol", colnames(anno), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(id_col) || is.na(gene_col)) {
    stop("GPL 注释表缺少探针或基因符号列")
  }

  anno <- anno[, c(id_col, gene_col)]
  colnames(anno) <- c("Probe", "Symbol")
  anno$Symbol <- trimws(anno$Symbol)
  anno$Symbol <- sub(" ///.*$", "", anno$Symbol)
  anno$Symbol <- sub(" //.*$", "", anno$Symbol)
  anno$Symbol <- sub(";.*$", "", anno$Symbol)

  probe <- rownames(expr_raw)
  matched <- anno[match(probe, anno$Probe), ]
  symbols <- matched$Symbol

  keep <- !is.na(symbols) & symbols != ""
  expr_raw <- expr_raw[keep, , drop = FALSE]
  symbols <- symbols[keep]

  expr_sum <- rowsum(expr_raw, symbols)
  expr_cnt <- rowsum(matrix(1, nrow(expr_raw), ncol(expr_raw)), symbols)
  expr_sum / expr_cnt
}

zscore_by_gene <- function(expr_mat) {
  z <- t(scale(t(expr_mat)))
  z[is.na(z)] <- 0
  z
}

maybe_log2p1 <- function(expr_mat) {
  q99 <- suppressWarnings(as.numeric(quantile(as.numeric(expr_mat), 0.99, na.rm = TRUE)))
  if (!is.na(q99) && q99 > 100) {
    # guard against negative values from normalized matrices
    return(log2(pmax(expr_mat, 0) + 1))
  }
  expr_mat
}

train_frozen_scaler <- function(model_genes) {
  train_expr_path <- file.path(proc_dir, "GSE14520_expr_symbol.rds")
  train_ids_path <- file.path(res_dir, "risk_score_data_ezhu.csv")
  if (!file.exists(train_expr_path) || !file.exists(train_ids_path)) {
    stop("缺少训练集表达或风险评分文件，无法构建冻结尺度参数")
  }
  train_expr <- readRDS(train_expr_path)
  train_expr <- maybe_log2p1(train_expr)
  train_ids <- read.csv(train_ids_path)$sample
  train_ids <- intersect(colnames(train_expr), train_ids)
  if (length(train_ids) < 20) stop("训练集样本过少，无法估计冻结尺度参数")
  genes <- intersect(rownames(train_expr), model_genes)
  if (length(genes) < 2) stop("训练集缺少模型基因，无法构建冻结尺度参数")
  mu <- rowMeans(train_expr[genes, train_ids, drop = FALSE], na.rm = TRUE)
  sdv <- apply(train_expr[genes, train_ids, drop = FALSE], 1, sd, na.rm = TRUE)
  sdv[is.na(sdv) | sdv == 0] <- 1
  list(mu = mu, sd = sdv, genes = genes)
}

zscore_by_train <- function(expr_mat, scaler) {
  expr_mat <- maybe_log2p1(expr_mat)
  common <- intersect(rownames(expr_mat), scaler$genes)
  if (length(common) < 1) stop("外部队列缺少模型基因，无法冻结尺度标准化")
  z <- expr_mat[common, , drop = FALSE]
  z <- (z - scaler$mu[common]) / scaler$sd[common]
  z[is.na(z)] <- 0
  z
}

extract_characteristic <- function(pdat, pattern) {
  char_cols <- grep("characteristics", colnames(pdat), ignore.case = TRUE, value = TRUE)
  if (length(char_cols) == 0) return(rep(NA_character_, nrow(pdat)))
  pdat[char_cols] <- lapply(pdat[char_cols], as.character)
  apply(pdat[, char_cols, drop = FALSE], 1, function(x) {
    hit <- x[grepl(pattern, x, ignore.case = TRUE)]
    if (length(hit) == 0) return(NA_character_)
    hit[1]
  })
}

extract_last_number <- function(x) {
  if (is.na(x)) return(NA_real_)
  m <- sub(".*:\\s*([0-9]+)\\s*$", "\\1", x)
  if (identical(m, x)) return(NA_real_)
  suppressWarnings(as.numeric(m))
}

extract_yes_no <- function(x) {
  if (is.na(x)) return(NA_character_)
  m <- sub(".*:\\s*([Yy]es|[Nn]o)\\b.*", "\\1", x)
  if (identical(m, x)) return(NA_character_)
  ifelse(tolower(m) == "yes", "Yes", "No")
}

collapse_characteristics <- function(pdat) {
  char_cols <- grep("^characteristics", colnames(pdat), ignore.case = TRUE, value = TRUE)
  if (length(char_cols) == 0) {
    return(rep("", nrow(pdat)))
  }
  apply(pdat[, char_cols, drop = FALSE], 1, function(x) {
    paste(na.omit(as.character(x)), collapse = ";")
  })
}

filter_hcc_samples <- function(pdat) {
  text <- tolower(collapse_characteristics(pdat))
  is_hcc <- grepl("hcc|hepatocellular|hepatocarcinoma|tumou?r|cancer", text)
  exclude <- grepl("non[- ]?tumor|nontumor|normal|hepatitis|cirrhosis|fibrosis", text)
  keep <- is_hcc & !exclude
  keep
}

summarize_validation <- function(df, label) {
  df <- df %>%
    filter(!is.na(time_months), !is.na(status), !is.na(risk_score))

  df$risk_group <- ifelse(df$risk_score >= median(df$risk_score), "High", "Low")

  surv_obj <- Surv(df$time_months, df$status)
  cox_fit <- coxph(surv_obj ~ risk_score, data = df)
  cindex <- summary(cox_fit)$concordance[1]

  logrank <- survdiff(surv_obj ~ risk_group, data = df)
  logrank_p <- 1 - pchisq(logrank$chisq, length(logrank$n) - 1)

  max_time <- max(df$time_months, na.rm = TRUE)
  roc_times <- c(12, 36, 60)
  roc_times <- roc_times[roc_times < max_time]

  auc_1y <- NA
  auc_3y <- NA
  auc_5y <- NA
  roc_obj <- NULL

  if (length(roc_times) > 0) {
    roc_obj <- timeROC(
      T = df$time_months,
      delta = df$status,
      marker = df$risk_score,
      cause = 1,
      weighting = "marginal",
      times = roc_times,
      iid = FALSE
    )

    if (12 %in% roc_times) auc_1y <- roc_obj$AUC[roc_obj$times == 12]
    if (36 %in% roc_times) auc_3y <- roc_obj$AUC[roc_obj$times == 36]
    if (60 %in% roc_times) auc_5y <- roc_obj$AUC[roc_obj$times == 60]
  }

  stats <- data.frame(
    Dataset = label,
    Sample_N = nrow(df),
    Events = sum(df$status == 1, na.rm = TRUE),
    AUC_1y = auc_1y,
    AUC_3y = auc_3y,
    AUC_5y = auc_5y,
    Logrank_P = logrank_p,
    C_index = cindex
  )

  list(stats = stats, data = df, roc = roc_obj)
}

summarize_validation_dualcut <- function(df, label, train_cutoff) {
  df <- df %>%
    filter(!is.na(time_months), !is.na(status), !is.na(risk_score))

  df$risk_group_median <- ifelse(df$risk_score >= median(df$risk_score), "High", "Low")
  df$risk_group_train <- ifelse(df$risk_score >= train_cutoff, "High", "Low")

  surv_obj <- Surv(df$time_months, df$status)
  cox_fit <- coxph(surv_obj ~ risk_score, data = df)
  cindex <- summary(cox_fit)$concordance[1]

  logrank_median <- survdiff(surv_obj ~ risk_group_median, data = df)
  logrank_p_median <- 1 - pchisq(logrank_median$chisq, length(logrank_median$n) - 1)

  logrank_p_train <- NA_real_
  if (length(unique(df$risk_group_train)) > 1) {
    logrank_train <- survdiff(surv_obj ~ risk_group_train, data = df)
    logrank_p_train <- 1 - pchisq(logrank_train$chisq, length(logrank_train$n) - 1)
  }

  max_time <- max(df$time_months, na.rm = TRUE)
  roc_times <- c(12, 36, 60)
  roc_times <- roc_times[roc_times < max_time]

  auc_1y <- NA
  auc_3y <- NA
  auc_5y <- NA
  roc_obj <- NULL

  if (length(roc_times) > 0) {
    roc_obj <- timeROC(
      T = df$time_months,
      delta = df$status,
      marker = df$risk_score,
      cause = 1,
      weighting = "marginal",
      times = roc_times,
      iid = FALSE
    )

    if (12 %in% roc_times) auc_1y <- roc_obj$AUC[roc_obj$times == 12]
    if (36 %in% roc_times) auc_3y <- roc_obj$AUC[roc_obj$times == 36]
    if (60 %in% roc_times) auc_5y <- roc_obj$AUC[roc_obj$times == 60]
  }

  stats <- data.frame(
    Dataset = label,
    Sample_N = nrow(df),
    Events = sum(df$status == 1, na.rm = TRUE),
    AUC_1y = auc_1y,
    AUC_3y = auc_3y,
    AUC_5y = auc_5y,
    Logrank_P_median = logrank_p_median,
    Logrank_P_traincut = logrank_p_train,
    C_index = cindex
  )

  list(stats = stats, data = df, roc = roc_obj)
}

plot_km <- function(df, title) {
  km_fit <- survfit(Surv(time_months, status) ~ risk_group, data = df)
  ggsurvplot(
    km_fit,
    data = df,
    pval = TRUE,
    risk.table = FALSE,
    palette = c("#2E9FDF", "#E7B800"),
    xlab = "Time (months)",
    ylab = "Overall Survival",
    legend.title = "",
    legend.labs = c("Low Risk", "High Risk"),
    ggtheme = theme_bw(base_size = 10)
  )$plot +
    labs(title = title)
}

plot_roc <- function(df, roc_obj, title) {
  if (is.null(roc_obj)) {
    return(ggplot() + theme_void() + labs(title = paste0(title, " (ROC NA)")))
  }
  roc_df <- data.frame()
  for (t in roc_obj$times) {
    idx <- which(roc_obj$times == t)
    tmp <- data.frame(
      FPR = roc_obj$FP[, idx],
      TPR = roc_obj$TP[, idx],
      Time = paste0(t / 12, "-year (AUC=", round(roc_obj$AUC[idx], 3), ")")
    )
    roc_df <- rbind(roc_df, tmp)
  }

  ggplot(roc_df, aes(x = FPR, y = TPR, color = Time)) +
    geom_path(linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(title = title, x = "1 - Specificity", y = "Sensitivity", color = "") +
    theme_bw(base_size = 10) +
    theme(legend.position = c(0.65, 0.25))
}

message("[外部验证] 开始 TCGA-LIHC 验证...")
expr_tcga <- readRDS(file.path(proc_dir, "TCGA_LIHC_expr.rds"))
clin_tcga <- readRDS(file.path(proc_dir, "TCGA_LIHC_clinical.rds"))

# Frozen-scaling sensitivity analysis (train-mean/train-sd; log2(x+1) if needed)
train_cutoff <- median(read.csv(file.path(res_dir, "risk_score_data_ezhu.csv"))$risk_score, na.rm = TRUE)
scaler <- train_frozen_scaler(model_genes)

patient_id <- substr(colnames(expr_tcga), 1, 12)
expr_tcga_t <- t(expr_tcga)
expr_tcga_sum <- rowsum(expr_tcga_t, patient_id)
expr_tcga_avg <- expr_tcga_sum / as.vector(table(patient_id))
expr_tcga_avg <- t(expr_tcga_avg)

risk_tcga <- calc_risk_score(zscore_by_gene(expr_tcga_avg), coef_vec)
tcga_df <- data.frame(
  patient_id = names(risk_tcga),
  risk_score = as.numeric(risk_tcga)
)

clin_tcga$patient_id <- clin_tcga$bcr_patient_barcode
clin_tcga$time_days <- ifelse(
  !is.na(clin_tcga$days_to_death),
  clin_tcga$days_to_death,
  clin_tcga$days_to_last_follow_up
)
clin_tcga$status <- ifelse(tolower(clin_tcga$vital_status) == "dead", 1, 0)
clin_tcga$time_months <- as.numeric(clin_tcga$time_days) / 30.4

val_tcga <- merge(tcga_df, clin_tcga, by = "patient_id")
val_tcga <- val_tcga %>% filter(!is.na(time_months), !is.na(status))

write.csv(val_tcga, file.path(res_dir, "TCGA_LIHC_risk_score.csv"), row.names = FALSE)

res_tcga <- summarize_validation(val_tcga, "TCGA-LIHC")

risk_tcga_frozen <- calc_risk_score(zscore_by_train(expr_tcga_avg, scaler), coef_vec)
tcga_df_frozen <- data.frame(
  patient_id = names(risk_tcga_frozen),
  risk_score = as.numeric(risk_tcga_frozen)
)
val_tcga_frozen <- merge(tcga_df_frozen, clin_tcga, by = "patient_id") %>%
  filter(!is.na(time_months), !is.na(status))
res_tcga_frozen <- summarize_validation_dualcut(val_tcga_frozen, "TCGA-LIHC", train_cutoff)

message("[外部验证] 开始 GSE76427 验证...")
expr_76427_path <- file.path(proc_dir, "GSE76427_expr_symbol.rds")
clin_76427_path <- file.path(proc_dir, "GSE76427_clinical.rds")

if (!file.exists(expr_76427_path) || !file.exists(clin_76427_path)) {
  g <- getGEO("GSE76427", GSEMatrix = TRUE, AnnotGPL = TRUE)
  es <- if (is.list(g)) g[[1]] else g
  pheno <- pData(es)
  expr_avg <- map_expression_to_symbols(es)

  pheno$time_months <- as.numeric(pheno[["duryears_os:ch1"]]) * 12
  pheno$status <- as.numeric(pheno[["event_os:ch1"]])

  saveRDS(expr_avg, expr_76427_path)
  saveRDS(pheno, clin_76427_path)
}

expr_76427 <- readRDS(expr_76427_path)
clin_76427 <- readRDS(clin_76427_path)

common_samples <- intersect(colnames(expr_76427), rownames(clin_76427))
expr_76427 <- expr_76427[, common_samples, drop = FALSE]
clin_76427 <- clin_76427[common_samples, , drop = FALSE]

risk_76427 <- calc_risk_score(zscore_by_gene(expr_76427), coef_vec)
val_76427 <- data.frame(
  sample_id = names(risk_76427),
  risk_score = as.numeric(risk_76427),
  time_months = as.numeric(clin_76427$time_months),
  status = as.numeric(clin_76427$status)
)
val_76427 <- val_76427 %>% filter(!is.na(time_months), !is.na(status))

write.csv(val_76427, file.path(res_dir, "GSE76427_risk_score.csv"), row.names = FALSE)

res_76427 <- summarize_validation(val_76427, "GSE76427")

risk_76427_frozen <- calc_risk_score(zscore_by_train(expr_76427, scaler), coef_vec)
val_76427_frozen <- data.frame(
  sample_id = names(risk_76427_frozen),
  risk_score = as.numeric(risk_76427_frozen),
  time_months = as.numeric(clin_76427$time_months[match(names(risk_76427_frozen), rownames(clin_76427))]),
  status = as.numeric(clin_76427$status[match(names(risk_76427_frozen), rownames(clin_76427))])
)
val_76427_frozen <- val_76427_frozen %>% filter(!is.na(time_months), !is.na(status))
res_76427_frozen <- summarize_validation_dualcut(val_76427_frozen, "GSE76427", train_cutoff)

message("[外部验证] 开始 GSE10143 验证...")
expr_10143_path <- file.path(proc_dir, "GSE10143_expr_symbol.rds")
clin_10143_path <- file.path(proc_dir, "GSE10143_clinical.rds")
expr_10143_hcc_path <- file.path(proc_dir, "GSE10143_expr_symbol_hcc.rds")
clin_10143_hcc_path <- file.path(proc_dir, "GSE10143_clinical_hcc.rds")

expr_10143 <- readRDS(expr_10143_path)
clin_10143 <- readRDS(clin_10143_path)

keep_hcc <- filter_hcc_samples(clin_10143)
message("[外部验证] GSE10143 HCC 样本数: ", sum(keep_hcc), "/", nrow(clin_10143))
clin_10143 <- clin_10143[keep_hcc, , drop = FALSE]

common_samples <- intersect(colnames(expr_10143), rownames(clin_10143))
expr_10143 <- expr_10143[, common_samples, drop = FALSE]
clin_10143 <- clin_10143[common_samples, , drop = FALSE]

saveRDS(expr_10143, expr_10143_hcc_path)
saveRDS(clin_10143, clin_10143_hcc_path)

risk_10143 <- calc_risk_score(zscore_by_gene(expr_10143), coef_vec)
val_10143 <- data.frame(
  sample_id = names(risk_10143),
  risk_score = as.numeric(risk_10143),
  time_months = as.numeric(clin_10143$time_months),
  status = as.numeric(clin_10143$status)
)
val_10143 <- val_10143 %>% filter(!is.na(time_months), !is.na(status))

write.csv(val_10143, file.path(res_dir, "GSE10143_HCC_risk_score.csv"), row.names = FALSE)

res_10143 <- summarize_validation(val_10143, "GSE10143-HCC")

risk_10143_frozen <- calc_risk_score(zscore_by_train(expr_10143, scaler), coef_vec)
val_10143_frozen <- data.frame(
  sample_id = names(risk_10143_frozen),
  risk_score = as.numeric(risk_10143_frozen),
  time_months = as.numeric(clin_10143$time_months[match(names(risk_10143_frozen), rownames(clin_10143))]),
  status = as.numeric(clin_10143$status[match(names(risk_10143_frozen), rownames(clin_10143))])
)
val_10143_frozen <- val_10143_frozen %>% filter(!is.na(time_months), !is.na(status))
res_10143_frozen <- summarize_validation_dualcut(val_10143_frozen, "GSE10143-HCC", train_cutoff)

stats <- bind_rows(res_tcga$stats, res_76427$stats, res_10143$stats)
write.csv(stats, file.path(res_dir, "external_validation_stats.csv"), row.names = FALSE)

stats_frozen <- bind_rows(
  transform(res_tcga_frozen$stats, Scaling = "Train-frozen"),
  transform(res_76427_frozen$stats, Scaling = "Train-frozen"),
  transform(res_10143_frozen$stats, Scaling = "Train-frozen")
)
write.csv(stats_frozen, file.path(res_dir, "external_validation_stats_frozen_scaling.csv"), row.names = FALSE)
write.csv(stats_frozen, file.path(res_dir, "Supplementary_Table_ExternalValidation_FrozenScaling.csv"), row.names = FALSE)

p_tcga_km <- plot_km(res_tcga$data, "TCGA-LIHC")
p_tcga_roc <- plot_roc(res_tcga$data, res_tcga$roc, "TCGA-LIHC ROC")

p_76427_km <- plot_km(res_76427$data, "GSE76427")
p_76427_roc <- plot_roc(res_76427$data, res_76427$roc, "GSE76427 ROC")

p_10143_km <- plot_km(res_10143$data, "GSE10143-HCC")
p_10143_roc <- plot_roc(res_10143$data, res_10143$roc, "GSE10143-HCC ROC")

fig <- plot_grid(
  p_tcga_km, p_tcga_roc,
  p_76427_km, p_76427_roc,
  p_10143_km, p_10143_roc,
  ncol = 2,
  labels = c("A", "B", "C", "D", "E", "F")
)

ggsave(file.path(plot_dir, "FigureS_external_validation.pdf"), fig, width = 12, height = 13)

message("[外部验证] 完成")
message("  - results/TCGA_LIHC_risk_score.csv")
message("  - results/GSE76427_risk_score.csv")
message("  - results/GSE10143_HCC_risk_score.csv")
message("  - results/external_validation_stats.csv")
message("  - plots/supplementary/FigureS_external_validation.pdf")
