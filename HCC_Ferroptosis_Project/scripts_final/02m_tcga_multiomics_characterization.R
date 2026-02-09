#!/usr/bin/env Rscript

# 02m_tcga_multiomics_characterization.R
# TCGA-LIHC 多组学表征（突变 + CNV；探索性）
# 输出：results/tcga_multiomics_*.csv + plots/supplementary/FigureS_tcga_multiomics.(png/pdf)

options(device = pdf)
graphics.off()
Sys.setenv(DISPLAY = "")

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(UCSCXenaTools)
  library(cowplot)
})

if (!dir.exists("data/processed")) {
  stop("请在项目根目录运行此脚本")
}

raw_xena_dir <- file.path("data", "raw", "xena")
res_dir <- "results"
plot_dir <- file.path("plots", "supplementary")

dir.create(raw_xena_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

ensure_xena_dataset <- function(dataset, hostName = "tcgaHub", destdir = raw_xena_dir) {
  destfile <- file.path(destdir, paste0(dataset, ".gz"))
  if (file.exists(destfile)) return(destfile)
  hub <- XenaHub(hostName = hostName, datasets = dataset)
  q <- XenaQuery(hub)
  XenaDownload(q, destdir = destdir, force = FALSE)
  if (!file.exists(destfile)) stop("Xena 下载失败：", dataset)
  destfile
}

message("[TCGA多组学] 读取 TCGA 风险评分...")
tcga <- read.csv(file.path(res_dir, "TCGA_LIHC_risk_score.csv"))
tcga <- tcga %>%
  select(patient_id, risk_score, time_months, status) %>%
  mutate(
    risk_group = ifelse(risk_score >= median(risk_score, na.rm = TRUE), "High", "Low")
  )

message("[TCGA多组学] 下载/读取 MC3 突变数据...")
mc3_path <- ensure_xena_dataset("mc3/LIHC_mc3.txt")
mc3_gene_level_path <- ensure_xena_dataset("mc3_gene_level/LIHC_mc3_gene_level.txt")

mc3 <- fread(mc3_path)
if (!all(c("sample", "gene", "effect") %in% colnames(mc3))) {
  stop("MC3 文件缺少必要列：sample/gene/effect")
}
mc3$patient_id <- substr(mc3$sample, 1, 12)

silent_effects <- c(
  "Silent", "Intron", "3'UTR", "5'UTR", "IGR", "RNA", "lincRNA", "Targeted_Region",
  "5'Flank", "3'Flank"
)
mc3$non_silent <- !(mc3$effect %in% silent_effects) & !is.na(mc3$effect) & mc3$effect != ""

mut_burden <- mc3[, .(
  mut_total = .N,
  mut_nonsilent = sum(non_silent, na.rm = TRUE)
), by = patient_id]

message("[TCGA多组学] 读取基因层面突变矩阵...")
mc3_gl <- fread(mc3_gene_level_path, data.table = FALSE)
gene_col <- colnames(mc3_gl)[1]
colnames(mc3_gl)[1] <- "gene"

key_genes <- c("TP53", "CTNNB1", "AXIN1", "ALB", "ARID1A", "ARID2", "KEAP1", "NFE2L2")
mc3_gl_key <- mc3_gl %>% filter(gene %in% key_genes)
if (nrow(mc3_gl_key) == 0) stop("未在 MC3 gene-level 文件中找到关键基因行")

mut_key_long <- mc3_gl_key %>%
  pivot_longer(cols = -gene, names_to = "sample", values_to = "mutated") %>%
  mutate(
    patient_id = substr(sample, 1, 12),
    mutated = as.integer(as.numeric(mutated) > 0)
  ) %>%
  group_by(patient_id, gene) %>%
  summarise(mutated = max(mutated, na.rm = TRUE), .groups = "drop")

message("[TCGA多组学] 下载/读取 GISTIC2 CNV 数据...")
gistic_path <- ensure_xena_dataset("TCGA.LIHC.sampleMap/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")
cnv_dt <- fread(gistic_path, data.table = FALSE, check.names = FALSE)
if (ncol(cnv_dt) < 10) stop("GISTIC2 CNV 文件列数异常，可能下载/解析失败")

gene_sym <- cnv_dt[[1]]
cnv_mat <- as.matrix(cnv_dt[, -1, drop = FALSE])
storage.mode(cnv_mat) <- "numeric"
colnames(cnv_mat) <- colnames(cnv_dt)[-1]
rownames(cnv_mat) <- make.unique(as.character(gene_sym))

cnv_burden_sample <- colMeans(abs(cnv_mat) >= 1, na.rm = TRUE)
cnv_burden <- data.frame(
  sample = names(cnv_burden_sample),
  cnv_burden = as.numeric(cnv_burden_sample)
) %>%
  mutate(patient_id = substr(sample, 1, 12)) %>%
  group_by(patient_id) %>%
  summarise(cnv_burden = mean(cnv_burden, na.rm = TRUE), .groups = "drop")

sig_genes <- read.csv(file.path(res_dir, "prognostic_model_coef_ezhu.csv"))$Gene
sig_genes <- intersect(sig_genes, rownames(cnv_mat))
cnv_sig <- NULL
if (length(sig_genes) >= 2) {
  cnv_sig_score <- colMeans(abs(cnv_mat[sig_genes, , drop = FALSE]), na.rm = TRUE)
  cnv_sig <- data.frame(
    sample = names(cnv_sig_score),
    cnv_sig_absmean = as.numeric(cnv_sig_score)
  ) %>%
    mutate(patient_id = substr(sample, 1, 12)) %>%
    group_by(patient_id) %>%
    summarise(cnv_sig_absmean = mean(cnv_sig_absmean, na.rm = TRUE), .groups = "drop")
}

message("[TCGA多组学] 合并多组学指标到 TCGA 风险评分...")
tcga_om <- tcga %>%
  left_join(mut_burden, by = "patient_id") %>%
  left_join(cnv_burden, by = "patient_id")
if (!is.null(cnv_sig)) tcga_om <- tcga_om %>% left_join(cnv_sig, by = "patient_id")

write.csv(tcga_om, file.path(res_dir, "tcga_multiomics_burden.csv"), row.names = FALSE)

message("[TCGA多组学] 关键基因突变与风险分组关联...")
mut_assoc <- mut_key_long %>%
  left_join(tcga %>% select(patient_id, risk_group), by = "patient_id") %>%
  filter(!is.na(risk_group)) %>%
  group_by(gene) %>%
  group_modify(~{
    x <- .x
    tab <- table(x$risk_group, x$mutated)
    # ensure 2x2
    if (!all(c("Low", "High") %in% rownames(tab))) {
      miss <- setdiff(c("Low", "High"), rownames(tab))
      for (m in miss) tab <- rbind(tab, setNames(c(0, 0), colnames(tab)))
      rownames(tab)[nrow(tab)] <- miss[1]
    }
    if (!all(c("0", "1") %in% colnames(tab))) {
      for (m in setdiff(c("0", "1"), colnames(tab))) tab <- cbind(tab, setNames(rep(0, nrow(tab)), m))
    }
    tab <- tab[c("Low", "High"), c("0", "1"), drop = FALSE]
    ft <- fisher.test(tab)
    data.frame(
      n_low = sum(tab["Low", ]),
      n_high = sum(tab["High", ]),
      mut_low = tab["Low", "1"],
      mut_high = tab["High", "1"],
      fisher_p = ft$p.value,
      odds_ratio = unname(ft$estimate)
    )
  }) %>%
  ungroup() %>%
  arrange(fisher_p)

write.csv(mut_assoc, file.path(res_dir, "tcga_multiomics_keygene_mutation_assoc.csv"), row.names = FALSE)

message("[TCGA多组学] 生成补充图...")
plot_df <- tcga_om %>% filter(!is.na(mut_nonsilent), !is.na(cnv_burden))

tp53 <- mut_key_long %>% filter(gene == "TP53") %>% select(patient_id, tp53_mut = mutated)
ctnnb1 <- mut_key_long %>% filter(gene == "CTNNB1") %>% select(patient_id, ctnnb1_mut = mutated)
plot_df <- plot_df %>%
  left_join(tp53, by = "patient_id") %>%
  left_join(ctnnb1, by = "patient_id") %>%
  mutate(
    tp53_mut = factor(ifelse(tp53_mut == 1, "Mut", "WT"), levels = c("WT", "Mut")),
    ctnnb1_mut = factor(ifelse(ctnnb1_mut == 1, "Mut", "WT"), levels = c("WT", "Mut"))
  )

p1 <- ggplot(plot_df, aes(x = tp53_mut, y = risk_score, fill = tp53_mut)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.35, size = 1) +
  scale_fill_manual(values = c("WT" = "#4C78A8", "Mut" = "#F58518")) +
  theme_bw(base_size = 10) +
  theme(legend.position = "none") +
  labs(title = "TP53 mutation vs risk score", x = "", y = "Risk score")

p2 <- ggplot(plot_df, aes(x = ctnnb1_mut, y = risk_score, fill = ctnnb1_mut)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.35, size = 1) +
  scale_fill_manual(values = c("WT" = "#4C78A8", "Mut" = "#F58518")) +
  theme_bw(base_size = 10) +
  theme(legend.position = "none") +
  labs(title = "CTNNB1 mutation vs risk score", x = "", y = "Risk score")

sp1 <- suppressWarnings(cor.test(plot_df$risk_score, plot_df$mut_nonsilent, method = "spearman"))
p3 <- ggplot(plot_df, aes(x = mut_nonsilent, y = risk_score)) +
  geom_point(alpha = 0.55, size = 1.2, color = "#4C78A8") +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
  theme_bw(base_size = 10) +
  labs(
    title = sprintf("Mutation burden vs risk (Spearman r=%.2f, P=%.2g)", sp1$estimate, sp1$p.value),
    x = "Non-silent mutation count (MC3)",
    y = "Risk score"
  )

sp2 <- suppressWarnings(cor.test(plot_df$risk_score, plot_df$cnv_burden, method = "spearman"))
p4 <- ggplot(plot_df, aes(x = cnv_burden, y = risk_score)) +
  geom_point(alpha = 0.55, size = 1.2, color = "#54A24B") +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
  theme_bw(base_size = 10) +
  labs(
    title = sprintf("CNV burden vs risk (Spearman r=%.2f, P=%.2g)", sp2$estimate, sp2$p.value),
    x = "CNV burden (fraction |GISTIC|≥1)",
    y = "Risk score"
  )

fig <- plot_grid(p1, p2, p3, p4, ncol = 2, labels = c("A", "B", "C", "D"))
ggsave(file.path(plot_dir, "FigureS_tcga_multiomics.pdf"), fig, width = 12, height = 8)
ggsave(file.path(plot_dir, "FigureS_tcga_multiomics.png"), fig, width = 12, height = 8, dpi = 300, bg = "white")

message("[TCGA多组学] 完成")
message("  ✅ results/tcga_multiomics_burden.csv")
message("  ✅ results/tcga_multiomics_keygene_mutation_assoc.csv")
message("  ✅ plots/supplementary/FigureS_tcga_multiomics.pdf")
message("  ✅ plots/supplementary/FigureS_tcga_multiomics.png")
