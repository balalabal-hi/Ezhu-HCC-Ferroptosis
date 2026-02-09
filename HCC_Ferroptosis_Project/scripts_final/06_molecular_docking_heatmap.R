#!/usr/bin/env Rscript
# =============================================================================
# Boltz-2 分子对接亲和力热图
#
# 输入: results/boltz2_affinity_results.csv
# 输出: plots/Figure11_molecular_docking_heatmap.pdf
# =============================================================================

# 加载包
library(ggplot2)
library(reshape2)

# 路径兼容：允许在仓库根目录或项目根目录运行
base_dir <- "."
if (!file.exists(file.path(base_dir, "results", "boltz2_affinity_results.csv"))) {
  if (file.exists(file.path("Ezhu_HCC_Project", "results", "boltz2_affinity_results.csv"))) {
    base_dir <- "Ezhu_HCC_Project"
  } else {
    stop("找不到 results/boltz2_affinity_results.csv；请在 Ezhu_HCC_Project/ 目录运行或提供完整路径。")
  }
}

# 创建输出目录
dir.create(file.path(base_dir, "plots"), showWarnings = FALSE, recursive = TRUE)

# 读取数据
boltz <- read.csv(file.path(base_dir, "results", "boltz2_affinity_results.csv"), stringsAsFactors = FALSE)
boltz <- subset(boltz, Status == "OK")
if (nrow(boltz) == 0) {
  stop("No Boltz-2 results found (Status == OK).")
}

# Prefer pIC50 (more interpretable; higher is better). Fallback to raw score.
score_col <- if ("Affinity_pIC50" %in% colnames(boltz)) "Affinity_pIC50" else "Affinity_Pred_Value"
score_higher_better <- score_col == "Affinity_pIC50"

# 转换为矩阵格式
docking_matrix <- dcast(boltz, Protein ~ Ligand, value.var = score_col)
rownames(docking_matrix) <- docking_matrix$Protein
docking_matrix$Protein <- NULL

# 转换为长格式用于ggplot
docking_long <- melt(as.matrix(docking_matrix))
colnames(docking_long) <- c("Protein", "Ligand", "Score")

# 设置蛋白和配体顺序（按平均亲和力排序）
protein_order <- aggregate(Score ~ Protein, docking_long, mean, na.rm = TRUE)
protein_order <- if (score_higher_better) {
  protein_order[order(-protein_order$Score), "Protein"]
} else {
  protein_order[order(protein_order$Score), "Protein"]
}
ligand_order <- aggregate(Score ~ Ligand, docking_long, mean, na.rm = TRUE)
ligand_order <- if (score_higher_better) {
  ligand_order[order(-ligand_order$Score), "Ligand"]
} else {
  ligand_order[order(ligand_order$Score), "Ligand"]
}

docking_long$Protein <- factor(docking_long$Protein, levels = protein_order)
docking_long$Ligand <- factor(docking_long$Ligand, levels = ligand_order)

# 生成热图
score_min <- min(docking_long$Score, na.rm = TRUE)
score_max <- max(docking_long$Score, na.rm = TRUE)
midpoint <- (score_min + score_max) / 2
legend_name <- if (score_higher_better) "pIC50 (higher is better)" else "Affinity score (lower is better)"

p <- ggplot(docking_long, aes(x = Ligand, y = Protein, fill = Score)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", Score)),
            color = "white", size = 5, fontface = "bold") +
  scale_fill_gradient2(
    low = "#4575B4",
    mid = "#FEE090",
    high = "#D73027",
    midpoint = midpoint,
    name = legend_name,
    limits = c(score_min, score_max)
  ) +
  labs(title = "Molecular Docking Results (Boltz-2)",
       subtitle = "Predicted affinity between candidate ligands and hub proteins",
       x = "Ligand (Active Compounds)",
       y = "Protein (Hub Genes)") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray40"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 13),
    legend.position = "right",
    legend.title = element_text(size = 11),
    panel.grid = element_blank()
  )

# 保存图片
ggsave(file.path(base_dir, "plots", "Figure11_molecular_docking_heatmap.pdf"), p,
       width = 10, height = 6, dpi = 300)

cat("✅ 热图已保存: ", file.path(base_dir, "plots", "Figure11_molecular_docking_heatmap.pdf"), "\n", sep = "")

# 输出统计摘要
cat("\n=== 分子对接结果摘要 ===\n")
cat("已筛选 Status == OK 的 Boltz-2 结果。\n\n")

# 各蛋白最佳配体
for (prot in unique(boltz$Protein)) {
  best <- boltz[boltz$Protein == prot, ]
  if (score_higher_better) {
    best <- best[which.max(best[[score_col]]), ]
    cat(sprintf("%s 最佳配体(pIC50最高): %s (%.2f)\n",
                prot, best$Ligand, best[[score_col]]))
  } else {
    best <- best[which.min(best[[score_col]]), ]
    cat(sprintf("%s 最佳配体(亲和力分数最低): %s (%.2f)\n",
                prot, best$Ligand, best[[score_col]]))
  }
}

if (score_higher_better) {
  cat("\n整体最佳结合(pIC50最高): ")
  best_all <- boltz[which.max(boltz[[score_col]]), ]
} else {
  cat("\n整体最佳结合(亲和力分数最低): ")
  best_all <- boltz[which.min(boltz[[score_col]]), ]
}
cat(sprintf("%s + %s (%.2f)\n",
            best_all$Protein, best_all$Ligand, best_all[[score_col]]))
