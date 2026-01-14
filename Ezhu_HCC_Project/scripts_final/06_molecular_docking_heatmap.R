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

# 转换为矩阵格式
docking_matrix <- dcast(boltz, Protein ~ Ligand, value.var = "Affinity_Pred_Value")
rownames(docking_matrix) <- docking_matrix$Protein
docking_matrix$Protein <- NULL

# 转换为长格式用于ggplot
docking_long <- melt(as.matrix(docking_matrix))
colnames(docking_long) <- c("Protein", "Ligand", "Affinity")

# 设置蛋白和配体顺序（按平均亲和力排序）
protein_order <- aggregate(Affinity ~ Protein, docking_long, mean, na.rm = TRUE)
protein_order <- protein_order[order(-protein_order$Affinity), "Protein"]
ligand_order <- aggregate(Affinity ~ Ligand, docking_long, mean, na.rm = TRUE)
ligand_order <- ligand_order[order(-ligand_order$Affinity), "Ligand"]

docking_long$Protein <- factor(docking_long$Protein, levels = protein_order)
docking_long$Ligand <- factor(docking_long$Ligand, levels = ligand_order)

# 生成热图
score_min <- min(docking_long$Affinity, na.rm = TRUE)
score_max <- max(docking_long$Affinity, na.rm = TRUE)
midpoint <- (score_min + score_max) / 2

p <- ggplot(docking_long, aes(x = Ligand, y = Protein, fill = Affinity)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", Affinity)),
            color = "white", size = 5, fontface = "bold") +
  scale_fill_gradient2(low = "#313695", mid = "#74ADD1", high = "#FEE090",
                       midpoint = midpoint,
                       name = "Affinity Score",
                       limits = c(score_min, score_max)) +
  labs(title = "Molecular Docking Results (Boltz-2)",
       subtitle = "Predicted affinity between Ezhu ligands and hub proteins",
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

cat("热图已保存: ", file.path(base_dir, "plots", "Figure11_molecular_docking_heatmap.pdf"), "\n", sep = "")

# 输出统计摘要
cat("\n=== 分子对接结果摘要 ===\n")
cat("已筛选 Status == OK 的 Boltz-2 结果。\n\n")

# 各蛋白最佳配体
for (prot in unique(boltz$Protein)) {
  best <- boltz[boltz$Protein == prot, ]
  best <- best[which.max(best$Affinity_Pred_Value), ]
  cat(sprintf("%s 最佳配体: %s (%.2f)\n",
              prot, best$Ligand, best$Affinity_Pred_Value))
}

cat("\n整体最佳结合: ")
best_all <- boltz[which.max(boltz$Affinity_Pred_Value), ]
cat(sprintf("%s + %s (%.2f)\n",
            best_all$Protein, best_all$Ligand, best_all$Affinity_Pred_Value))
