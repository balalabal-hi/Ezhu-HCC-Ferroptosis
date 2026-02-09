#!/usr/bin/env Rscript

# 07a_fig2_jimr_deg_ferroptosis.R
# JIMR-friendly Figure 2: DEG + heatmap + DEG∩Ferroptosis (no TCM/compound layer)
#
# Outputs:
# - plots/publication/JIMR_Figure2_deg_ferroptosis.pdf
# - plots/publication/JIMR_Figure2_deg_ferroptosis.png
#
# Usage:
#   Rscript scripts_final/07a_fig2_jimr_deg_ferroptosis.R

if (!interactive()) pdf(NULL)

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
  library(cowplot)
})

detect_base_dir <- function() {
  candidates <- c(".", "Ezhu_HCC_Project")
  for (d in candidates) {
    if (dir.exists(file.path(d, "results")) && dir.exists(file.path(d, "plots")) && dir.exists(file.path(d, "data/processed"))) return(d)
  }
  stop("Cannot locate project root.")
}

base_dir <- detect_base_dir()
res_dir <- file.path(base_dir, "results")
proc_dir <- file.path(base_dir, "data/processed")
ref_dir <- file.path(base_dir, "data/references")
pub_dir <- file.path(base_dir, "plots/publication")
dir.create(pub_dir, showWarnings = FALSE, recursive = TRUE)

theme_pub <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0.5),
      axis.title = element_text(size = base_size, face = "bold"),
      axis.text = element_text(size = base_size - 1, color = "#2B3A3F"),
      axis.line = element_line(color = "#9FB1B6", linewidth = 0.4),
      panel.grid.major = element_line(color = "#E6ECEE"),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )
}

deg_all <- read.csv(file.path(res_dir, "deg_GSE14520_all.csv"))
deg_all$Sig <- "NS"
deg_all$Sig[deg_all$adj.P.Val < 0.05 & deg_all$logFC > 1] <- "Up"
deg_all$Sig[deg_all$adj.P.Val < 0.05 & deg_all$logFC < -1] <- "Down"

p2a <- ggplot(deg_all, aes(x = logFC, y = -log10(adj.P.Val + 1e-300), color = Sig)) +
  geom_point(alpha = 0.3, size = 0.7) +
  scale_color_manual(values = c(Up = "#6FA8B5", Down = "#2F4F5F", NS = "#C8D7DB")) +
  geom_text_repel(
    data = head(deg_all %>% arrange(adj.P.Val), 10),
    aes(label = Gene),
    size = 3,
    fontface = "bold",
    color = "black"
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "#9FB1B6", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#9FB1B6", linewidth = 0.4) +
  theme_pub() +
  labs(title = "Volcano Plot (Tumor vs Non-tumor)", y = "-log10(FDR)")

expr_mat <- readRDS(file.path(proc_dir, "GSE14520_expr_symbol.rds"))
clinical_all <- readRDS(file.path(proc_dir, "GSE14520_clinical.rds"))
top_genes <- deg_all %>% arrange(adj.P.Val) %>% head(60) %>% pull(Gene)
top_genes <- intersect(top_genes, rownames(expr_mat)) %>% head(20)
expr_sub <- expr_mat[top_genes, , drop = FALSE]

group_info <- clinical_all %>%
  select(geo_accession, `Tissue:ch1`) %>%
  mutate(Group = ifelse(grepl("Tumor", `Tissue:ch1`, ignore.case = TRUE), "Tumor", "Non-Tumor")) %>%
  select(geo_accession, Group)

sample_order <- group_info$geo_accession[group_info$geo_accession %in% colnames(expr_sub)]
expr_long <- as.data.frame(expr_sub[, sample_order, drop = FALSE]) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "Sample", values_to = "Expression") %>%
  group_by(Gene) %>%
  mutate(Z = as.numeric(scale(Expression))) %>%
  ungroup() %>%
  left_join(group_info, by = c("Sample" = "geo_accession"))

sample_order <- group_info %>%
  filter(geo_accession %in% sample_order) %>%
  arrange(Group) %>%
  pull(geo_accession)

expr_long$Sample <- factor(expr_long$Sample, levels = sample_order)
expr_long$Gene <- factor(expr_long$Gene, levels = rev(top_genes))

p2b <- ggplot(expr_long, aes(x = Sample, y = Gene, fill = Z)) +
  geom_tile() +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B") +
  theme_pub() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(title = "Top DEGs Heatmap", x = NULL, y = NULL, fill = "Z-score")

ferro_context_path <- file.path(res_dir, "ferroptosis_genes_hcc_context.csv")
ferro_ref_path <- file.path(ref_dir, "ferroptosis_genes_expanded.csv")
ferro_path <- if (file.exists(ferro_context_path)) ferro_context_path else ferro_ref_path
if (!file.exists(ferro_path)) stop("Missing ferroptosis gene set.")
ferro_genes <- read.csv(ferro_path)$Gene

deg_sig <- deg_all %>% filter(adj.P.Val < 0.05 & abs(logFC) > 1) %>% pull(Gene)
deg_set <- unique(deg_sig)
ferro_set <- unique(ferro_genes)

deg_only <- length(setdiff(deg_set, ferro_set))
ferro_only <- length(setdiff(ferro_set, deg_set))
deg_ferro <- length(intersect(deg_set, ferro_set))

circle_points <- function(x0, y0, r, n = 200) {
  t <- seq(0, 2 * pi, length.out = n)
  data.frame(x = x0 + r * cos(t), y = y0 + r * sin(t))
}

deg_circle <- circle_points(-0.8, 0, 1.7)
ferro_circle <- circle_points(0.8, 0, 1.7)

p2c <- ggplot() +
  geom_polygon(data = deg_circle, aes(x = x, y = y), fill = "#3C5488FF", alpha = 0.20, color = "#3C5488FF") +
  geom_polygon(data = ferro_circle, aes(x = x, y = y), fill = "#E64B35FF", alpha = 0.20, color = "#E64B35FF") +
  annotate("text", x = -2.3, y = -1.4, label = paste0("DEG\n(n=", length(deg_set), ")"), size = 3.6, color = "#3C5488FF") +
  annotate("text", x = 2.3, y = -1.4, label = paste0("Ferroptosis\n(n=", length(ferro_set), ")"), size = 3.6, color = "#E64B35FF") +
  annotate("text", x = -1.8, y = 0.1, label = deg_only, size = 4) +
  annotate("text", x = 1.8, y = 0.1, label = ferro_only, size = 4) +
  annotate("text", x = 0, y = 0.1, label = deg_ferro, size = 5, fontface = "bold") +
  coord_equal() +
  theme_void() +
  labs(title = "HCC-context ferroptosis candidates (DEG ∩ Ferroptosis)")

candidate_df <- deg_all %>%
  filter(Gene %in% intersect(deg_set, ferro_set)) %>%
  group_by(Gene) %>%
  slice_min(order_by = adj.P.Val, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    NegLogFDR = -log10(adj.P.Val + 1e-300),
    Direction = ifelse(logFC >= 0, "Up", "Down")
  ) %>%
  arrange(desc(NegLogFDR)) %>%
  head(20) %>%
  arrange(logFC) %>%
  mutate(Gene = factor(Gene, levels = Gene))

p2d <- ggplot(candidate_df, aes(x = logFC, y = Gene)) +
  geom_segment(aes(x = 0, xend = logFC, y = Gene, yend = Gene),
               color = "#C8D7DB", linewidth = 0.6) +
  geom_point(aes(size = NegLogFDR, color = Direction), alpha = 0.9) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#9FB1B6", linewidth = 0.4) +
  scale_color_manual(values = c(Up = "#6FA8B5", Down = "#2F4F5F")) +
  scale_size_continuous(range = c(2.2, 5.2)) +
  theme_pub() +
  labs(
    title = "Top HCC-context ferroptosis candidates",
    x = "logFC (tumor vs non-tumor)",
    y = NULL,
    color = "Direction",
    size = "-log10(FDR)"
  )

fig2 <- plot_grid(p2a, p2b, p2c, p2d, labels = c("A", "B", "C", "D"), ncol = 2, rel_widths = c(1.1, 0.9))

ggsave(file.path(pub_dir, "JIMR_Figure2_deg_ferroptosis.pdf"), fig2, width = 12, height = 9, device = cairo_pdf)
ggsave(file.path(pub_dir, "JIMR_Figure2_deg_ferroptosis.png"), fig2, width = 12, height = 9, dpi = 300, bg = "white")
message("✅ Saved: plots/publication/JIMR_Figure2_deg_ferroptosis (PDF+PNG)")
