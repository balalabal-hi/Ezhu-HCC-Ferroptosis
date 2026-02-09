#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(cowplot)
})

detect_base_dir <- function() {
  candidates <- c(".", "Ezhu_HCC_Project")
  for (d in candidates) {
    if (file.exists(file.path(d, "results", "boltz2_affinity_results_modelgenes.csv"))) {
      return(d)
    }
  }
  stop("Cannot locate project base directory with model-gene Boltz-2 results.")
}

base_dir <- detect_base_dir()
res_file <- file.path(base_dir, "results", "boltz2_affinity_results_modelgenes.csv")
out_dir <- file.path(base_dir, "plots", "publication")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

boltz <- read.csv(res_file, stringsAsFactors = FALSE)
boltz$Affinity_pIC50 <- suppressWarnings(as.numeric(boltz$Affinity_pIC50))
boltz$Confidence <- suppressWarnings(as.numeric(boltz$Confidence))

score_col <- if ("Affinity_pIC50" %in% names(boltz)) "Affinity_pIC50" else "Affinity_Pred_Value"
ok <- boltz %>% filter(Status == "OK")

ligand_order <- ok %>%
  group_by(Ligand) %>%
  summarize(mean_score = mean(.data[[score_col]], na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_score)) %>%
  pull(Ligand)

heat_dat <- boltz %>%
  mutate(
    Protein = factor(Protein, levels = c("SLC27A5", "ENO1")),
    Ligand = factor(Ligand, levels = ligand_order),
    Label = ifelse(is.na(.data[[score_col]]), "NA", sprintf("%.2f", .data[[score_col]]))
  )

p1 <- ggplot(heat_dat, aes(x = Ligand, y = Protein, fill = .data[[score_col]])) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = Label), size = 2.8, color = "white", fontface = "bold") +
  scale_fill_gradient2(
    low = "#4575B4",
    mid = "#FEE090",
    high = "#D73027",
    midpoint = mean(ok[[score_col]], na.rm = TRUE),
    na.value = "#D0D5DA"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "Model-gene Boltz-2 predicted pIC50",
    x = "Phytochemical ligand",
    y = "Model gene",
    fill = "pIC50"
  )

top_pairs <- ok %>%
  group_by(Protein) %>%
  slice_max(order_by = .data[[score_col]], n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(.data[[score_col]]))

p2 <- ggplot(top_pairs, aes(x = .data[[score_col]], y = Protein, size = Confidence, color = Protein)) +
  geom_point(alpha = 0.9) +
  geom_text(aes(label = Ligand), nudge_x = 0.02, size = 3.2, hjust = 0) +
  scale_size(range = c(4, 8)) +
  scale_color_manual(values = c("SLC27A5" = "#6FA8B5", "ENO1" = "#2F4F5F")) +
  xlim(min(top_pairs[[score_col]], na.rm = TRUE) * 0.95, max(top_pairs[[score_col]], na.rm = TRUE) * 1.2) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  ) +
  labs(
    title = "Best ligand per model gene",
    x = "Predicted pIC50 (higher is better)",
    y = NULL,
    size = "Confidence"
  )

fig <- plot_grid(p1, p2, labels = c("A", "B"), ncol = 1, rel_heights = c(2.2, 1))

png_out <- file.path(out_dir, "Supp_Figure_S8_ModelGene_Docking.png")
pdf_out <- file.path(out_dir, "Supp_Figure_S8_ModelGene_Docking.pdf")
ggsave(png_out, fig, width = 13, height = 10, dpi = 300, bg = "white")
ggsave(pdf_out, fig, width = 13, height = 10)

top_out <- file.path(base_dir, "results", "boltz2_modelgene_top_pairs.csv")
write.csv(top_pairs, top_out, row.names = FALSE)

message("Saved: ", png_out)
message("Saved: ", pdf_out)
message("Saved: ", top_out)
