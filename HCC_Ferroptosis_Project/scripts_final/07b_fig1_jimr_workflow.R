#!/usr/bin/env Rscript

# 07b_fig1_jimr_workflow.R
# JIMR-friendly workflow schematic (no TCM / docking layer).
#
# Outputs:
# - plots/publication/JIMR_Figure1_workflow.pdf
# - plots/publication/JIMR_Figure1_workflow.png
#
# Usage:
#   Rscript scripts_final/07b_fig1_jimr_workflow.R

if (!interactive()) pdf(NULL)

suppressPackageStartupMessages({
  library(ggplot2)
})

detect_base_dir <- function() {
  candidates <- c(".", "Ezhu_HCC_Project")
  for (d in candidates) {
    if (dir.exists(file.path(d, "results")) && dir.exists(file.path(d, "plots")) && dir.exists(file.path(d, "data/processed"))) return(d)
  }
  stop("Cannot locate project root.")
}

base_dir <- detect_base_dir()
pub_dir <- file.path(base_dir, "plots/publication")
dir.create(pub_dir, showWarnings = FALSE, recursive = TRUE)

mk_box <- function(xc, yc, w, h, label, fill, border, text_size = 3.4, text_face = "plain") {
  list(
    rect = data.frame(
      xmin = xc - w / 2,
      xmax = xc + w / 2,
      ymin = yc - h / 2,
      ymax = yc + h / 2,
      fill = fill,
      border = border
    ),
    text = data.frame(
      x = xc,
      y = yc,
      label = label,
      size = text_size,
      face = text_face
    )
  )
}

mk_arrow <- function(x1, y1, x2, y2) {
  data.frame(x = x1, y = y1, xend = x2, yend = y2)
}

boxes <- list()

# Colors
fill_title <- "#E6ECEE"
fill_body <- "#F4F7F8"
border <- "#A8B6BA"
arrow_col <- "#6B7C81"

# 1. Data & Resources
boxes[[1]] <- mk_box(0.20, 0.935, 0.34, 0.08, "1. Data & Resources", fill_title, border, text_size = 4.2, text_face = "bold")
boxes[[2]] <- mk_box(
  0.20, 0.80, 0.34, 0.22,
  "Discovery:\nGEO GSE14520 (tumor vs non-tumor; OS)\n\nExternal validation:\nTCGA-LIHC; GEO GSE76427; GSE10143-HCC; GSE27150\n\nFerroptosis reference:\nFerrDb (curated)",
  fill_body, border, text_size = 3.2
)

# 2. HCC-context landscape
boxes[[3]] <- mk_box(0.73, 0.935, 0.54, 0.08, "2. HCC-context Ferroptosis Landscape", fill_title, border, text_size = 4.2, text_face = "bold")
boxes[[4]] <- mk_box(
  0.73, 0.80, 0.54, 0.22,
  "Differential expression (limma):\ntumor vs non-tumor in GSE14520\n\nHCC-context ferroptosis candidates:\nFerrDb genes filtered by DEG (FDR < 0.05)\n→ 125 candidates",
  fill_body, border, text_size = 3.2
)

# 3. Prognostic signature
boxes[[5]] <- mk_box(0.50, 0.655, 0.90, 0.08, "3. Prognostic Signature (Discovery → Validation)", fill_title, border, text_size = 4.2, text_face = "bold")
boxes[[6]] <- mk_box(
  0.50, 0.51, 0.90, 0.20,
  "Univariate Cox screening → LASSO-Cox (10-fold CV) → 9-gene signature\n\nPerformance: KM survival; time-dependent ROC (1/3/5 years); bootstrap validation\n\nExternal validation: TCGA-LIHC and independent GEO cohorts",
  fill_body, border, text_size = 3.2
)

# 4. Clinical utility
boxes[[7]] <- mk_box(0.20, 0.385, 0.34, 0.08, "4. Clinical Utility", fill_title, border, text_size = 4.2, text_face = "bold")
boxes[[8]] <- mk_box(
  0.20, 0.26, 0.34, 0.16,
  "Univariate & multivariate Cox\n\nNomogram:\nrisk score + available clinical variables\n\nDecision curve analysis (3-year)",
  fill_body, border, text_size = 3.2
)

# 5. Immune context
boxes[[9]] <- mk_box(0.73, 0.385, 0.54, 0.08, "5. Immune Context (Exploratory)", fill_title, border, text_size = 4.2, text_face = "bold")
boxes[[10]] <- mk_box(
  0.73, 0.26, 0.54, 0.16,
  "Immune infiltration:\nssGSEA (immune signatures)\n\nImmune checkpoints:\ndifferential expression (High vs Low risk)\n\nTherapeutic response inference:\nTIDE (exploratory)",
  fill_body, border, text_size = 3.2
)

# 6. Calibration & robustness
boxes[[11]] <- mk_box(0.50, 0.115, 0.90, 0.08, "6. Calibration & Robustness", fill_title, border, text_size = 4.2, text_face = "bold")
boxes[[12]] <- mk_box(
  0.50, 0.035, 0.90, 0.10,
  "Calibration: internal bootstrap; external 3-year calibration (discovery baseline hazard)\nPrediction error: Brier score (IPCW)\nDiagnostics: PH test; RMST (60 months); random-signature sanity check",
  fill_body, border, text_size = 3.2
)

rect_df <- do.call(rbind, lapply(boxes, `[[`, "rect"))
text_df <- do.call(rbind, lapply(boxes, `[[`, "text"))

arrows <- rbind(
  mk_arrow(0.37, 0.80, 0.46, 0.80),
  mk_arrow(0.73, 0.69, 0.60, 0.60),
  mk_arrow(0.50, 0.41, 0.37, 0.34),
  mk_arrow(0.50, 0.41, 0.63, 0.34),
  mk_arrow(0.50, 0.41, 0.50, 0.16)
)

p <- ggplot() +
  geom_rect(
    data = rect_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = rect_df$fill,
    color = rect_df$border,
    linewidth = 0.6
  ) +
  geom_segment(
    data = arrows,
    aes(x = x, y = y, xend = xend, yend = yend),
    color = arrow_col,
    linewidth = 0.8,
    arrow = grid::arrow(length = grid::unit(0.02, "npc"), type = "closed")
  ) +
  geom_text(
    data = text_df,
    aes(x = x, y = y, label = label),
    size = text_df$size,
    fontface = text_df$face,
    lineheight = 1.05
  ) +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", color = NA))

ggsave(file.path(pub_dir, "JIMR_Figure1_workflow.pdf"), p, width = 13.33, height = 7.5, device = cairo_pdf)
ggsave(file.path(pub_dir, "JIMR_Figure1_workflow.png"), p, width = 13.33, height = 7.5, dpi = 300, bg = "white")
message("✅ Saved: plots/publication/JIMR_Figure1_workflow (PDF+PNG)")
