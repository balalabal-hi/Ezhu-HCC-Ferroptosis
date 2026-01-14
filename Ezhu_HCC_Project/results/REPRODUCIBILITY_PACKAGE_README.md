# Reproducibility Package (Ezhu–Ferroptosis–HCC)

This repository contains analysis scripts and derived results for the Ezhu (Curcuma phaeocaulis) – ferroptosis – HCC study.

## What is included in the reproducibility ZIP

The ZIP package (generated at the project root) is intended for journal submission and peer review. It contains:

- `Ezhu_HCC_Project/scripts_final/`: end-to-end scripts used in the current pipeline
- `Ezhu_HCC_Project/results/`: key result tables used in the manuscript (model stats, external validation, docking/AI outputs, etc.)
- `Ezhu_HCC_Project/plots/publication/`: final Figure 2–6 outputs (PDF/PNG/TIFF)
- `Ezhu_HCC_Project/plots/supplementary/`: supplementary plots (e.g., external validation panels)
- `Ezhu_HCC_Project/data/processed/`: processed expression/clinical objects used to reproduce the results
- `Ezhu_HCC_Project/data/references/`: curated reference inputs used by the pipeline (excluding large third-party pharmacogenomic matrices)

## What is NOT redistributed (and why)

- `Ezhu_HCC_Project/data/references/GDSC/` is intentionally excluded because it is large and should be obtained from the original provider/source.

Drug sensitivity results needed for the paper are already included as derived outputs under `Ezhu_HCC_Project/results/`.

## How to reproduce main results (Figure 2–6)

From the project root:

1. Run the pipeline:
   - `Rscript Ezhu_HCC_Project/scripts_final/run_complete_pipeline.R`

2. Re-generate final Figure 2–6 panels:
   - `Rscript Ezhu_HCC_Project/scripts_final/07_final_publication_figures.R`

## Notes on API keys

- `scripts_final/09_boltz2_affinity.py` requires `NVIDIA_API_KEY` to be set in the environment.
- API keys must never be embedded in scripts or shared archives.

