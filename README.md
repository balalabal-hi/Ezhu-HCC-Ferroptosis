# Ezhu–Ferroptosis–HCC (Reproducible Analysis Package)

This repository contains the analysis scripts and derived outputs for the study of *Curcuma phaeocaulis* (Ezhu) target-guided ferroptosis landscape and prognostic modeling in hepatocellular carcinoma (HCC).

## What is in this repository

- `Ezhu_HCC_Project/scripts_final/`: end-to-end analysis scripts (R + a small amount of Python for docking task preparation)
- `Ezhu_HCC_Project/results/`: derived tables used in the manuscript (model coefficients, risk scores, external validation summaries, docking summaries, etc.)
- `Ezhu_HCC_Project/plots/publication/`: final Figure 2–6 panels (PDF/PNG/TIFF as available)
- `Ezhu_HCC_Project/data/references/`: small, curated reference inputs required by the pipeline (gene lists, target lists, docking ligand metadata)
- `docs/`: notes on data acquisition and reproducibility

Raw downloads (e.g., GEO series matrix files, TCGA raw downloads, large pharmacogenomic matrices) are not included due to size and/or redistribution constraints. See `docs/data-acquisition.md`.

## Quick orientation

- For a read-only verification of the main outputs, start from:
  - `Ezhu_HCC_Project/plots/publication/`
  - `Ezhu_HCC_Project/results/`
- For a full rerun (network + compute required), see `docs/reproducibility.md`.
- For dataset placement and expected filenames, use:
  - `Ezhu_HCC_Project/DATA_MANIFEST.md`
  - `Ezhu_HCC_Project/data/references/README_DATA_SOURCES.md`

## License

Code in this repository is released under the MIT License (see `LICENSE`).
