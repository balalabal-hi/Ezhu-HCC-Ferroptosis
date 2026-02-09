# HCC-Ferroptosis-Signature

This repository provides a reproducible analysis pipeline for a ferroptosis-related prognostic signature in hepatocellular carcinoma (HCC), including multi-cohort validation, immune-context analysis, therapeutic response inference, and exploratory docking context.

## Repository structure

- `HCC_Ferroptosis_Project/scripts_final/`: end-to-end analysis scripts (R + small Python utilities).
- `HCC_Ferroptosis_Project/results/`: derived result tables used in the manuscript.
- `HCC_Ferroptosis_Project/plots/publication/`: main figure panels.
- `HCC_Ferroptosis_Project/plots/supplementary/`: supplementary figure panels.
- `HCC_Ferroptosis_Project/data/references/`: curated reference inputs required by the pipeline.
- `docs/`: reproducibility, data acquisition, and repository scope notes.

Large raw downloads (GEO/TCGA) and large third-party matrices are intentionally not versioned. See `docs/data-acquisition.md`.

## Quick start

1. Review outputs directly:
   - `HCC_Ferroptosis_Project/results/`
   - `HCC_Ferroptosis_Project/plots/publication/`
2. Re-run from source data:
   - `docs/reproducibility.md`
3. Check expected input layout:
   - `HCC_Ferroptosis_Project/DATA_MANIFEST.md`
   - `HCC_Ferroptosis_Project/data/references/README_DATA_SOURCES.md`

## Notes

- Some filenames and column names retain legacy suffixes (for backward compatibility with earlier scripts), but the current analysis branch is the authoritative version.

## License

MIT License (`LICENSE`).
