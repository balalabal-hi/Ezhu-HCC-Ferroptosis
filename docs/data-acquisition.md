# Data acquisition (what to download, and where)

This repository intentionally does not redistribute large raw downloads (GEO/TCGA) or large third‑party pharmacogenomic matrices. The analysis scripts expect those resources to be prepared under `Ezhu_HCC_Project/` using the same relative paths as documented below.

For the authoritative “what goes where” manifest (including expected filenames), see:
- `Ezhu_HCC_Project/DATA_MANIFEST.md`
- `Ezhu_HCC_Project/data/references/README_DATA_SOURCES.md`

## GEO cohorts

Primary discovery cohort:
- **GSE14520** (tumor vs non-tumor; survival available)

External validation cohorts:
- **GSE76427**
- **GSE10143**
- **GSE27150** (two-channel microarray; expression is a log2 ratio and should be interpreted cautiously)

Recommended workflow:
- Follow the manifest in `Ezhu_HCC_Project/DATA_MANIFEST.md`.
- Run `Ezhu_HCC_Project/scripts_final/01b_download_GSE14520.R` to download/parse GSE14520 into `Ezhu_HCC_Project/data/processed/`.
- Run `Ezhu_HCC_Project/scripts_final/02e_external_validation.R` to assemble external validation inputs (some cohorts may still require manual preparation depending on GEO access and platform files).

If you prefer manual downloads (e.g., restricted network), the manifest documents the expected landing paths, including:
- `Ezhu_HCC_Project/data/raw/GSE14520/GSE14520-GPL3921_series_matrix.txt.gz`
- `Ezhu_HCC_Project/data/raw/GSE14520/GSE14520-GPL571_series_matrix.txt.gz` (optional, if used)
- `Ezhu_HCC_Project/data/raw/GSE14520/GPL571.annot.gz` (platform annotation required by immune scripts)

## TCGA-LIHC

The pipeline supports TCGA-LIHC as an external validation cohort.
- If you already have processed objects, place them under `Ezhu_HCC_Project/data/processed/` as described in `Ezhu_HCC_Project/DATA_MANIFEST.md`.
- Otherwise, obtain the data via the original provider/workflow (e.g., GDC portal / TCGAbiolinks) and convert into the expected RDS formats.

## Reference inputs

Small curated reference tables required for the pipeline are included under:
- `Ezhu_HCC_Project/data/references/`

For provenance and required columns, see:
- `Ezhu_HCC_Project/data/references/README_DATA_SOURCES.md`

## Drug sensitivity (GDSC2)

GDSC2 expression/response matrices can be large and are not redistributed here.
- See `Ezhu_HCC_Project/DATA_MANIFEST.md` for the expected path under `Ezhu_HCC_Project/data/references/GDSC/`.
- The paper‑level derived outputs (IC50 predictions and downstream summaries) are already included under `Ezhu_HCC_Project/results/`.
