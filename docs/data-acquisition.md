# Data acquisition (what to download, and where)

This repository intentionally does not redistribute:
- large raw downloads (GEO/TCGA),
- rebuilt intermediate objects (`data/processed/`),
- large third-party pharmacogenomic matrices.

All required source datasets can be downloaded from official portals and then prepared under `HCC_Ferroptosis_Project/` using the paths below.

For the authoritative “what goes where” manifest (including expected filenames), see:
- `HCC_Ferroptosis_Project/DATA_MANIFEST.md`
- `HCC_Ferroptosis_Project/data/references/README_DATA_SOURCES.md`

## GEO cohorts

Primary discovery cohort:
- **GSE14520** (tumor vs non-tumor; survival available)
- GEO page: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE14520

External validation cohorts:
- **GSE76427**
- **GSE10143**
- **GSE27150** (two-channel microarray; expression is a log2 ratio and should be interpreted cautiously)
- GEO pages:
  - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76427
  - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10143
  - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE27150

Recommended workflow:
- Follow the manifest in `HCC_Ferroptosis_Project/DATA_MANIFEST.md`.
- Run `HCC_Ferroptosis_Project/scripts_final/01b_download_GSE14520.R` to download/parse GSE14520 into `HCC_Ferroptosis_Project/data/processed/`.
- Run `HCC_Ferroptosis_Project/scripts_final/02e_external_validation.R` to assemble external validation inputs (some cohorts may still require manual preparation depending on GEO access and platform files).

If you prefer manual downloads (e.g., restricted network), the manifest documents the expected landing paths, including:
- `HCC_Ferroptosis_Project/data/raw/GSE14520/GSE14520-GPL3921_series_matrix.txt.gz`
- `HCC_Ferroptosis_Project/data/raw/GSE14520/GSE14520-GPL571_series_matrix.txt.gz` (optional, if used)
- `HCC_Ferroptosis_Project/data/raw/GSE14520/GPL571.annot.gz` (platform annotation required by immune scripts)

## TCGA-LIHC

The pipeline supports TCGA-LIHC as an external validation cohort.
- GDC project page: https://portal.gdc.cancer.gov/projects/TCGA-LIHC
- If you already have processed objects, place them under `HCC_Ferroptosis_Project/data/processed/` as described in `HCC_Ferroptosis_Project/DATA_MANIFEST.md`.
- Otherwise, obtain the data via the original provider/workflow (e.g., GDC portal / TCGAbiolinks) and convert into the expected RDS formats.

## ICGC-LIRI-JP / HCCDB

- HCCDB portal (for HCCDB18 package): http://lifeome.net/database/hccdb/home.html

## Reference inputs

Small curated reference tables required for the pipeline are included under:
- `HCC_Ferroptosis_Project/data/references/`

For provenance and required columns, see:
- `HCC_Ferroptosis_Project/data/references/README_DATA_SOURCES.md`

## Drug sensitivity (GDSC2)

GDSC2 expression/response matrices can be large and are not redistributed here.
- GDSC portal: https://www.cancerrxgene.org/
- See `HCC_Ferroptosis_Project/DATA_MANIFEST.md` for the expected path under `HCC_Ferroptosis_Project/data/references/GDSC/`.
- The paper-level derived outputs (IC50 predictions and downstream summaries) are already included under `HCC_Ferroptosis_Project/results/`.

## FerrDb reference

- FerrDb portal: http://www.zhounan.org/ferrdb/current/
