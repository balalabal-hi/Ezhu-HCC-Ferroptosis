# Repository scope (tracked vs omitted)

## Scientific scope

- Transcriptome-based prognostic modeling in hepatocellular carcinoma (HCC).
- Ferroptosis-related signature development and multi-cohort validation.
- Immune-context and therapeutic-response analyses.
- Exploratory docking context as hypothesis-generating support.

## Cohorts

- Discovery: GEO `GSE14520`.
- External validation: `TCGA-LIHC`, `GSE76427`, `GSE10143-HCC`, `GSE27150`, `ICGC-LIRI-JP (HCCDB18)`.

## Robustness checks

- PH diagnostics (Schoenfeld tests).
- RMST at 60 months as a PH-robust complement.
- Random-signature sanity check.
- External calibration and IPCW Brier scores.
- Incremental value over available clinical covariates in TCGA.
- Exploratory TCGA multi-omics characterization.

## Provenance rule

- Curated reference inputs: `HCC_Ferroptosis_Project/data/references/` (immutable inputs).
- Derived outputs: `HCC_Ferroptosis_Project/results/` and `HCC_Ferroptosis_Project/plots/`.

## Tracked in Git

- `HCC_Ferroptosis_Project/scripts_final/`
- `HCC_Ferroptosis_Project/results/` (derived tables)
- `HCC_Ferroptosis_Project/plots/` (rendered figures)
- `HCC_Ferroptosis_Project/data/references/` (small curated references)
- `HCC_Ferroptosis_Project/DATA_MANIFEST.md`
- `docs/`

## Omitted from Git

- Large raw downloads: `HCC_Ferroptosis_Project/data/raw/`, `HCC_Ferroptosis_Project/GDCdata/`
- Large third-party matrices: `HCC_Ferroptosis_Project/data/references/GDSC/`
- Local curation dumps: `HCC_Ferroptosis_Project/data/references/raw/`
- Local API response caches: `HCC_Ferroptosis_Project/results/boltz2_responses/`, `HCC_Ferroptosis_Project/results/boltz2_failures/`

See `docs/data-acquisition.md` and `HCC_Ferroptosis_Project/DATA_MANIFEST.md` for reproducible input placement.
