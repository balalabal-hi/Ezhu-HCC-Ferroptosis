# Reproducibility Package Notes

## Included in the repository

- `HCC_Ferroptosis_Project/scripts_final/`: finalized analysis scripts.
- `HCC_Ferroptosis_Project/results/`: derived result tables for reporting.
- `HCC_Ferroptosis_Project/plots/publication/`: main figure outputs.
- `HCC_Ferroptosis_Project/plots/supplementary/`: supplementary figure outputs.
- `HCC_Ferroptosis_Project/data/processed/`: processed cohort objects required by scripts.
- `HCC_Ferroptosis_Project/data/references/`: curated reference inputs.

## Not redistributed in Git

- Large third-party pharmacogenomic matrices under `data/references/GDSC/`.
- Large raw downloads under `data/raw/` and `GDCdata/`.

Downstream derived outputs depending on these resources are included under `results/`.

## Re-run instructions

From repository root:

1. `cd HCC_Ferroptosis_Project`
2. `Rscript scripts_final/00_setup_env.R`
3. `Rscript scripts_final/run_complete_pipeline.R`

To refresh figure panels only:

- `Rscript scripts_final/07_final_publication_figures.R`
- `Rscript scripts_final/07a_fig2_jimr_deg_ferroptosis.R`
- `Rscript scripts_final/07b_fig1_jimr_workflow.R`
- `Rscript scripts_final/12_modelgene_boltz2_supplement.R`

## Optional docking affinity calls

- `scripts_final/09_boltz2_affinity.py` requires `NVIDIA_API_KEY` in the environment.
