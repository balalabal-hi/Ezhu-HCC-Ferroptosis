# Reproducibility guide

## Mode A: output audit (no re-run)

Use this mode to verify reported results from included outputs.

- Main figures: `HCC_Ferroptosis_Project/plots/publication/`
- Supplementary figures: `HCC_Ferroptosis_Project/plots/supplementary/`
- Derived tables: `HCC_Ferroptosis_Project/results/`

## Mode B: full pipeline re-run

1. Enter project directory:
- `cd HCC_Ferroptosis_Project`

2. Install required R packages:
- `Rscript scripts_final/00_setup_env.R`

3. Prepare raw/reference inputs:
- Follow `DATA_MANIFEST.md` and `docs/data-acquisition.md`.

4. Run end-to-end pipeline:
- `Rscript scripts_final/run_complete_pipeline.R`

Expected outputs:
- `results/`
- `plots/publication/`
- `plots/supplementary/`

## Optional utilities

- Regenerate Figure 1 workflow: `Rscript scripts_final/07b_fig1_jimr_workflow.R`
- Regenerate Figure 2 panel: `Rscript scripts_final/07a_fig2_jimr_deg_ferroptosis.R`
- Regenerate model-gene docking supplement: `Rscript scripts_final/12_modelgene_boltz2_supplement.R`

## Notes

- Docking affinity API calls are optional and require `NVIDIA_API_KEY` in the shell environment.
