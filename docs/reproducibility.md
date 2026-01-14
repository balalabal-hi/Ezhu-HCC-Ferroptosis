# Reproducibility notes

Two practical ways to use this repository are:

## Mode A: Audit / read-only verification

Use this mode to inspect the exact derived tables and figure panels used for reporting.

- Main figures: `Ezhu_HCC_Project/plots/publication/`
- Derived tables: `Ezhu_HCC_Project/results/`
- Figure-to-file mapping: `Ezhu_HCC_Project/FIGURE_LIST.md`

## Mode B: Full rerun (from downloads)

This mode re-runs the pipeline end-to-end and regenerates the figure panels.

1) Enter the project directory:
- `cd Ezhu_HCC_Project`

2) Install required R packages:
- `Rscript scripts_final/00_setup_env.R`

3) Prepare required raw/reference inputs:
- Follow `DATA_MANIFEST.md` and `docs/data-acquisition.md`.

4) Run the full pipeline:
- `Rscript scripts_final/run_complete_pipeline.R`

Expected outputs:
- `results/` (tables)
- `plots/publication/` (Figure 2–6)
- `plots/supplementary/` (supplementary panels)

