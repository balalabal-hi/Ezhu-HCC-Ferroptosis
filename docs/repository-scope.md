# Repository scope (what is tracked vs omitted)

## Scientific scope (high-level)

- This repository supports a transcriptome-based analysis in hepatocellular carcinoma (HCC), centered on a ferroptosis-related prognostic signature derived from public cohorts.
- The Ezhu (*Curcuma phaeocaulis*) target/docking layer is presented as hypothesis-generating (pharmacologically actionable nodes), and is not framed as experimental proof of mechanism or efficacy.

## Cohorts used

- Discovery / training: GEO `GSE14520`.
- External validation: `TCGA-LIHC`, GEO `GSE76427`, `GSE10143-HCC`, `GSE27150`.
- Supplementary validation: MVI cohort from `PMC8692135` (Table S2; file `ijbsv18p0261s2.csv`).

## Robustness checks (Supplementary)

- Cox proportional hazards (PH) diagnostics (Schoenfeld residual test).
- RMST-based comparison at 60 months as a PH-robust complement.
- Random-signature sanity check (1000 iterations) to contextualize the observed C-index.

## Provenance rule (reference vs derived)

- Curated reference inputs live under `Ezhu_HCC_Project/data/references/` and should not be overwritten by analysis scripts.
- All computed outputs (DEG, risk scores, validation metrics, diagnostics, docking summaries) are written to `Ezhu_HCC_Project/results/` and/or `Ezhu_HCC_Project/plots/`.

## Figure 1 prompt disclosure

- Figure 1 (workflow schematic) was generated with assistance from Google Gemini ("Nano Banana" model) based on a text prompt, then reviewed and edited by the authors.
- The prompt used for Figure 1 is documented in `docs/figure1-generation-prompt.md` (a copy may also be provided as Supplementary material during submission).

## Tracked in Git

- Analysis scripts: `Ezhu_HCC_Project/scripts_final/`
- Derived outputs used for reporting: `Ezhu_HCC_Project/results/`
- Figure panels: `Ezhu_HCC_Project/plots/`
- Small curated reference inputs: `Ezhu_HCC_Project/data/references/`
- Documentation: `README.md`, `Ezhu_HCC_Project/DATA_MANIFEST.md`, `docs/`

## Omitted from Git

The following are intentionally not tracked here:

- Large raw inputs (GEO/TCGA downloads): `Ezhu_HCC_Project/data/raw/`, `Ezhu_HCC_Project/GDCdata/`
- Large third‑party matrices (GDSC2): `Ezhu_HCC_Project/data/references/GDSC/`
- Local exports/snapshots used during data curation: `Ezhu_HCC_Project/data/references/raw/`
- Optional API response archives: `Ezhu_HCC_Project/results/boltz2_responses/`

See `docs/data-acquisition.md` for how to obtain missing inputs, and `Ezhu_HCC_Project/DATA_MANIFEST.md` for the expected file locations.
