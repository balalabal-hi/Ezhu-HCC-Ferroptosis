# Project Status (Ezhu–Ferroptosis–HCC)

Last updated: 2026-01-29

This folder is the "openspec" layer for the project: concise, human-readable statements of scope, datasets, and provenance rules that are easy to audit.

## Scope (what this project does and does not claim)

- Builds a ferroptosis-related prognostic signature for hepatocellular carcinoma (HCC) using public transcriptomic cohorts.
- Uses Ezhu (*Curcuma phaeocaulis*) target evidence + docking as a hypothesis-generating "actionable node" layer (not experimental proof of mechanism or efficacy).

## Cohorts used in the analysis

- Discovery / training: GEO `GSE14520`.
- External validation: `TCGA-LIHC`, GEO `GSE76427`, `GSE10143-HCC`, `GSE27150`.
- Supplementary validation: MVI cohort from `PMC8692135` (Table S2; file `ijbsv18p0261s2.csv`).

## Key robustness checks (Supplementary)

- Proportional hazards (PH) diagnostics (Schoenfeld residual test).
- RMST-based comparison at 60 months as a PH-robust complement.
- Random-signature sanity check (1000 iterations) to contextualize the observed C-index.

## Provenance rules (reference vs derived)

- `Ezhu_HCC_Project/data/references/` contains only curated reference inputs (e.g., FerrDb gene list snapshots, curated target lists, small metadata tables). These files should not be overwritten by analysis scripts.
- All computed outputs (DEG, risk scores, validation metrics, diagnostics, docking summaries) must be written to `Ezhu_HCC_Project/results/` and/or `Ezhu_HCC_Project/plots/`.

## AI-assisted Figure 1 disclosure

- The workflow schematic (Figure 1) was generated with assistance from Google Gemini ("Nano Banana" model) using a text prompt, then reviewed/edited by the authors.
- The current canonical prompt is `openspec/FIGURE1_PROMPT_EZHU_HCC.md` (a copy may be included in Supplementary Materials for submission).

