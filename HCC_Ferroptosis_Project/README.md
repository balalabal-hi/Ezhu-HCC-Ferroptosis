# HCC Ferroptosis Signature Project

## Overview
This directory contains the finalized reproducible pipeline for a ferroptosis-related prognostic signature in hepatocellular carcinoma (HCC), including:
- model development in discovery data,
- multi-cohort external validation,
- immune-context analysis,
- therapeutic response inference,
- exploratory docking context.

## Directory layout
- `scripts_final/`: finalized analysis scripts.
- `results/`: derived tables used by the manuscript.
- `plots/publication/`: main figure outputs.
- `plots/supplementary/`: supplementary figure outputs.
- `data/references/`: curated reference inputs.
- `DATA_MANIFEST.md`: required input/output map.

## Reproducibility
Run from `HCC_Ferroptosis_Project/`:

```bash
Rscript scripts_final/00_setup_env.R
Rscript scripts_final/run_complete_pipeline.R
```

Optional figure-only refresh:

```bash
Rscript scripts_final/07b_fig1_jimr_workflow.R
Rscript scripts_final/07a_fig2_jimr_deg_ferroptosis.R
Rscript scripts_final/07_final_publication_figures.R
Rscript scripts_final/12_modelgene_boltz2_supplement.R
```

## Cohorts in the finalized branch
- Discovery: `GSE14520`
- External validation: `TCGA-LIHC`, `GSE76427`, `GSE10143-HCC`, `GSE27150`, `ICGC-LIRI-JP (HCCDB18)`

## Notes
- Some file names retain legacy suffixes (for backward compatibility with earlier script interfaces).
- The authoritative execution entrypoint is `scripts_final/run_complete_pipeline.R`.

## License
MIT License (`../LICENSE`).
