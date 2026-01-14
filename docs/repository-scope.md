# Repository scope (what is tracked vs omitted)

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

