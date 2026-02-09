# Preflight checklist (before pushing to GitHub)

## Security

- Confirm no `.env` is tracked and no API keys appear in the repository history.
- Run a quick secret scan (examples):
  - `rg -n "API_KEY|SECRET|TOKEN|PASSWORD|sk-|AIza" -S .`

## Reproducibility

- Confirm the manifest files are present:
  - `HCC_Ferroptosis_Project/DATA_MANIFEST.md`
  - `HCC_Ferroptosis_Project/data/references/README_DATA_SOURCES.md`
- Confirm key derived tables exist under `HCC_Ferroptosis_Project/results/` (risk score, model coefficients, external validation).
- Confirm publication figures exist under `HCC_Ferroptosis_Project/plots/publication/`.

## Size / licensing

- Do not push large raw downloads (`data/raw/`, `GDCdata/`), rebuilt intermediates (`data/processed/`), or large third-party matrices (e.g., GDSC).
- If you need to share large archives, use GitHub Releases or an archival repository (Zenodo/OSF) and link from `README.md`.
