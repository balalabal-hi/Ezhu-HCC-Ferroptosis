# Figure 1 generation prompt (Ezhu–HCC)

This prompt was used as a starting point to generate a vector-style workflow schematic (Figure 1) consistent with the analysis pipeline and the terminology used in Figures 2–6.

## Prompt (English, copy-paste)

Create a clean, publication-quality **vector-style workflow schematic (Figure 1)** for an SCI paper.  
Topic: *Curcuma phaeocaulis* (Ezhu)–anchored ferroptosis landscape and prognostic modeling in hepatocellular carcinoma (HCC).

### Style requirements (match modern SCI figures)
- Flat vector infographic, **white background**, subtle light-gray separators, minimal shadows.
- Color palette: low-saturation **blue–teal** with neutral grays; avoid bright colors.
- Consistent sans-serif font (Helvetica/Arial-like), high readability, no decorative fonts.
- Rounded rectangles, thin lines, clear arrows, consistent spacing.
- No logos, no photos, no 3D, no cartoon style.
- Layout: **16:9** wide canvas, high resolution (e.g., 4800×2700). Export-ready.

### Figure structure (left-to-right, modular, with arrows)
1) **Module 1 — Data & Resources**
   - GEO: **GSE14520** (HCC; tumor vs non-tumor; OS available)
   - External validation: **TCGA-LIHC**, **GSE76427**, **GSE10143-HCC**
   - Supplement: **MVI cohort (PMC8692135, Table S2: ijbsv18p0261s2.csv)**
   - Ferroptosis genes: **FerrDb V2**
   - Ezhu target evidence: **HERB (reference mining + statistical inference) + ChEMBL (known targets)**

2) **Module 2 — HCC-context Ferroptosis Landscape**
   - Differential expression (limma): tumor vs non-tumor in GSE14520
   - Build **HCC-context ferroptosis** list: FerrDb filtered by DEG (FDR<0.05)
   - Define two gene sets (show as a callout box with two bullet points):
     - **Intersection anchors (strict 3-way)**: DEG (FDR<0.05 & |log2FC|>1) ∩ HCC-context ferroptosis ∩ Ezhu targets → **5 genes**
     - **Expanded candidate panel (HERB-supported)**: DEG (FDR<0.05) ∩ HCC-context ferroptosis ∩ HERB-supported Ezhu targets → **11 genes**

3) **Module 3 — Prognostic Signature (Discovery → Validation)**
   - From “ferroptosis landscape genes” → univariate Cox (P<0.10) → **LASSO-Cox (10-fold CV, lambda.min)** → **9-gene signature**
   - Outputs: risk score, KM survival, time-dependent ROC (1/3/5-year), bootstrap validation (n=500)
   - External validation: TCGA-LIHC + GEO cohorts (report as “external validation” without claiming strong performance)

4) **Module 4 — Clinical Utility**
   - Univariate & multivariate Cox
   - **Nomogram**: risk score + tumor size + AFP + TNM group
   - **DCA (3-year; threshold 0–0.6)** and **C-index comparison**

5) **Module 5 — Tumor Microenvironment & Therapeutic Clues**
   - Immune infiltration: **ssGSEA (28 immune cell types)** → correlation & group differences
   - Immune checkpoints: differential expression (log2FC High–Low risk)
   - Drug sensitivity: **GDSC2-based predicted IC50** (batch adjusted by ComBat) → correlation & group comparison

6) **Module 6 — Structural Evidence (In silico)**
   - **CB-Dock2 / AutoDock Vina** blind docking scores
   - **Boltz-2 AI affinity** (note: “Affinity score: lower is better; NA = no valid output after retries”)
   - Consistency check: Spearman correlation between -Vina (higher is better) vs AI affinity (lower is better)

### Design notes
- Put a title at top: “Ezhu-anchored ferroptosis landscape and prognostic modeling in HCC”.
- Use simple icons (database cylinder, gene list, Kaplan–Meier curve, ROC curve, forest plot, immune cell, pill, protein-ligand).
- Ensure all text is spelled correctly and concise; avoid any claims like “proved therapeutic effect”.
- **Do NOT include** GO/KEGG, WGCNA, CHB datasets, network pharmacology topology, or any unused modules.

## Checklist
- Datasets mentioned: GSE14520, TCGA-LIHC, GSE76427, GSE10143-HCC, PMC8692135 S2 (MVI supplement).
- Immune cell types: **28**; strict intersection: **5**; expanded panel: **11**; signature: **9-gene**.
- Avoid any modules not used in the main pipeline.

