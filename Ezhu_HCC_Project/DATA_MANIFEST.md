# Data Manifest（手动下载与放置路径说明）

> 目的：把“论文实际研究依赖了哪些外部数据、哪些是手动下载、放到哪里、被哪些脚本读取”写清楚，方便后续整理对外复现材料。
>
> 注意：由于体积/访问限制，本项目在本地可能包含部分已处理的 `.rds` / `.csv`，但对外复现包不一定包含这些文件；因此建议把“**需要手动准备**”与“**脚本会生成**”严格区分。

---

## 0) 两种复现模式（建议对审稿人明确）

### Mode A：对账式复核（不重新计算）

适用：审稿人只需要确认“图版/关键数字”来源清晰、文件可追溯。  
所需文件：`results/` + `plots/publication/`（以及参考基因表/靶点表）。

- 图版目录：`plots/publication/`
- 关键结果目录：`results/`
- 最终图板脚本：`scripts_final/07_final_publication_figures.R`

### Mode B：全流程重跑（从原始数据再生成）

适用：需要从 GEO/TCGA 等源头重算中间结果（更耗时、对网络/环境要求更高）。  
除 Mode A 外，还需要准备下述“原始/参考数据”。

---

## 1) GEO 原始数据（通常需要手动下载后放置）

### 1.1 GSE14520（HCC 主队列）

脚本：`scripts_final/01b_download_GSE14520.R`

需要手动准备（放到 `data/raw/GSE14520/`）：
- `GSE14520-GPL3921_series_matrix.txt.gz`（主平台；脚本优先解析）
- `GSE14520-GPL571_series_matrix.txt.gz`（次平台；若存在则一并解析）

说明：
- 这些文件可从 NCBI GEO 的 GSE14520 页面下载（Series Matrix File(s)）。
- 若网络受限，可先人工下载，再放置到上述路径，然后运行脚本解析为 `data/processed/*.rds`。

### 1.2 GSE83148（CHB，已移出主线）

CHB 相关队列已从论文主线移除，不再作为复现必需数据。

### 1.3 外部验证 GEO 队列（主线使用）

脚本：`scripts_final/02e_external_validation.R`

需要联网下载或手动准备：
- `GSE76427`（外部验证，已支持自动下载与解析）
- `GSE10143`（外部验证，需准备 processed 文件）
- `GSE27150`（外部验证；two-channel 平台，表达量为 log2 比值；脚本支持自动下载 + GPL 映射，但需对结果口径保持谨慎）

---

## 2) 平台注释（会被免疫脚本读取，通常需要手动准备）

脚本：`scripts_final/03a_immune_infiltration.R`

需要手动准备（放到 `data/raw/GSE14520/`）：
- `GPL571.annot.gz`

用途：
- 将 GSE14520（GPL571 平台的探针 ID）映射到基因符号，生成 `data/processed/GSE14520_expr_symbol.rds`。

---

## 3) TCGA-LIHC（外部验证；可能需要手动/联网）

本地已有（可直接用于复跑外部验证脚本）：
- `data/processed/TCGA_LIHC_expr.rds`
- `data/processed/TCGA_LIHC_clinical.rds`

说明：
- 若对外复现包不包含上述文件，则必须明确：需要联网通过 `TCGAbiolinks` 下载，或从 GDC 门户手动下载并转换为脚本可读格式（对外包需提供具体落点与转换脚本）。

---

## 4) 参考数据（通常需要手动准备/整理）

### 4.1 莪术靶点（HERB + ChEMBL）

文件：`data/references/tcm_targets_ezhu.csv`（最低需要 `Target,Source` 两列）

可选补充明细：
- `data/references/tcm_targets_ezhu_chembl.csv`（ChEMBL 已知靶点扩展明细）

获取指南：`data/references/README_DATA_SOURCES.md`

### 4.2 铁死亡基因集

参考基因池（用于溯源/审计，“分母”）：
- `data/references/ferroptosis_genes_expanded.csv`（FerrDb V2 整理后的静态列表）

衍生结果（用于主线分析/作图，“分子”）：
- `results/ferroptosis_genes_hcc_context.csv`（FerrDb ∩ GSE14520-DEG，含 logFC/adjP）
  - 生成脚本：`scripts_final/00c_prepare_hcc_ferro_genes.R`

说明：
- 主线脚本默认优先使用 `results/ferroptosis_genes_hcc_context.csv`；若缺失则回退到参考池 `data/references/ferroptosis_genes_expanded.csv`。
- 参考池不应被分析脚本写入覆盖。

---

## 5) 药物敏感性（GDSC2；通常需要手动准备）

脚本：`scripts_final/05c_drug_sensitivity.R`

需要准备（放到 `data/references/GDSC/`）：
- `GDSC2_Expr.rds`
- `GDSC2_Res.rds`

本地还存在：
- `data/references/GDSC/DataFiles.zip`（体积很大；对外包通常不直接包含）

---

## 6) 分子对接（存在“手动平台结果”与“API 计算结果”两类）

### 6.1 CB-Dock2（手动在线平台）

原始任务与结果（本地保留）：
- `data/references/docking/CB_DOCK2_TASKS.md`
- `data/references/docking/cbdock2_results/`（平台输出的 txt）

论文/图表使用的汇总结果：
- `results/cbdock2_docking_results.csv`

配体清单（用于对接准备）：
- `data/references/docking/ligands_smiles_ezhu.csv`

### 6.2 Boltz-2（NVIDIA API；需要 API Key）

论文使用的最终结果（由外部平台/脚本产出）：
- `results/molecular_docking/boltz2_docking_results_final.csv`

---

## 7) 最小检查清单（建议对外材料也保持一致）

- [ ] `results/risk_score_data_ezhu.csv` 存在（并作为免疫/药敏/外部验证的唯一风险分群来源）
- [ ] `plots/publication/` 下 Figure 1–6 文件齐全
- [ ] `data/references/tcm_targets_ezhu.csv` 与 `data/references/ferroptosis_genes_expanded.csv` 存在（用于溯源/审计参考池）
- [ ] `results/ferroptosis_genes_hcc_context.csv` 存在（用于主线 hub genes 交集与图板动态计算）
