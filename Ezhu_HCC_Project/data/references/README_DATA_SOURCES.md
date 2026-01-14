# 参考数据获取指南

本项目需要从权威数据库获取以下参考数据。**禁止使用虚拟数据或硬编码数据**。

---

## 1. 莪术 (Curcumae Rhizoma) 靶点数据 (tcm_targets_ezhu.csv)

### 数据来源选项

#### 选项A: ETCM数据库 (推荐)
1. 访问 http://www.tcmip.cn/ETCM/
2. 搜索 "莪术" 或 "Curcuma phaeocaulis" 或 "Ezhu"
3. 下载活性成分及其靶点数据
4. 保留可追溯的成分-靶点关系

#### 选项B: TCMSP数据库 (可选补充)
1. 访问 https://tcmsp-e.com/
2. 搜索 "莪术" 或 "Curcuma phaeocaulis"
3. 下载活性成分及其靶点数据
4. 筛选条件: OB ≥ 30%, DL ≥ 0.18

#### 选项C: HERB数据库
1. 访问 http://herb.ac.cn/
2. 搜索 "莪术" 或 "Curcuma phaeocaulis"
3. 下载成分-靶点数据

#### 选项D: ChEMBL数据库 (已知靶点补充)
1. 使用 ChEMBL API 获取化合物-靶点活性记录（Human）
2. 建议筛选 IC50/EC50/Ki/Kd 类型，保留已知靶点
3. 作为已知靶点补充，保留 Source=ChEMBL_known
4. 作为已知靶点补充来源记录在 `tcm_targets_ezhu.csv`

#### 选项E: BATMAN-TCM数据库
1. 访问 http://bionet.ncpsb.org.cn/batman-tcm/
2. 搜索莪术相关成分
3. 下载预测靶点

### 文件格式要求
当前主线使用 **HERB 导出靶点**，并可叠加 **ChEMBL 已知靶点**，最低要求是 **Target 列**：
```csv
Target,Source,Reference_ID,PubMed_ID
AKT1,HERB_reference_mining,HBREF002152,33866197
CASP3,HERB_statistical_inference,,
EGFR,ChEMBL_known,,
```

- **Target**: 靶点基因符号 (HGNC标准)
- **Source**: HERB_reference_mining / HERB_statistical_inference / ChEMBL_known

### 保存位置
`data/references/tcm_targets_ezhu.csv`

### 可选补充文件
`data/references/tcm_targets_ezhu_chembl.csv`（ChEMBL 扩展靶点明细）

---

## 1.2 莪术配体信息（ligands_smiles_ezhu.csv）

### 数据来源
- HERB 导出的活性成分列表（含 SMILES）
- 可选补充：TCMSP 成分表（用于 OB/DL 筛选）

### 文件格式要求
```csv
Ligand,PubChem_CID,SMILES,SDF_File
Curcuma_Compound_A,123456,SMILES_STRING,Curcuma_Compound_A.sdf
Curcuma_Compound_B,789012,SMILES_STRING,Curcuma_Compound_B.sdf
```

- **Ligand**: 活性成分英文名（建议用下划线替代空格）
- **PubChem_CID**: 可选
- **SMILES**: 必填，用于 Boltz-2 对接
- **SDF_File**: 可选，用于传统对接（缺省则按 `Ligand.sdf` 推断）

### 保存位置
`data/references/docking/ligands_smiles_ezhu.csv`

---

## 2. 铁死亡基因集

### 数据来源（多来源整合）

本项目整合了三个权威来源的铁死亡相关基因：

#### 来源1: FerrDb V2 核心基因 (15个)
- **数据库**: FerrDb V2 (http://www.zhounan.org/ferrdb/)
- **引用**: Zhou N, et al. FerrDb V2: update of the manually curated database of ferroptosis regulators and ferroptosis-disease associations. Nucleic Acids Res. 2023;51(D1):D571-D582.
- **类别**: Core_Ferroptosis

#### 来源2: KEGG铁死亡通路 (2个独特基因)
- **通路ID**: hsa04216 (Ferroptosis)
- **数据库**: KEGG Pathway Database (https://www.kegg.jp/)
- **引用**: Kanehisa M, et al. KEGG: Kyoto Encyclopedia of Genes and Genomes. Nucleic Acids Res. 2000;28(1):27-30.
- **类别**: KEGG_Ferroptosis

#### 来源3: MSigDB基因集 (91个)
- **数据库**: Molecular Signatures Database (https://www.gsea-msigdb.org/)
- **基因集**: 包含4个铁死亡相关基因集
- **引用**: Liberzon A, et al. The Molecular Signatures Database Hallmark Gene Set Collection. Cell Syst. 2015;1(6):417-425.
- **类别**: MSigDB_Ferroptosis

### 可用文件

#### ferroptosis_genes_CHB.csv (原始版本, 45个基因)
- 包含核心铁死亡调控基因
- 适用于初步筛选

#### ferroptosis_genes_expanded.csv（FerrDb V2 静态参考池，n=1290；推荐用于溯源/审计）
**当前复现包主线口径**：以 **FerrDb V2** 作为铁死亡参考基因池（便于审计“分母”），并在分析阶段再结合 HCC 队列 DEG 生成 HCC-context 列表（用于主线分析“分子”）。

- `data/references/ferroptosis_genes_expanded.csv`：FerrDb V2 整理后的静态参考池（不应被分析脚本覆盖）
- `results/ferroptosis_genes_hcc_context.csv`：HCC-context 铁死亡基因（FerrDb ∩ GSE14520-DEG，含 logFC/adjP/方向），由 `scripts_final/00c_prepare_hcc_ferro_genes.R` 生成

#### ferroptosis_genes_simple.csv (简化版本)
- 仅包含基因符号列
- 适用于简单交集分析

### 文件格式
```csv
Gene,FerrDb_Category,Source
GPX4,Suppressor,FerrDb_V2_early_preview_20231231
PTGS2,Marker,FerrDb_V2_early_preview_20231231
...
```

### 保存位置
- `data/references/ferroptosis_genes_expanded.csv` - FerrDb V2 静态参考池（推荐用于溯源）
- `results/ferroptosis_genes_hcc_context.csv` - HCC-context 衍生结果（主线分析使用）

### 生成脚本（HCC-context）
`scripts_final/00c_prepare_hcc_ferro_genes.R` - 使用 GSE14520 DEG 对 FerrDb 参考池进行过滤并输出到 `results/`

---

## 3. 免疫检查点基因 (immune_checkpoints.csv)

### 数据来源
文献来源，需标注引用:
- Pardoll DM. The blockade of immune checkpoints in cancer immunotherapy. Nat Rev Cancer. 2012.
- 或使用 ImmPort 数据库: https://www.immport.org/

### 文件格式要求
```csv
Gene,Type
PDCD1,Inhibitory
CD274,Inhibitory
CTLA4,Inhibitory
...
```

### 保存位置
`data/references/immune_checkpoints.csv`

---

## 4. 免疫细胞标志物 (immune_signatures.rds)

### 数据来源
- CIBERSORT LM22 signature matrix
- 或 xCell signatures
- 或文献定义的marker基因集

### 保存位置
`data/references/immune_signatures.rds`

---

## 数据验证清单

在运行分析前，请确认以下文件存在且格式正确:

- [ ] `data/references/tcm_targets_ezhu.csv` - 莪术靶点 (来自HERB)
- [ ] `data/references/tcm_targets_ezhu_chembl.csv` - ChEMBL 已知靶点扩展明细 (可选)
- [ ] `data/references/ferroptosis_genes_expanded.csv` - FerrDb V2 静态参考池
- [ ] `results/ferroptosis_genes_hcc_context.csv` - HCC-context 铁死亡基因（由 00c 生成）
- [ ] `data/references/immune_checkpoints.csv` - 免疫检查点基因
- [ ] `data/references/immune_signatures.rds` - 免疫细胞标志物
- [ ] `data/references/docking/ligands_smiles_ezhu.csv` - 莪术配体 (用于对接)
- [ ] `data/references/raw/herb_ingredient_2025_12_21.xlsx` - HERB 成分原始导出
- [ ] `data/references/raw/herb_reference_target_2025_12_21.xlsx` - HERB 靶点（参考挖掘）
- [ ] `data/references/raw/herb_target_2025_12_21.xlsx` - HERB 靶点（统计推断）
- [ ] `data/references/raw/tcmsp_curcumae_rhizoma_ingredients.md` - TCMSP 成分原始导出

---

## 引用要求

在论文中必须引用数据来源:

### HERB (莪术靶点)
> Fang S, et al. HERB: a high-throughput experiment- and reference-guided database of traditional Chinese medicine. Nucleic Acids Res. 2021;49(D1):D1197-D1206.

### TCMSP (可选补充)
> Ru J, et al. TCMSP: a database of systems pharmacology for drug discovery from herbal medicines. J Cheminform. 2014;6:13.

### ChEMBL (已知靶点补充)
> Gaulton A, et al. The ChEMBL database in 2017. Nucleic Acids Res. 2017;45(D1):D945-D954.

### FerrDb V2 (铁死亡核心基因)
> Zhou N, et al. FerrDb V2: update of the manually curated database of ferroptosis regulators and ferroptosis-disease associations. Nucleic Acids Res. 2023;51(D1):D571-D582.

### KEGG (铁死亡通路)
> Kanehisa M, et al. KEGG: Kyoto Encyclopedia of Genes and Genomes. Nucleic Acids Res. 2000;28(1):27-30.

### MSigDB (铁死亡基因集)
> Liberzon A, et al. The Molecular Signatures Database Hallmark Gene Set Collection. Cell Syst. 2015;1(6):417-425.

### GEO数据集
> Barrett T, et al. NCBI GEO: archive for functional genomics data sets. Nucleic Acids Res. 2013;41:D991-D995.

---

## 常见问题

### Q: TCMSP网站打不开怎么办?
A: 尝试使用VPN，或使用HERB数据库作为替代。

### Q: HERB靶点数量很少怎么办?
A: 可在 HERB 成分基础上使用 ChEMBL API 生成已知靶点，并在方法中明确标注为“已知来源”。建议保留原始 HERB 靶点作为高可信子集。

### Q: 下载的数据格式不对怎么办?
A: 使用Excel或R整理为上述标准格式。

### Q: 靶点基因符号不是HGNC标准怎么办?
A: 使用 biomaRt 或 org.Hs.eg.db 进行基因符号转换。
