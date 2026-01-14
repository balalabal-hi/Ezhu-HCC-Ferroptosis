# 项目说明文档 (README)

## 项目简介
本项目基于公开转录组数据，构建 **莪术 (Curcuma phaeocaulis, Ezhu)** 介导的肝细胞癌 (HCC) 铁死亡调控网络，并开展预后模型、免疫微环境、药敏与分子对接验证的整合分析。

## 目录结构
- `data/`: 数据
  - `raw/`: 原始下载数据 (GEO/TCGA)
  - `processed/`: 清洗后的中间数据 (RDS/CSV)
  - `references/`: 参考基因集与中药靶点
- `scripts_final/`: 论文主线脚本 (可复跑)
- `results/`: 结果表格
- `plots/`: 图版输出
- `manuscript/`: 手稿生成脚本与文档

## 主线运行方式

在 `Ezhu_HCC_Project/` 目录执行：

```bash
# 1. 环境准备 (安装所需 R 包)
Rscript scripts_final/00_setup_env.R

# 2. 运行完整分析管道
Rscript scripts_final/run_complete_pipeline.R
```

如只重生成最终图板：

```bash
Rscript scripts_final/07_final_publication_figures.R
```

## 可复现性说明 (Reproducibility)
本项目提供了全量复现包 `Ezhu_HCC_Reproducibility_Package.zip`，包含：
- 清洗后的表达矩阵与临床数据 (RDS 格式)
- 铁死亡与药靶参考基因集
- 完整的分析与绘图脚本
- 分子对接模型与 Boltz-2 AI 原始响应结果

## 引用建议 (Citation)
如果您在研究中使用了本项目的代码或数据，请引用：
> *Ferroptosis-based prognostic signature in hepatocellular carcinoma with Curcuma phaeocaulis target-guided mechanistic insights. (Submission v11)*

## 许可证
本项目采用 [MIT License](LICENSE) 授权。


## 关键手动数据

- **莪术靶点**: `data/references/tcm_targets_ezhu.csv`
  - 来源: HERB（高可信）+ ChEMBL（已知靶点补充）
  - 格式: `Target,Source`
- **莪术配体 (对接)**: `data/references/docking/ligands_smiles_ezhu.csv`
  - 格式: `Ligand,PubChem_CID,SMILES,SDF_File`
- **铁死亡参考基因池（FerrDb V2）**: `data/references/ferroptosis_genes_expanded.csv`
  - 用途：作为“分母/参考池”，便于审计与溯源（不应被分析脚本覆盖）
- **HCC-context 铁死亡基因（由 DEG 过滤后的衍生结果）**: `results/ferroptosis_genes_hcc_context.csv`
  - 生成脚本：`scripts_final/00c_prepare_hcc_ferro_genes.R`

详细说明见：`data/references/README_DATA_SOURCES.md` 与 `DATA_MANIFEST.md`。

## 外部验证
默认外部验证为 TCGA-LIHC、GSE76427 与 GSE10143，由 `scripts_final/02e_external_validation.R` 统一接入与更新。

## 说明
- `scripts/` 与 CHB 相关内容为历史遗留，不属于当前论文主线。
- 请以 `scripts_final/` 与 `results/` 为主线权威来源。
