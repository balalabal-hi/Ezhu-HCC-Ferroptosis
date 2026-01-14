# 论文图表清单

**更新时间**: 2024-12-06  
**项目**: 莪术 (Curcuma phaeocaulis) 调控HCC铁死亡的网络药理学与机制研究

---

## 主图 (Main Figures)

| 图号 | 文件名 | 内容 | 状态 |
|------|--------|------|------|
| Figure 1 | figure1.jpeg | 研究流程示意图（生成提示词见 `docs/figure1-generation-prompt.md`） | 需更新 |
| Figure 2 | Figure2_DEG_analysis.pdf | DEG分析 (火山图+Venn图+Hub基因表达) | 已更新 |
| Figure 3 | Figure3_prognostic_model.pdf | 预后模型 (KM曲线+ROC+Risk Score+系数) | 已更新 |
| Figure 4 | Figure4_clinical_utility.pdf | 临床应用 (单/多因素Cox+DCA+C-index) | 已更新 |
| Figure 5 | Figure5_immune_analysis.pdf | 免疫分析 (相关性+分组差异+检查点+散点) | 已更新 |
| Figure 6 | Figure6_drug_docking.pdf | 药物敏感性 + 分子对接 (CB-Dock2 + Boltz-2) | 已更新 |

---

## 补充图 (Supplementary Figures)

| 图号 | 文件名 | 内容 | 状态 |
|------|--------|------|------|
| Figure S1 | FigureS_external_validation.pdf | 外部验证 (TCGA-LIHC + GSE76427 + GSE10143) | 已更新 |

---

## 分子对接结果详情

> 莪术 (Curcuma phaeocaulis) 的分子对接结果需在靶点与配体重建后重新生成，本节待更新。

---

## 数据文件位置

| 类型 | 路径 |
|------|------|
| 图表 | Ezhu_HCC_Project/plots/ |
| 结果 | Ezhu_HCC_Project/results/ |
| 脚本 | Ezhu_HCC_Project/scripts_final/ |
| 对接原始数据 | Ezhu_HCC_Project/data/references/docking/cbdock2_results/ |

---

## 方法学引用

### CB-Dock2 分子对接
- **平台**: CB-Dock2 (https://cadd.labshare.cn/cb-dock2/)
- **算法**: AutoDock Vina
- **引用**: Liu Y, et al. CB-Dock2: improved protein-ligand blind docking by integrating cavity detection, docking and homologous template fitting. Nucleic Acids Res. 2022.

### NVIDIA Boltz-2 (补充)
- **平台**: NVIDIA NIM (https://build.nvidia.com/mit/boltz2)
- **引用**: Wohlwend J, et al. Boltz-2: Towards Accurate and Efficient Biomolecular Structure Prediction. bioRxiv. 2024.
