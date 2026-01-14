#!/usr/bin/env Rscript

# run_complete_pipeline.R
# HCC 主线分析管道 - 从头到尾的可复现分析
# 说明：仅保留主线脚本，支持一键复现
# 
# 使用方法: Rscript scripts_final/run_complete_pipeline.R

# ============================================================================
# 禁用交互式图形设备
# ============================================================================
options(device = pdf)
graphics.off()
Sys.setenv(DISPLAY = "")

message("================================================================================")
message("HCC 莪术 (Ezhu)-铁死亡-免疫研究 - 完整分析管道")
message("================================================================================")
message("开始时间: ", Sys.time())

# 检查工作目录
if (!dir.exists("data/processed")) {
  stop("错误: 请在项目根目录运行此脚本")
}

# ============================================================================
# Stage 1: 多队列整合与验证
# ============================================================================
message("\n[Stage 1] 数据准备与基础分析...")
source("scripts_final/00_setup_env.R")
# 准备 HCC 特异 ferroptosis 基因集
source("scripts_final/00c_prepare_hcc_ferro_genes.R")
source("scripts_final/01b_download_GSE14520.R")
source("scripts_final/02b_multi_cohort_DEG.R")

# ============================================================================
# Stage 2: 临床关联分析 (核心更新部分)
# ============================================================================
message("\n[Stage 2] 预后模型构建 (Ezhu 9基因)...")
source("scripts_final/02c_prognostic_model_ezhu.R")
source("scripts_final/02d_nomogram_calibration_dca.R")
source("scripts_final/02e_external_validation.R")
source("scripts_final/02f_mvi_validation.R")

# ============================================================================
# Stage 3: 免疫微环境分析
# ============================================================================
message("\n[Stage 3] 免疫微环境分析...")
source("scripts_final/03a_immune_infiltration.R")
source("scripts_final/03b_immune_checkpoint.R")

# ============================================================================
# Stage 4: 莪术网络药理学与分子对接
# ============================================================================
message("\n[Stage 4] 分子对接与药物敏感性...")

# 药物敏感性 (V2)
source("scripts_final/05c_drug_sensitivity.R")

# 分子对接热图绘制（仅在有结果时运行）
dock_result <- "results/cbdock2_docking_results.csv"
if (file.exists(dock_result) && file.exists("scripts_final/06_molecular_docking_heatmap.R")) {
  source("scripts_final/06_molecular_docking_heatmap.R")
} else {
  message("[Stage 4] 跳过分子对接热图（缺少对接结果）")
}

# ============================================================================
# Stage 5: 最终出版级图表生成
# ============================================================================
message("\n[Stage 5] 生成最终出版级图表 (Figures 2-6)...")
source("scripts_final/07_final_publication_figures.R")

message("\n================================================================================")
message("分析管道全部完成。")
message("请查看 'plots/publication/' 目录获取最终图表。")
message("================================================================================")
