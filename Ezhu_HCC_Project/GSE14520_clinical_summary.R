#!/usr/bin/env Rscript
# GSE14520_clinical_summary.R
# 输出一份临床信息摘要：Tissue, Disease state, 等字段分布
library(GEOquery)

clin <- readRDS('data/processed/GSE14520_clinical.rds')
expr <- readRDS('data/processed/GSE14520_expr.rds')

cat('--- GSE14520 临床信息摘要 ---\n')
cat('总样本数: ', nrow(clin), '\n')
cat('基因数: ', nrow(expr), '\n\n')

# 解析 characteristics（可能分布在多个列：characteristics_ch1, .1, .2 ...）
char_cols <- grep('^characteristics_ch1', colnames(clin), value = TRUE)
chars <- apply(clin[, char_cols, drop = FALSE], 1, function(row) paste(na.omit(as.character(row)), collapse = '\n'))
# 尝试解析出 Tissue, Disease state, Individual
extract_field <- function(txt, key) {
  m <- regmatches(txt, regexpr(paste0(key, ':\\s*[^\\n]+'), txt, ignore.case = TRUE))
  ifelse(length(m) > 0 & nchar(m) > 0, sub(paste0('.*', key, ':\\s*'), '', m, ignore.case = TRUE), NA)
}
tissue <- sapply(chars, extract_field, key = 'Tissue')
disease <- sapply(chars, extract_field, key = 'Disease state')
individual <- sapply(chars, extract_field, key = 'Individual')

table_tissue <- table(tissue, useNA='ifany')
table_disease <- table(disease, useNA='ifany')
cat('组织类型:\n')
print(table_tissue)
cat('\n疾病状态:\n')
print(table_disease)

# 按 Tissue/Disease 做简单分组
group <- paste(tissue, disease, sep='; ')
group_tidy <- gsub('NA;.*|.*; NA', 'Unknown', group, ignore.case=TRUE)
tbl_group <- table(group_tidy, useNA='ifany')
cat('\n分组合并结果:\n')
print(tbl_group)

# 检查是否有配对样本（相同的 Individual 既有 Tumor 又有 Non-Tumor）
ind_table <- table(individual, tissue)
col_needed <- intersect(colnames(ind_table), c('Liver Tumor Tissue', 'Liver Non-Tumor Tissue', 'Liver'))
paired_count <- 0
if (length(col_needed) >= 2) {
  paired <- ind_table[, col_needed, drop = FALSE]
  paired_pairs <- rowSums(paired > 0) >= 2
  paired_count <- sum(paired_pairs)
  cat('\n配对样本数 (同一患者多组织类型): ', paired_count, '\n')
  if (paired_count > 0) {
    cat('配对患者列表 (前10):\n')
    print(head(names(paired_pairs)[paired_pairs], 10))
  }
} else {
  cat('\n未检测到可用于配对的两类组织列（例如 Tumor/Non-Tumor）。\n')
}

# 输出简易 CSV 以备手动扩展临床信息
out_df <- data.frame(
  Sample = clin$geo_accession,
  Tissue = tissue,
  Disease = disease,
  Individual = individual,
  Group = group_tidy
)
write.csv(out_df, file='data/processed/GSE14520_summary.csv', row.names=FALSE, quote=FALSE)
cat('\n已保存 GSE14520 临床摘要到: data/processed/GSE14520_summary.csv\n')
cat('建议下一步: 若能找到原始论文中的生存时间表格，可按 Individual 匹配补充。\n')
cat('目前可用于：肿瘤 vs 非肿瘤配对分析；若加 GSE81252 引入分期，可在 Stage 2 做分层。\n')