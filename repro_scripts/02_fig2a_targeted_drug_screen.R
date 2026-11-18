# ==============================================================================
# 02_fig2a_targeted_drug_screen.R
# 复现 ZimmermannLab/transwell_2026 Figure 2A (AB005 靶向药物菌群代谢筛选):
#
# 上游源码对应:
#   upstream_transwell_2026/Fig2A_tardrug_analysis_pipeline.R
# 输入文件:
#   - input_folder/Fig2A_drugs_all_raw_with_metadata.rds (载入 my_df: 14,008 × 18)
#   - input_folder/NEW_pooling_scheme_info_with_IS.csv (481 药 + 6 内标 × pool1-26)
# 输出产物:
#   - input_folder/AB005_results_p02_10percent_241126.csv (同步供 Step 03 识别药物转化产物)
#   - data/processed/02_fig2a_norm_df.rds (内标校正后的完整靶向峰面积表)
#   - repro_results/02_fig2a_degraded_drugs_results.tsv & 02_fig2a_strain_hit_counts.tsv
#   - repro_results/02_fig2a_drug_degradation_heatmap.png & 02_fig2a_example_timecourses.png
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(RColorBrewer)
  library(ComplexHeatmap)
  library(grid)
})

get_project_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE))
    return(normalizePath(file.path(script_dir, ".."), mustWork = TRUE))
  }
  getwd()
}

proj_root <- get_project_root()
in_dir <- file.path(proj_root, "input_folder")
res_dir <- file.path(proj_root, "repro_results")
proc_dir <- file.path(proj_root, "data", "processed")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)

cat(">>> [1/5] 载入 Fig2A 原始质谱积分表并按最小 deltaRT 挑选特征峰...\n")
load(file.path(in_dir, "Fig2A_drugs_all_raw_with_metadata.rds")) # 载入对象 my_df
pooling_scheme <- read.csv(file.path(in_dir, "NEW_pooling_scheme_info_with_IS.csv"))
pooling_scheme <- na.omit(pooling_scheme)

tmp <- left_join(my_df, pooling_scheme[, c(1, 30:35)], by = c("drug" = "drug"))
# 理论保留时间 RT.NR 原单位为分钟 (min)，转为秒 (s) 与原始积分 rt 对齐
tmp$RT.NR <- tmp$RT.NR * 60

# 若同一采样管中积分到多个候选异构/杂峰，仅保留距离理论保留时间最近的一个主峰
test <- tmp %>%
  mutate(deltaRT = abs(rt - RT.NR)) %>%
  group_by(sample_name, drug) %>%
  filter(deltaRT == min(deltaRT, na.rm = TRUE)) %>%
  slice(1) %>%
  ungroup()

cat(sprintf("    匹配标准库药物数: %d | 过滤后记录数: %d\n", n_distinct(test$Product.Name), nrow(test)))

cat(">>> [2/5] 执行内标 (Internal Standards) 中位数比例校正...\n")
my_prep_df <- as.data.table(test)
my_prep_df[grep("AB005", sample_name, fixed = TRUE), Status := "SAMPLE"]
my_prep_df[grep("IS_", drug, fixed = TRUE), Status := "INT_STD"]

# 遵照原著质控逻辑剔除 IS_sulfamethoxazole：部分肠道菌本身会还原代谢磺胺类内标，不适合做归一化参考
IS_drugs <- grep("IS_", unique(my_prep_df$drug), fixed = TRUE, value = TRUE)
IS_only_df <- my_prep_df[drug %in% IS_drugs & drug != "IS_sulfamethoxazole"] %>%
  as.data.frame()

# 计算每个样品中剩余内标的相对中位数偏差均值 intb_FC_FC_ISmedian，用以校正进样体积/信号漂移
median_IS <- IS_only_df %>%
  group_by(interaction(strain, drug)) %>%
  mutate(intb_Median = median(intb, na.rm = TRUE), intb_FC = intb / median(intb, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(interaction(strain, pool)) %>%
  mutate(intb_FC_FC_ISmedian = mean(intb_FC, na.rm = TRUE)) %>%
  ungroup() %>%
  dplyr::select(sample_name, intb_FC_FC_ISmedian) %>%
  distinct()

my_norm_df <- left_join(as.data.frame(my_prep_df), median_IS, by = "sample_name")
my_norm_df$correctedArea <- my_norm_df$intb / my_norm_df$intb_FC_FC_ISmedian

cat(">>> [3/5] 计算 TP8 与 TP2 细菌药物降解 log2FC 与 BH-FDR 显著性...\n")
# Efaecalis 在 AB005 靶向批次中由独立批次补测，遵照原著剔除后仅统计至少有 >=3 个平行重复 (N >= 3) 的组合
my_cols <- c("sample_name", "Product.Name", "rt", "intb", "pool", "strain", "TPs", "well_location", "correctedArea")
df <- my_norm_df[, my_cols] %>%
  filter(strain != "Efaecalis") %>%
  na.omit()

setDT(df)
to_keep <- df[, .N, by = .(strain, Product.Name, TPs)][N >= 3, ]
df <- df[na.omit(df[to_keep, on = .(strain, Product.Name, TPs), which = TRUE])]
df <- as.data.frame(df) %>% filter(!grepl("IS_", Product.Name))

# 单药单时间点降解判定函数：单尾 t 检验 (alternative = "less"，检验孵育后浓度是否显著低于 TP0)，
# 并按药物内部对所有菌株做 Benjamini-Hochberg FDR 校正，阈值 padj_strain < 0.2
calc_tp_degradation <- function(df_in, tp_target) {
  df1 <- df_in %>% filter(TPs == "TP0") %>% dplyr::select(sample_name, Product.Name, pool, strain, correctedArea)
  df2 <- df_in %>% filter(TPs == tp_target) %>% dplyr::select(sample_name, Product.Name, pool, strain, correctedArea)
  names(df1)[names(df1) == "correctedArea"] <- "t0"
  names(df2)[names(df2) == "correctedArea"] <- "t_end"

  df3 <- full_join(df1, df2, by = c("Product.Name", "pool", "strain")) %>%
    filter(!is.na(sample_name.x) & !is.na(sample_name.y))

  df_fc <- df3 %>%
    group_by(Product.Name, strain) %>%
    mutate(
      FC = log2(mean(t_end, na.rm = TRUE) / mean(t0, na.rm = TRUE)),
      FC_pval = tryCatch(t.test(t_end, t0, alternative = "less")$p.value, error = function(cond) 1)
    ) %>%
    ungroup()

  df_padj <- df_fc %>%
    dplyr::select(Product.Name, strain, FC, FC_pval) %>%
    distinct() %>%
    group_by(Product.Name) %>%
    mutate(padj_strain = p.adjust(FC_pval, method = "BH")) %>%
    ungroup()

  res <- df_padj %>%
    filter(FC < 0 & padj_strain < 0.2) %>%
    mutate(TP = tp_target)
  res
}

result_tp8 <- calc_tp_degradation(df, "TP8")
result_tp2 <- calc_tp_degradation(df, "TP2")
results <- bind_rows(result_tp8, result_tp2) %>%
  mutate(percentage = 100 - (2^FC) * 100) %>%
  filter(percentage > 10)

results_nobug <- results %>%
  filter(strain == "nobug") %>%
  dplyr::select(Product.Name, control = percentage, TP)

fil_results <- results %>%
  left_join(results_nobug, by = c("Product.Name", "TP")) %>%
  mutate(PassesFilter = ifelse(strain != "nobug" & (!is.na(control) & percentage > control), "yes", "no"))

# 仅保留细菌组本身通过滤膜（或无显著自发降解对照）的命中
results_final <- fil_results %>%
  filter(strain != "nobug") %>%
  filter(is.na(control) | PassesFilter == "yes") %>%
  dplyr::select(Product.Name, strain, FC, FC_pval, padj_strain, TP, percentage)

# 兼容 Fig2B-E 预期输入位置
write.csv(results_final, file.path(in_dir, "AB005_results_p02_10percent_241126.csv"), row.names = FALSE)
write.table(results_final, file.path(res_dir, "02_fig2a_degraded_drugs_results.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
saveRDS(df, file.path(proc_dir, "02_fig2a_norm_df.rds"))

cat(sprintf("    检出受菌群代谢降解药物-菌株对: %d 条 (涉及独立药物 %d 种)\n",
            nrow(results_final), n_distinct(results_final$Product.Name)))

cat(">>> [4/5] 绘制 Fig 2A 菌群-药物降解率热图与代表药物时程曲线...\n")
results_mat <- results_final %>%
  dplyr::select(Product.Name, strain, percentage, TP) %>%
  group_by(Product.Name, TP, strain) %>%
  summarise(percentage = max(percentage), .groups = "drop") %>%
  pivot_wider(names_from = strain, values_from = percentage, values_fill = 0)

strain_order_cols <- setdiff(colnames(results_mat), c("Product.Name", "TP"))
non_zero_counts <- colSums(results_mat[, strain_order_cols, drop = FALSE] != 0)
strain_order_cols <- names(sort(non_zero_counts, decreasing = TRUE))
row_sums <- rowSums(results_mat[, strain_order_cols, drop = FALSE])
results_mat <- results_mat[order(row_sums, decreasing = TRUE), ]

tps_vec <- results_mat$TP
mat_num <- as.matrix(results_mat[, strain_order_cols, drop = FALSE])
rownames(mat_num) <- results_mat$Product.Name

strain_pretty <- c(
  "Cramosum" = "C. ramosum", "Ccomes" = "C. comes", "Buni" = "B. uniformis",
  "Buniformis" = "B. uniformis",
  "Ghaemolysans" = "G. haemolysans", "Btheta" = "B. thetaiotaomicron",
  "Ef2" = "E. faecalis", "Anaeslundi" = "A. naeslundii",
  "Aomnicolens" = "A. omnicolens", "Dformi" = "D. formicigenerans", "Cscindens" = "C. scindens"
)
colnames(mat_num) <- ifelse(colnames(mat_num) %in% names(strain_pretty), strain_pretty[colnames(mat_num)], colnames(mat_num))

ht <- Heatmap(
  mat_num,
  name = "% Degradation",
  col = colorRampPalette(c("#F7F7F7", "#F06292", "#1A1A1A"))(100),
  right_annotation = rowAnnotation(TP = tps_vec, col = list(TP = c("TP8" = "#C2185B", "TP2" = "#F06292"))),
  row_title = sprintf("Degraded Drugs (%d drug-TP entries)", nrow(mat_num)),
  column_title = "Gut Bacterial Strains (Fig 2A)",
  show_row_names = nrow(mat_num) <= 60,
  row_names_gp = gpar(fontsize = 5),
  column_names_gp = gpar(fontsize = 9, fontface = "italic"),
  column_names_rot = 45,
  cluster_columns = FALSE,
  cluster_rows = FALSE
)

png(file.path(res_dir, "02_fig2a_drug_degradation_heatmap.png"), width = 8.5, height = 7.5, units = "in", res = 300)
draw(ht)
dev.off()

# 代表性药物降解时程验证 (Azilsartan, Roxatidine Acetate HCl, Carfilzomib)
example_drugs <- intersect(c("Azilsartan", "Roxatidine Acetate HCl", "Carfilzomib (PR-171)"), unique(df$Product.Name))
if (length(example_drugs) > 0) {
  p_ex <- df %>%
    filter(Product.Name %in% example_drugs, strain %in% c("Btheta", "Buni", "Ccomes", "Cramosum", "nobug")) %>%
    mutate(time_hr = as.numeric(sub("TP", "", TPs))) %>%
    ggplot(aes(x = time_hr, y = correctedArea, color = strain)) +
    geom_point(size = 1.1, alpha = 0.8) +
    stat_summary(fun = mean, geom = "line", linewidth = 0.8) +
    facet_grid(Product.Name ~ strain, scales = "free_y") +
    scale_y_continuous(labels = scales::scientific) +
    labs(title = "Representative Bacterial Drug Degradation Timecourses", x = "Incubation Time (h)", y = "IS-Normalized Peak Area") +
    theme_bw(base_size = 9) +
    theme(legend.position = "none")
  ggsave(file.path(res_dir, "02_fig2a_example_timecourses.png"), p_ex, width = 10, height = 5.5, dpi = 300)
}

cat(">>> [5/5] 输出各菌株降解药物统计:\n")
strain_counts <- results_final %>%
  group_by(strain, TP) %>%
  summarise(Drugs = n_distinct(Product.Name), .groups = "drop") %>%
  arrange(desc(Drugs))
print(as.data.frame(strain_counts))
write.table(strain_counts, file.path(res_dir, "02_fig2a_strain_hit_counts.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
cat("02_fig2a_targeted_drug_screen.R 运行成功。\n")
