# ==============================================================================
# 03_fig2be_untargeted_mzmine.R
# 复现 ZimmermannLab/transwell_2026 Figure 2B-E (11 株菌非靶向 MZmine 转化产物挖掘):
#
# 上游源码对应:
#   upstream_transwell_2026/Fig2B_E_untar_analysis.Rmd
# 输入文件:
#   - input_folder/MZMINE_<strain>__iimn_gnps_quant.csv (11 个菌株非靶向矩阵)
#   - input_folder/NEW_pooling_scheme_info_with_IS.csv (药物池分配方案)
#   - input_folder/AB005_results_p02_10percent_241126.csv (由 Step 02 生成的靶向降解清单)
# 输出产物:
#   - data/processed/TP8_ttest_results.RDS (2.27 GB 向量化检验结果缓存，修复原著断链)
#   - data/processed/240607_mzmine_candidates.csv & repro_results/03_fig2be_mzmine_candidates.tsv
#   - repro_results/03_fig2be_untargeted_summary.png
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(matrixStats)
  library(patchwork)
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

cat(">>> [1/5] 载入 11 组 MZmine 定量矩阵与药物池分配方案...\n")
data_list <- sort(list.files(in_dir, pattern = "^MZMINE_.*_quant\\.csv$", full.names = TRUE))
if (length(data_list) == 0) stop("input_folder 中未找到 MZMINE_*_quant.csv 文件")

mzmine_peaks <- list()
mzmine_areas <- list()
for (fp in data_list) {
  df <- fread(fp)
  df$feature_ID <- paste0("FT", df$`row ID`)
  mzmine_peak <- as.data.frame(df[, 1:13])
  mzmine_peak$feature_ID <- df$feature_ID
  mzmine_area <- df %>% dplyr::select(contains("pool"))
  mzmine_area <- cbind(feature_ID = df$feature_ID, as.data.frame(mzmine_area))
  pool_cols <- grep("pool", colnames(df), value = TRUE)
  strain <- sapply(strsplit(as.character(pool_cols[1]), "_"), "[", 2)
  mzmine_peaks[[strain]] <- mzmine_peak
  mzmine_areas[[strain]] <- mzmine_area
}
cat("    已加载菌株矩阵:", paste(names(mzmine_peaks), collapse = ", "), "\n")

pools <- read_csv(file.path(in_dir, "NEW_pooling_scheme_info_with_IS.csv"), show_col_types = FALSE) %>%
  filter(!grepl("IS_", Product.Name)) %>%
  na.omit()

# 转置单株菌的宽表面积矩阵为 "行 = 样本, 列 = 特征" 并拆解 sample_name 元信息
prep_allpeaks_t <- function(complete_df) {
  m_t <- setDT(data.table::transpose(complete_df, keep.names = "rn"))
  col_nms <- as.character(m_t[1, ])
  m_t <- m_t[-1, ]
  colnames(m_t) <- c("sample_name", col_nms[-1])
  feat_cols <- 2:ncol(m_t)
  for (j in feat_cols) set(m_t, j = j, value = as.numeric(m_t[[j]]))
  m_df <- as.data.frame(m_t)
  m_df$strain <- sapply(strsplit(as.character(m_df$sample_name), "_"), "[", 2)
  m_df$pool <- sapply(strsplit(as.character(m_df$sample_name), "_"), "[", 3)
  m_df$TP <- sapply(strsplit(as.character(m_df$sample_name), "_"), "[", 4)
  m_df
}

allpeaks_t <- lapply(mzmine_areas, prep_allpeaks_t)

# 仅保留 TP8 (8 小时终末孵育)，将缺失/零峰面积置底为 1 后进行 log2 转换供 Welch t 检验
prep_tp8_log <- function(m_df) {
  tp8 <- m_df %>% filter(TP == "TP8")
  feat_names <- setdiff(colnames(tp8), c("sample_name", "strain", "pool", "TP"))
  mat <- as.matrix(tp8[, feat_names, drop = FALSE])
  mat[mat == 0] <- 1
  list(mat = log2(mat), pool = tp8$pool, features = feat_names)
}
log_tp8 <- lapply(allpeaks_t, prep_tp8_log)

# 构建每种药物对应的 4 个含药池索引 (每个化合物按正交设计被分入 26 池中的 4 个唯一池子)
pool_cols_bin <- colnames(pools)[apply(pools == 1, 2, any)]
to_include <- apply(pools, 1, function(x) colnames(pools)[which(x == 1)])
names(to_include) <- pools$Product.Name

cat(">>> [2/5] 执行向量化 Welch t 检验并校验与 stats::t.test 数值一致性...\n")
rds_cache <- file.path(proc_dir, "TP8_ttest_results.RDS")

# 矩阵向量化不等方差 Welch t 检验：针对“含某药的 4 个池 (in_idx) vs 其余 22 个不含该药池 (out_idx)”一次性对全部 11k-25k 特征求解
vec_welch_ttest <- function(mat, in_idx, out_idx) {
  df1 <- mat[in_idx, , drop = FALSE]
  df2 <- mat[out_idx, , drop = FALSE]
  n1 <- colSums(!is.na(df1))
  n2 <- colSums(!is.na(df2))
  m1 <- colMeans(df1, na.rm = TRUE)
  m2 <- colMeans(df2, na.rm = TRUE)
  v1 <- colVars(df1, na.rm = TRUE)
  v2 <- colVars(df2, na.rm = TRUE)
  se1 <- v1 / n1
  se2 <- v2 / n2
  se_tot <- sqrt(se1 + se2)
  t_stat <- (m1 - m2) / se_tot
  df_w <- (se1 + se2)^2 / ((se1^2 / pmax(n1 - 1, 1)) + (se2^2 / pmax(n2 - 1, 1)))
  pval <- 2 * pt(-abs(t_stat), df = df_w)
  bad <- is.na(se_tot) | se_tot == 0 | n1 < 2 | n2 < 2
  pval[bad] <- 1
  t_stat[bad] <- 0
  padj <- p.adjust(pval, method = "BH")
  fc <- m1 - m2
  cbind(pvalue = pval, statistic = t_stat, padj = padj, FoldChange = fc)
}

if (!file.exists(rds_cache)) {
  TP8_pools_results <- list()
  for (st in names(log_tp8)) {
    obj <- log_tp8[[st]]
    res_st <- list()
    for (d in names(to_include)) {
      in_idx <- obj$pool %in% to_include[[d]]
      out_idx <- !in_idx
      if (sum(in_idx) >= 2 && sum(out_idx) >= 2) {
        m <- vec_welch_ttest(obj$mat, in_idx, out_idx)
        rownames(m) <- obj$features
        res_st[[d]] <- m
      }
    }
    TP8_pools_results[[st]] <- res_st
    cat(sprintf("    完成菌株 %s (%d 种药物 × %d 特征)\n", st, length(res_st), ncol(obj$mat)))
  }
  saveRDS(TP8_pools_results, rds_cache)
} else {
  TP8_pools_results <- readRDS(rds_cache)
}

# 数值一致性验证 (随机抽取第1菌株第1药前5个非零方差特征对比 stats::t.test)
s1 <- names(log_tp8)[1]
d1 <- names(to_include)[1]
idx_in1 <- log_tp8[[s1]]$pool %in% to_include[[d1]]
max_diff <- max(vapply(1:5, function(k) {
  ref <- t.test(log_tp8[[s1]]$mat[idx_in1, k], log_tp8[[s1]]$mat[!idx_in1, k], paired = FALSE)$p.value
  abs(ref - TP8_pools_results[[s1]][[d1]][k, "pvalue"])
}, numeric(1)), na.rm = TRUE)
cat(sprintf("    向量化 Welch t 检验与 stats::t.test 最大 p 值误差: %.2e\n", max_diff))

cat(">>> [3/5] 筛选显著池特异性离子、deltaMZ / deltaRT 与去 nobug 伪阳性...\n")
select_ions <- function(res_list) {
  hits <- lapply(names(res_list), function(d) {
    df <- as.data.frame(res_list[[d]])
    df$value <- rownames(df)
    df <- df %>% filter(!is.na(padj), padj < 0.01, abs(FoldChange) > 2)
    if (nrow(df) > 0) df$drug <- d
    df
  })
  bind_rows(hits)
}

ions_raw <- lapply(TP8_pools_results, select_ions)
ions_info <- lapply(names(ions_raw), function(st) {
  df <- ions_raw[[st]]
  if (nrow(df) == 0) return(NULL)
  peaks <- mzmine_peaks[[st]][, c(2, 3, 9, 13, 14)]
  df %>%
    left_join(peaks, by = c("value" = "feature_ID")) %>%
    left_join(pools[, c("Product.Name", "M.SD", "RT.NR")], by = c("drug" = "Product.Name")) %>%
    mutate(strain = st)
})
my_df <- bind_rows(ions_info) %>%
  mutate(
    deltaRT = RT.NR - `row retention time`,
    deltaMZ = `row m/z` - M.SD,
    round_mz = round(`row m/z`, 2)
  )

# 4 级假阳性漏斗过滤：
#   Step A: |deltaMZ| > 1 Da 且 |deltaRT| > 0.3 min -> 排除母药自身同位素峰 (M+1) 及离子源加合物
#   Step B: ion_count == 1 -> 同一 m/z 仅特异性富集于某单一母药的 4 个池中
#   Step C: anti_join(nobug) -> 剔除无菌培养基中本身也会产生的非酶促自发降解产物
my_df_fil <- my_df %>% filter(abs(deltaMZ) > 1, abs(deltaRT) > 0.3)
drug_candidates <- my_df_fil %>%
  group_by(strain, round_mz) %>%
  mutate(ion_count = n_distinct(drug)) %>%
  ungroup() %>%
  filter(ion_count == 1)

nobug_mzs <- drug_candidates %>% filter(strain == "nobug") %>% distinct(round_mz)
candidates <- anti_join(drug_candidates, nobug_mzs, by = "round_mz") %>%
  filter(strain != "nobug")

cat(sprintf("    显著特征 %d -> deltaMZ/RT过滤 %d -> 单药+非nobug特异性候选: %d\n",
            nrow(my_df), nrow(my_df_fil), nrow(candidates)))

cat(">>> [4/5] 验证 TP0->TP8 时程积聚与药物池覆盖 (>=3/4 池)，并与 Fig2A 降解药物交集对齐...\n")
fc_rows <- vector("list", nrow(candidates))
for (i in seq_len(nrow(candidates))) {
  val <- candidates$value[i]
  st <- candidates$strain[i]
  drg <- candidates$drug[i]
  df_st <- allpeaks_t[[st]]
  in_pools <- to_include[[drg]]
  df_sub <- df_st[, c("pool", "TP", val), drop = FALSE]
  colnames(df_sub)[3] <- "area"
  df_sub$in_pool <- df_sub$pool %in% in_pools
  mean_tp0 <- mean(df_sub$area[df_sub$in_pool & df_sub$TP == "TP0"], na.rm = TRUE)
  tp8_in <- df_sub$area[df_sub$in_pool & df_sub$TP == "TP8"]
  tp8_out <- df_sub$area[!df_sub$in_pool & df_sub$TP == "TP8"]
  mean_tp8 <- mean(tp8_in, na.rm = TRUE)
  med_tp8 <- median(tp8_in, na.rm = TRUE)
  mean_non <- mean(tp8_out, na.rm = TRUE)
  med_non <- median(tp8_out, na.rm = TRUE)
  n_pools_det <- sum(tapply(tp8_in, df_sub$pool[df_sub$in_pool & df_sub$TP == "TP8"], function(x) any(x > 0, na.rm = TRUE)))
  fc_rows[[i]] <- data.frame(
    value = val, strain = st, drug = drg,
    meanTP0 = mean_tp0, meanTP8 = mean_tp8, mean_non = mean_non,
    pool_percentage = ifelse(mean_tp8 > 0, (mean_non / mean_tp8) * 100, Inf),
    pool_median = ifelse(med_tp8 > 0, (med_non / med_tp8) * 100, Inf),
    n_pools_det = n_pools_det,
    stringsAsFactors = FALSE
  )
}

# Step D: 验证该特征在 TP8 的非目标池均值/中位数仅占目标池 <= 10% (<= 100%) 且至少在 >= 3/4 个含药池检出
fc_df <- bind_rows(fc_rows) %>%
  filter(pool_percentage <= 10, pool_median <= 100, n_pools_det >= 3)

integration_df <- inner_join(fc_df, candidates, by = c("value", "strain", "drug")) %>%
  dplyr::select(value, strain, drug, padj, FoldChange, `row m/z`, `row retention time`, deltaMZ, deltaRT, M.SD, RT.NR, pool_percentage, n_pools_det)

# 修复菌株缩写错配：Fig2A 使用 Buniformis/Efaecalis，而 MZMINE 矩阵列名为 Buni/Ef2
degraded_path <- file.path(in_dir, "AB005_results_p02_10percent_241126.csv")
if (!file.exists(degraded_path)) stop("缺少 Fig2A 生成文件: ", degraded_path)
degraded_drugs <- read.csv(degraded_path) %>%
  mutate(
    strain_norm = recode(strain, "Buniformis" = "Buni", "Efaecalis" = "Ef2"),
    combo = paste0(Product.Name, "_", strain_norm)
  )

integration_df <- integration_df %>%
  mutate(
    combo = paste0(drug, "_", strain),
    category = ifelse(combo %in% degraded_drugs$combo, "Drug_Metabolite", "Physiological_Metabolite")
  )

write.table(integration_df, file.path(res_dir, "03_fig2be_mzmine_candidates.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.csv(integration_df, file.path(proc_dir, "240607_mzmine_candidates.csv"), row.names = FALSE)

cat(">>> [5/5] 绘制 Fig 2B-E 菌株转换产物计数与 deltaMZ/deltaRT 分布图...\n")
strain_pretty <- c(
  "Anaeslundi" = "A. naeslundii", "Aomnicolens" = "A. omnicolens", "Btheta" = "B. thetaiotaomicron",
  "Buni" = "B. uniformis", "Ccomes" = "C. comes", "Cramosum" = "C. ramosum",
  "Cscindens" = "C. scindens", "Dformi" = "D. formicigenerans", "Ef2" = "E. faecalis",
  "Ghaemolysans" = "G. haemolysans"
)

drug_mets <- integration_df %>%
  filter(category == "Drug_Metabolite") %>%
  mutate(strain_label = recode(strain, !!!strain_pretty))

p_count <- drug_mets %>%
  group_by(strain_label) %>%
  summarise(count = n_distinct(value), .groups = "drop") %>%
  ggplot(aes(x = reorder(strain_label, count), y = count)) +
  geom_col(fill = "#C2185B") +
  coord_flip() +
  labs(title = "Confirmed Bacterial Drug Metabolite Ions (Fig 2B)", x = NULL, y = "No. of Drug Metabolite Features") +
  theme_bw(base_size = 9) +
  theme(axis.text.y = element_text(face = "italic"))

p_dmz <- ggplot(drug_mets, aes(x = deltaMZ)) +
  geom_histogram(bins = 35, fill = "#5F9EA0", color = "white") +
  labs(title = "Mass Shift Distribution (delta m/z)", x = "delta m/z (Metabolite - Parent)", y = "Count") +
  theme_bw(base_size = 9)

p_drt <- ggplot(drug_mets, aes(x = category, y = deltaRT)) +
  geom_boxplot(width = 0.3, fill = "#DDA0DD", outlier.shape = NA) +
  geom_jitter(width = 0.08, alpha = 0.45, size = 1.0) +
  labs(title = "Retention Time Shift (delta RT)", x = NULL, y = "RT_parent - RT_metabolite (min)") +
  theme_bw(base_size = 9)

combined <- (p_count | (p_dmz / p_drt)) + plot_layout(widths = c(1.1, 1))
ggsave(file.path(res_dir, "03_fig2be_untargeted_summary.png"), combined, width = 10.5, height = 6, dpi = 300)

cat(sprintf("    最终入选池特异性候选离子: %d 条 (药物转化产物 Drug_Metabolite: %d | 生理/诱发变化: %d)\n",
            nrow(integration_df), sum(integration_df$category == "Drug_Metabolite"), sum(integration_df$category == "Physiological_Metabolite")))
cat("03_fig2be_untargeted_mzmine.R 运行成功。\n")
