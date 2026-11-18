# ==============================================================================
# 05_fig4_api_baso_eic_and_stats.R
# 第二阶段：复现 ZimmermannLab/transwell_2026 Figure 4 (DMSO 非靶向顶侧/基底侧统计与通透 EIC)
#
# 上游源码对应:
#   - Fig4_analysis.R + compare_tonnobug_mzmine.R (Apical / Basolateral 双侧统计检验与注释)
#   - Fig4CD_plots.R                              (Vanillic acid / 5-hydroxynorvaline 双室 EIC 曲线)
# 输入文件:
#   - input_folder/mzmine_dmso_api_baso__quant.csv (35,434 特征 × 80 DMSO 双室样本列)
#   - input_folder/FT7682_ion.csv & FT5225_ion.csv (预提 EIC 峰表)
#   - 本地缓存 MassBank CompDb (AnnotationHub AH116166)
# 输出产物:
#   - data/processed/AB007_api_ttest_candidates.csv & AB007_baso_ttest_candidates.csv
#   - repro_results/05_fig4_api_ttest_candidates.tsv & 05_fig4_baso_ttest_candidates.tsv
#   - repro_results/05_fig4_api_baso_volcano_and_counts.png & 05_fig4cd_vanillic_acid_and_5hnv_eic.png
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(matrixStats)
  library(ComplexHeatmap)
  library(AnnotationHub)
  library(CompoundDb)
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

cell_plot_theme <- theme(
  text = element_text(family = "sans", size = 9),
  axis.text = element_text(size = 7),
  axis.title = element_text(size = 8),
  legend.text = element_text(size = 7),
  legend.title = element_text(size = 7),
  plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
  panel.background = element_blank(),
  panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  plot.margin = margin(5, 5, 5, 5)
)

cat(">>> [1/4] 载入 mzmine_dmso_api_baso__quant.csv (35,434 特征 × 80 样本)...\n")
quant_raw <- fread(file.path(in_dir, "mzmine_dmso_api_baso__quant.csv"))
quant_df <- as.data.frame(quant_raw)
quant_df$feature_ID <- paste0("FT", quant_df$`row ID`)

hydrogen <- 1.007825
peak_info <- quant_df[, c(1:13)]
peak_info$feature_ID <- quant_df$feature_ID
# 当 MZmine 未推断出加合离子形式 (best ion 为空) 时，参照原著公式以质子化精确质量反推近似中性分子质量
peak_info$`neutral M mass` <- ifelse(
  is.na(peak_info$`best ion`) | peak_info$`best ion` == "",
  peak_info$`row m/z` + hydrogen,
  peak_info$`neutral M mass`
)

# 向量化双室对比函数 (替代上游 compare_tonnobug_mzmine.R 串行 t.test，将 ~31.8 万次单边检验加速至 <3 秒)
run_side_vs_nobug <- function(quant_tbl, side_tag) {
  pos_cols <- grep("pos", colnames(quant_tbl), value = TRUE)
  side_cols <- grep(paste0("_", side_tag, "_"), pos_cols, value = TRUE)
  mat <- as.matrix(quant_tbl[, side_cols, drop = FALSE])
  mat[is.na(mat) | mat == 0] <- 1
  log_mat <- t(log2(mat)) # 行=样本，列=特征
  colnames(log_mat) <- quant_tbl$feature_ID

  strains <- sapply(strsplit(side_cols, "_"), "[", 2)
  nobug_idx <- strains == "nobug"
  if (sum(nobug_idx) < 2) stop("nobug 样本数不足: ", side_tag)

  test_strains <- setdiff(unique(strains), "nobug")
  out_list <- lapply(test_strains, function(st) {
    st_idx <- strains == st
    df1 <- log_mat[st_idx, , drop = FALSE]
    df2 <- log_mat[nobug_idx, , drop = FALSE]
    n1 <- nrow(df1)
    n2 <- nrow(df2)
    m1 <- colMeans(df1, na.rm = TRUE)
    m2 <- colMeans(df2, na.rm = TRUE)
    v1 <- colVars(df1, na.rm = TRUE)
    v2 <- colVars(df2, na.rm = TRUE)
    se1 <- v1 / n1
    se2 <- v2 / n2
    se_tot <- sqrt(se1 + se2)
    t_stat <- (m1 - m2) / se_tot
    df_w <- (se1 + se2)^2 / ((se1^2 / max(n1 - 1, 1)) + (se2^2 / max(n2 - 1, 1)))
    pval <- 2 * pt(-abs(t_stat), df = df_w)
    bad <- is.na(se_tot) | se_tot == 0
    pval[bad] <- 1
    t_stat[bad] <- 0
    padj <- p.adjust(pval, method = "BH")
    data.frame(
      strain = st,
      ion = colnames(log_mat),
      pval = pval,
      statistic = t_stat,
      padj = padj,
      FC = m1 - m2,
      stringsAsFactors = FALSE
    )
  })
  bind_rows(out_list)
}

cat(">>> [2/4] 对 Apical 与 Basolateral 两侧执行菌株 vs nobug 向量化 Welch t 检验...\n")
api_res <- run_side_vs_nobug(quant_df, "api")
baso_res <- run_side_vs_nobug(quant_df, "baso")

# 遵照原著双侧差异阈值：Apical（菌群原生环境信号高）要求 log2FC > 5；Basolateral（跨膜后浓度稀释）放宽至 log2FC > 3
api_sig <- api_res %>% filter(FC > 5, padj < 0.05) %>%
  left_join(peak_info %>% dplyr::select(feature_ID, `row ID`, `row m/z`, `row retention time`, `neutral M mass`),
            by = c("ion" = "feature_ID"))
baso_sig <- baso_res %>% filter(FC > 3, padj < 0.05) %>%
  left_join(peak_info %>% dplyr::select(feature_ID, `row ID`, `row m/z`, `row retention time`, `neutral M mass`),
            by = c("ion" = "feature_ID"))

cat(sprintf("    Apical 显著菌源特征 (FC>5, FDR<0.05): %d 条 (%d 唯一特征)\n", nrow(api_sig), n_distinct(api_sig$ion)))
cat(sprintf("    Basolateral 显著菌源特征 (FC>3, FDR<0.05): %d 条 (%d 唯一特征)\n", nrow(baso_sig), n_distinct(baso_sig$ion)))

cat(">>> [3/4] 调用本地缓存 MassBank (AH116166) 执行候选特征 MS1 中性质量匹配...\n")
# 本地合规注释函数：直接加载 BiocFileCache 中的 MassBank AH116166 CompDb 对象，按 ±0.002 Da 或 ±20 ppm 双窗匹配中性质量
annotate_with_local_massbank <- function(sig_df) {
  ah_ok <- tryCatch({
    ah <- AnnotationHub(localHub = TRUE)
    mb_db <- ah[["AH116166"]]
    cmp <- compounds(mb_db, columns = c("compound_id", "name", "exactmass", "inchikey"))
    cmp <- cmp[!is.na(cmp$exactmass) & cmp$exactmass > 0, ]
    list(ok = TRUE, cmp = cmp)
  }, error = function(e) list(ok = FALSE, msg = conditionMessage(e)))

  if (!ah_ok$ok) {
    message("    [提示] 本地 AnnotationHub AH116166 不可用: ", ah_ok$msg)
    sig_df$massbank_ms1_match <- NA_character_
    return(sig_df)
  }

  cmp <- ah_ok$cmp
  uniq_ions <- sig_df %>% dplyr::select(ion, `neutral M mass`) %>% distinct()
  matches <- vapply(seq_len(nrow(uniq_ions)), function(j) {
    cm <- uniq_ions$`neutral M mass`[j]
    if (is.na(cm)) return(NA_character_)
    diff_m <- abs(cm - cmp$exactmass)
    hit <- which(diff_m <= 0.002 | diff_m <= (20e-6 * cm))
    if (length(hit) == 0) return(NA_character_)
    paste(head(unique(cmp$name[hit]), 3), collapse = "; ")
  }, character(1))
  uniq_ions$massbank_ms1_match <- matches
  left_join(sig_df, uniq_ions[, c("ion", "massbank_ms1_match")], by = "ion")
}

api_sig <- annotate_with_local_massbank(api_sig)
baso_sig <- annotate_with_local_massbank(baso_sig)

write.csv(api_sig, file.path(proc_dir, "AB007_api_ttest_candidates.csv"), row.names = FALSE)
write.csv(baso_sig, file.path(proc_dir, "AB007_baso_ttest_candidates.csv"), row.names = FALSE)
write.table(api_sig, file.path(res_dir, "05_fig4_api_ttest_candidates.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(baso_sig, file.path(res_dir, "05_fig4_baso_ttest_candidates.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

cat(sprintf("    其中命中本地 MassBank MS1 候选注释: Apical %d 唯一特征 | Basolateral %d 唯一特征\n",
            n_distinct(api_sig$ion[!is.na(api_sig$massbank_ms1_match)]),
            n_distinct(baso_sig$ion[!is.na(baso_sig$massbank_ms1_match)])))

# 绘制 Fig 4 火山图与菌株命中对比图
p_vol_api <- ggplot(api_res, aes(x = FC, y = -log10(pmax(padj, 1e-300)))) +
  geom_point(aes(color = FC > 5 & padj < 0.05), size = 0.6, alpha = 0.6) +
  scale_color_manual(values = c("FALSE" = "grey75", "TRUE" = "#4A90E2"), labels = c("FALSE" = "NS", "TRUE" = "Enriched")) +
  geom_vline(xintercept = 5, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  facet_wrap(~strain, nrow = 3) +
  labs(title = "Apical Metabolites vs Sterile Control (FC > 5, FDR < 0.05)", x = "log2 Fold Change", y = "-log10(BH FDR)", color = NULL) +
  theme_bw(base_size = 8) +
  theme(legend.position = "top")

both_summary <- bind_rows(
  api_sig %>% group_by(strain) %>% summarise(Features = n_distinct(ion), Side = "Apical (FC>5)"),
  baso_sig %>% group_by(strain) %>% summarise(Features = n_distinct(ion), Side = "Basolateral (FC>3)")
)
p_bar <- ggplot(both_summary, aes(x = reorder(strain, Features), y = Features, fill = Side)) +
  geom_col(position = "dodge") +
  coord_flip() +
  scale_fill_manual(values = c("Apical (FC>5)" = "#4A90E2", "Basolateral (FC>3)" = "#E75480")) +
  labs(title = "Significant Bacterial Features by Compartment", x = NULL, y = "Unique Features", fill = NULL) +
  theme_bw(base_size = 9) +
  theme(legend.position = "top")

ggsave(file.path(res_dir, "05_fig4_api_baso_volcano_and_counts.png"), p_vol_api | p_bar + plot_layout(widths = c(1.6, 1)), width = 12, height = 6.5, dpi = 300)

cat(">>> [4/4] 复现 Figure 4C-D Vanillic acid (FT7682) 与 5-hydroxynorvaline (FT5225) 双室 EIC 时程图...\n")
parse_eic_ion_csv <- function(csv_name) {
  df <- read.csv(file.path(in_dir, csv_name))
  df$side <- sapply(strsplit(as.character(df$sample_name), "_"), "[", 7)
  df$time_raw <- sapply(strsplit(as.character(df$sample_name), "_"), "[", 4)
  df$pool <- sapply(strsplit(as.character(df$sample_name), "_"), "[", 6)
  df$strain <- sapply(strsplit(as.character(df$sample_name), "_"), "[", 2)
  df$bac_time <- sapply(strsplit(as.character(df$sample_name), "_"), "[", 3)
  df$time <- suppressWarnings(as.numeric(gsub("T", "", df$time_raw)))
  df %>%
    filter(!is.na(time), !is.na(intb)) %>%
    group_by(sample_name, ion, column) %>%
    slice_max(sn, n = 1, with_ties = FALSE) %>%
    ungroup()
}

ft7682 <- parse_eic_ion_csv("FT7682_ion.csv")
ft5225 <- parse_eic_ion_csv("FT5225_ion.csv")

strain_labs_4 <- c("Cramosum" = "C. ramosum", "Dformi" = "D. formicigenerans", "Ghaemo" = "G. haemolysans")

p_van <- ft7682 %>%
  filter(strain %in% names(strain_labs_4), bac_time == "TP8") %>%
  ggplot(aes(x = time, y = intb, color = side)) +
  geom_smooth(se = FALSE, method = "loess", formula = y ~ x, linewidth = 0.9) +
  geom_point(size = 1.0) +
  facet_grid(side ~ strain, scales = "free_y", labeller = labeller(strain = strain_labs_4, side = c("api" = "Apical", "baso" = "Basolateral"))) +
  scale_color_manual(values = c("api" = "#4A90E2", "baso" = "#E75480"), labels = c("api" = "Apical", "baso" = "Basolateral")) +
  scale_x_continuous(breaks = c(0, 120, 240)) +
  scale_y_continuous(labels = scales::scientific, n.breaks = 3) +
  labs(title = "Fig 4C: Vanillic acid (FT7682) Transwell Permeation", x = "Time (minutes)", y = "EIC Peak Area", color = NULL) +
  theme_minimal() +
  cell_plot_theme +
  theme(strip.text = element_text(face = "italic"), legend.position = "top")

p_5hnv <- ft5225 %>%
  filter(strain %in% c("Cramosum", "Dformi", "Ghaemo", "Buni", "Ccomes"), bac_time == "TP8") %>%
  ggplot(aes(x = time, y = intb, color = side)) +
  geom_smooth(se = FALSE, method = "loess", formula = y ~ x, linewidth = 0.9) +
  geom_point(size = 0.9) +
  facet_grid(side ~ strain, scales = "free_y", labeller = labeller(side = c("api" = "Apical", "baso" = "Basolateral"))) +
  scale_color_manual(values = c("api" = "#4A90E2", "baso" = "#E75480"), labels = c("api" = "Apical", "baso" = "Basolateral")) +
  scale_x_continuous(breaks = c(0, 120, 240)) +
  scale_y_continuous(labels = scales::scientific, n.breaks = 3) +
  labs(title = "Fig 4D: 5-Hydroxynorvaline (FT5225) Transwell Permeation", x = "Time (minutes)", y = "EIC Peak Area", color = NULL) +
  theme_minimal() +
  cell_plot_theme +
  theme(strip.text = element_text(face = "italic"), legend.position = "top")

ggsave(file.path(res_dir, "05_fig4cd_vanillic_acid_and_5hnv_eic.png"), p_van / p_5hnv, width = 10.5, height = 7.5, dpi = 300)

cat("05_fig4_api_baso_eic_and_stats.R 运行成功。\n")
