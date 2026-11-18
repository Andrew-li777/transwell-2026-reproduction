# ==============================================================================
# 00_smoke_test.R
# ZimmermannLab/transwell_2026 本地环境与数据就绪冒烟测试 (< 10s)
#
# 上游映射: 对应 Fig1_TEER_viability_script.R 核心质控加载验证
# 输入文件: input_folder/Caco2_viability_data.csv, input_folder/Caco2_TEER.csv 等核心文件存在性校验
# 输出产物: repro_results/00_smoke_qc.png, repro_results/00_smoke_verification_summary.tsv
#
# 验证目标:
#   1. 核验核心 CRAN / Bioconductor 包在 R 4.6.0 下可正常加载 (跳过依赖 JDK 的 xlsx)
#   2. 检查 input_folder/ 下 Zenodo 22044327 核心输入表与 11 组 MZMINE 矩阵是否就绪
#   3. 基于 Caco2_viability_data.csv 与 Caco2_TEER.csv 生成冒烟质控预览图与定量摘要
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(readxl)
  library(ComplexHeatmap)
  library(patchwork)
  library(ggforce)
  library(fuzzyjoin)
})

# 解析 --file= 自动定位项目根目录
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
log_dir <- file.path(proj_root, "logs")
# 标准化工作区目录：所有产物统一写入 repro_results/ 与 data/processed/，不再创建上游旧目录 Plot/Plots/Env/output
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

cat(">>> [1/3] 核验 Zenodo 22044327 核心输入文件...\n")
required_core <- c(
  "Caco2_viability_data.csv",
  "Caco2_TEER.csv",
  "AB012I_TEER_24well.xlsx",
  "Fig2A_drugs_all_raw_with_metadata.rds",
  "NEW_pooling_scheme_info_with_IS.csv",
  "AB012_drugs_mets_table.csv",
  "Data_AB012I_transwell_drugs.csv"
)
present_flags <- file.exists(file.path(in_dir, required_core))
mzmine_files <- list.files(in_dir, pattern = "^MZMINE_.*_quant\\.csv$")
cat(sprintf("    核心基准表: %d/%d 就绪 | MZMINE 菌株定量矩阵: %d 个\n",
            sum(present_flags), length(required_core), length(mzmine_files)))
if (!all(present_flags)) {
  stop("缺少核心输入文件: ", paste(required_core[!present_flags], collapse = ", "))
}

cat(">>> [2/3] 读取 Caco-2 细胞活力与 TEER 冒烟校验数据集...\n")
df_viab <- read.csv(file.path(in_dir, "Caco2_viability_data.csv"))
# 遵照原著 Fig1 质控规则：剔除 C. ramosum (Cr) 在 T0 即已破损的异常孔 (TEER < 1000 Ohm)
df_teer <- read.csv(file.path(in_dir, "Caco2_TEER.csv")) %>%
  filter(!(condition == "Cr" & name == "T0" & TEER < 1000))

p1 <- ggplot(df_viab, aes(x = condition, y = rel_viab, fill = bacterial_time)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(group = bacterial_time), position = position_jitterdodge(jitter.width = 0.08), size = 0.9) +
  scale_fill_manual(values = c("16_hour" = "#FF7F7F", "24_hour" = "#E75480")) +
  labs(title = "Smoke Check: Caco-2 Relative Viability", x = NULL, y = "Relative Fluo 560/590", fill = "Incubation") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")

p2 <- ggplot(df_teer, aes(x = name, y = TEER, fill = name)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 0.7, alpha = 0.7) +
  facet_wrap(~condition, nrow = 2) +
  scale_fill_manual(values = c("T0" = "#FFC0CB", "T4_HBSS" = "#5F9EA0"), labels = c("T0" = "0 h", "T4_HBSS" = "4 h")) +
  labs(title = "Smoke Check: Caco-2 TEER Across Strains", x = "Exposure", y = "TEER (Ohm)", fill = "Time") +
  theme_bw(base_size = 10) +
  theme(legend.position = "top")

combined <- p1 / p2 + plot_annotation(title = "transwell_2026 Smoke Test: Barrier Integrity QC")
ggsave(file.path(res_dir, "00_smoke_qc.png"), combined, width = 10, height = 7, dpi = 200)

cat(">>> [3/3] 导出冒烟测试指标摘要...\n")
summary_df <- data.frame(
  Metric = c(
    "Core_Inputs_Present",
    "MZMINE_Matrices_Count",
    "Viability_Rows",
    "TEER_Rows_Filtered",
    "Median_Viability_Overall",
    "Median_TEER_T0",
    "Median_TEER_T4_HBSS"
  ),
  Value = c(
    sprintf("%d/%d", sum(present_flags), length(required_core)),
    as.character(length(mzmine_files)),
    as.character(nrow(df_viab)),
    as.character(nrow(df_teer)),
    sprintf("%.3f", median(df_viab$rel_viab, na.rm = TRUE)),
    sprintf("%.1f", median(df_teer$TEER[df_teer$name == "T0"], na.rm = TRUE)),
    sprintf("%.1f", median(df_teer$TEER[df_teer$name == "T4_HBSS"], na.rm = TRUE))
  )
)
write.table(summary_df, file.path(res_dir, "00_smoke_verification_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
print(summary_df)
cat("冒烟测试完成！图表已保存:", file.path(res_dir, "00_smoke_qc.png"), "\n")
