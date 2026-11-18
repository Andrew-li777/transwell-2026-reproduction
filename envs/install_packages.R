# ==============================================================================
# ZimmermannLab/transwell_2026 依赖检查与安装脚本 (Windows + R 4.6.0)
#
# 用法:
#   Rscript reproduction/metabolomics/transwell_2026/envs/install_packages.R --dry-run
#   Rscript reproduction/metabolomics/transwell_2026/envs/install_packages.R
#
# 说明:
#   - 显式跳过依赖 rJava / 本地 JDK 的 CRAN 包 `xlsx`；所有 Excel 读写统一使用 `readxl` / `openxlsx`。
#   - 补齐复现重构阶段新增的高性能向量化与本地注释依赖: matrixStats, scales, CompoundDb, AnnotationHub, Biobase。
# ==============================================================================

# 解析命令行参数，--dry-run 仅打印包状态而不触发实际安装
args <- commandArgs(trailingOnly = TRUE)
dry_run <- "--dry-run" %in% args

# 默认使用清华 CRAN 镜像提升国内/校园网下载稳定性
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

# CRAN 依赖列表 (数据处理、正则映射、矩阵计算、统计可视化与 Rmd 工具)
cran_pkgs <- c(
  "tidyverse", "data.table", "RColorBrewer", "ggforce", "readxl",
  "openxlsx", "hrbrthemes", "gplots", "fuzzyjoin", "gridExtra",
  "reshape2", "showtext", "plyr", "cowplot", "patchwork", "here",
  "rmarkdown", "jsonlite", "matrixStats", "scales"
)

# Bioconductor 3.23 依赖列表 (热图、质谱数据读取、并行与本地 MassBank CompDb 注释)
bioc_pkgs <- c(
  "Biobase", "ComplexHeatmap", "xcms", "MSnbase", "mzR",
  "BiocParallel", "CompoundDb", "AnnotationHub"
)

check_status <- function(pkgs) {
  sapply(pkgs, requireNamespace, quietly = TRUE)
}

cran_status <- check_status(cran_pkgs)
bioc_status <- check_status(bioc_pkgs)

cat("=== CRAN 依赖检查 ===\n")
print(cran_status)
cat("\n=== Bioconductor 依赖检查 ===\n")
print(bioc_status)

missing_cran <- names(cran_status)[!cran_status]
missing_bioc <- names(bioc_status)[!bioc_status]

cat(sprintf("\n待安装 CRAN (%d): %s\n", length(missing_cran), ifelse(length(missing_cran) == 0, "无", paste(missing_cran, collapse = ", "))))
cat(sprintf("待安装 Bioconductor (%d): %s\n", length(missing_bioc), ifelse(length(missing_bioc) == 0, "无", paste(missing_bioc, collapse = ", "))))

if (dry_run) {
  cat("\n[--dry-run 模式] 仅检查状态，未执行安装。\n")
  quit(save = "no", status = 0)
}

if (length(missing_cran) > 0) {
  cat("\n>>> 开始安装缺失 CRAN 包...\n")
  install.packages(missing_cran, dependencies = TRUE)
}

if (length(missing_bioc) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  cat("\n>>> 开始安装缺失 Bioconductor 包...\n")
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

cat("\n=== 安装后最终状态验证 ===\n")
print(check_status(c(cran_pkgs, bioc_pkgs)))
