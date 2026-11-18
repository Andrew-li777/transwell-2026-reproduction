# ==============================================================================
# ZimmermannLab/transwell_2026 全流程自动化复现调度器
#
# 设计说明:
#   - 按模块顺序串联 Fig 1 -> Fig 2A -> Fig 2B-E -> Fig 2E-F & Fig 3 -> Fig 4。
#   - 每步作为独立 Rscript 子进程执行并计时，任意一步退出码非 0 即刻中止并报错。
#   - 结束后汇总 CRAN 与 Bioconductor 3.23 全部实际调用软件包版本号写入 repro_results/environment_versions.tsv。
#
# 使用方式:
#   Rscript reproduction/metabolomics/transwell_2026/repro_scripts/run_all.R --mode smoke
#   Rscript reproduction/metabolomics/transwell_2026/repro_scripts/run_all.R --mode phase2
#   Rscript reproduction/metabolomics/transwell_2026/repro_scripts/run_all.R --mode full
# ==============================================================================

# 解析 --mode 参数 (默认 full 全流程 01-05)
args <- commandArgs(trailingOnly = TRUE)
mode <- "full"
if ("--mode" %in% args) {
  idx <- which(args == "--mode")
  if (idx < length(args)) mode <- args[idx + 1]
} else if (any(grepl("^--mode=", args))) {
  mode <- sub("^--mode=", "", args[grep("^--mode=", args)][1])
}

# 动态定位 transwell_2026/ 项目根目录
get_project_root <- function() {
  c_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", c_args, value = TRUE)
  if (length(file_arg) > 0) {
    script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE))
    return(normalizePath(file.path(script_dir, ".."), mustWork = TRUE))
  }
  getwd()
}

proj_root <- get_project_root()
setwd(proj_root)

cat("===================================================================\n")
cat("  ZimmermannLab/transwell_2026 全流程自动化复现调度器\n")
cat("  工作目录:", proj_root, "\n")
cat("  运行模式:", toupper(mode), "\n")
cat("===================================================================\n\n")

scripts_to_run <- if (mode == "smoke") {
  c("00_smoke_test.R")
} else if (mode == "phase2") {
  c("05_fig4_api_baso_eic_and_stats.R")
} else {
  c(
    "01_fig1_caco2_teer_markers.R",
    "02_fig2a_targeted_drug_screen.R",
    "03_fig2be_untargeted_mzmine.R",
    "04_fig3_transwell_transport.R",
    "05_fig4_api_baso_eic_and_stats.R"
  )
}

status_records <- list()
for (s in scripts_to_run) {
  s_path <- file.path("repro_scripts", s)
  cat(sprintf(">>> 正在执行: %s ...\n", s))
  t_start <- proc.time()
  res <- system2("Rscript", args = c(s_path), stdout = "", stderr = "")
  elapsed <- round((proc.time() - t_start)[3], 1)
  status_str <- if (res == 0) "SUCCESS" else "FAILED"
  cat(sprintf("<<< 完成: %s | 状态: %s | 耗时: %s 秒\n\n", s, status_str, elapsed))
  status_records[[s]] <- data.frame(Script = s, Status = status_str, Elapsed_sec = elapsed)
  if (res != 0) stop("调度流程中止于: ", s)
}

cat("=== 收集并导出 transwell_2026 环境软件包版本快照 ===\n")
# 覆盖 00-05 复现流程实际直接或间接加载的全部 CRAN 与 Bioconductor 核心依赖
bioc_set <- c("Biobase", "ComplexHeatmap", "xcms", "MSnbase", "mzR", "BiocParallel", "CompoundDb", "AnnotationHub")
pkgs <- c(
  "R", "tidyverse", "data.table", "RColorBrewer", "ggforce", "readxl",
  "openxlsx", "hrbrthemes", "gplots", "fuzzyjoin", "gridExtra",
  "reshape2", "showtext", "plyr", "cowplot", "patchwork", "here",
  "rmarkdown", "jsonlite", "matrixStats", "scales", bioc_set
)
ver_list <- lapply(pkgs, function(p) {
  v <- if (p == "R") paste0(R.version$major, ".", R.version$minor) else tryCatch(as.character(packageVersion(p)), error = function(e) "Not Installed")
  src <- if (p == "R") "Base R" else if (p %in% bioc_set) "Bioconductor 3.23" else "CRAN"
  data.frame(Package = p, Version = v, Repository = src, stringsAsFactors = FALSE)
})
ver_df <- do.call(rbind, ver_list)
dir.create("repro_results", showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)
write.table(ver_df, file.path("repro_results", "environment_versions.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

sink(file.path("logs", "sessionInfo.txt"))
cat("Generated at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Mode:", mode, "\n\n")
print(do.call(rbind, status_records))
cat("\n=== R sessionInfo() ===\n")
print(sessionInfo())
sink()

cat("调度运行全部成功完成！\n")
