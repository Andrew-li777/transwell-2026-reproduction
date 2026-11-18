# ==============================================================================
# 00_download_zenodo_data.R
# 数据获取脚本：从公开 Zenodo Record 22044327 拉取 26 个输入文件至 input_folder/
#
# 上游源码对应:
#   upstream_transwell_2026/download_zenodo_inputs.R + zenodo_utils.R
# 关键修复说明:
#   - 上游 download_zenodo_inputs.R 硬编码了审稿期草稿链接 `/draft?preview=1&token=...`，
#     随 Zenodo 记录正式公开后已返回 HTTP 404；本脚本显式改用公开记录 URL
#     `https://zenodo.org/records/22044327` 复用上游 download_zenodo_inputs() 工具函数完成下载。
# 用法:
#   Rscript reproduction/metabolomics/transwell_2026/repro_scripts/00_download_zenodo_data.R [--overwrite]
# ==============================================================================

# 解析 --file= 定位 transwell_2026/ 项目根目录
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
utils_path <- file.path(proj_root, "upstream_transwell_2026", "zenodo_utils.R")
if (!file.exists(utils_path)) {
  stop("未找到上游 zenodo_utils.R: ", utils_path)
}

# 载入上游 Zenodo API 解析与文件校验下载函数
source(utils_path)

dest_dir <- file.path(proj_root, "input_folder")
dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)

# 195.4 MB 全量拉取在常规网络环境下耗时较长，将 HTTP 超时阈值放宽至至少 600 秒
options(timeout = max(600, getOption("timeout")))

cat(">>> 正在从公开 Zenodo 记录 22044327 下载输入数据集至:", dest_dir, "\n")
downloaded <- download_zenodo_inputs(
  identifier = "https://zenodo.org/records/22044327",
  destination_dir = dest_dir,
  overwrite = "--overwrite" %in% commandArgs(trailingOnly = TRUE)
)

cat(sprintf("<<< 下载校验完成！input_folder/ 共就绪 %d 个文件。\n", length(downloaded)))
