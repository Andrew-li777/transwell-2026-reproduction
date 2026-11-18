# ==============================================================================
# 04_fig3_transwell_transport.R
# 复现 ZimmermannLab/transwell_2026 Figure 2E-F & Figure 3:
#
# 上游源码对应:
#   合并并修复上游 4 个独立脚本:
#   - Fig2EF_Drugassay.R     (供体粪便菌群体外代谢验证: Roxatidine / Nicergoline / Bisacodyl)
#   - Fig3_plots.R           (AB012 双室药物与代谢物转运)
#   - Fig3_bisacodyl.R       (Roxatidine acetate / Deacetyl-roxatidine 对比)
#   - Fig3D_transwell_drugs.R(Dexamethasone 菌群代谢后上皮渗透性变化)
# 输入文件:
#   - input_folder/Data_AB012I_drugassay3.csv
#   - input_folder/AB012_drugs_mets_table.csv
#   - input_folder/Data_AB012I_transwell_drugs.csv
# 输出产物:
#   - repro_results/04_fig2ef_microbiota_drug_assays.png
#   - repro_results/04_fig3_bisacodyl_and_dexamethasone.png
#   - repro_results/04_fig3_dexamethasone_summary.tsv
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(RColorBrewer)
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
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

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

# 动态靶向表头解析器：兼容不同批次导出的 CSV（例如 AB012 有 5 个前置元数据列，而 Data_AB012I 仅 1 列），
# 根据首见 "Results" 列反推目标物数量，并依据第 1 行 "Name" 标记精确定位样本名称列
parse_targeted_csv <- function(csv_path, sample_col_out = "Sample_name") {
  raw <- read.csv(csv_path, stringsAsFactors = FALSE)
  if (any(raw[[1]] == "!", na.rm = TRUE)) raw <- raw[, -c(1, 2)]
  first_res <- grep("Results", colnames(raw), fixed = TRUE)[1]
  t_drugs <- gsub(" ", "_", gsub(" Results", "", grep("Results", colnames(raw), fixed = TRUE, value = TRUE), fixed = TRUE), fixed = TRUE)
  n_meta <- first_res - 1
  meta_names <- if (n_meta == 1) {
    sample_col_out
  } else {
    c(sample_col_out, paste0("meta_", seq_len(n_meta - 1)))
  }
  # 若第 1 行标明 "Name" 在非第 1 列则自动对齐样本名列
  row1_chr <- as.character(unlist(raw[1, seq_len(n_meta)]))
  name_idx <- which(row1_chr == "Name")
  if (length(name_idx) == 1) meta_names[name_idx] <- sample_col_out
  colnames(raw) <- as.character(c(meta_names, as.vector(sapply(t_drugs, function(x) paste0(x, c("_RT", "_Area"))))))
  df <- raw %>%
    slice(-1) %>%
    mutate(across(everything(), ~ na_if(as.character(.), "")))
  area_cols <- grep("_Area$", colnames(df))
  df[, area_cols] <- lapply(df[, area_cols], function(x) suppressWarnings(as.numeric(x)))
  df
}

cat(">>> [1/3] 复现 Figure 2E-F 粪便菌群体外药代转化验证 (Roxatidine / Nicergoline / Bisacodyl)...\n")
df_assay <- parse_targeted_csv(file.path(in_dir, "Data_AB012I_drugassay3.csv"), sample_col_out = "Name")
df_assay$group <- sapply(strsplit(df_assay$Name, "_"), "[", 2)
df_assay$condition <- sapply(strsplit(df_assay$Name, "_"), "[", 3)
df_assay$TPs <- sapply(strsplit(df_assay$Name, "_"), "[", 4)
df_assay <- df_assay %>% filter(TPs != "all")

prep_assay <- df_assay %>%
  dplyr::select(Name, condition, TPs, group, ends_with("_Area")) %>%
  pivot_longer(cols = ends_with("_Area"), names_to = "drug", values_to = "Area") %>%
  mutate(TPs = gsub("TP", "", TPs))

strain_assay_labs <- c("GMM" = "Sterile control", "MB002" = "Donor 1", "MB003" = "Donor 2", "MB005" = "Donor 3", "mouse" = "Mouse")

plot_assay_pair <- function(drugs_vec, drug_labs_vec, title_txt) {
  prep_assay %>%
    filter(drug %in% drugs_vec, condition != "DMSO", !is.na(Area)) %>%
    ggplot(aes(x = TPs, y = Area, color = group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(position = position_jitter(width = 0.15), size = 0.8) +
    facet_grid(drug ~ group, scales = "free_y", labeller = labeller(drug = drug_labs_vec, group = strain_assay_labs)) +
    scale_y_continuous(labels = scales::scientific, n.breaks = 3) +
    scale_color_manual(values = c("#DDA0DD", "#5F9EA0", "#E75480", "#F08080", "#FFD1DC", "#FFA07A")) +
    labs(title = title_txt, x = "Time (h)", y = "Peak Area (a.u.)") +
    theme_minimal() +
    cell_plot_theme +
    theme(legend.position = "none")
}

p_rox <- plot_assay_pair(
  c("drug_roxatidine.Results_Area", "met_roxatidine.Results_Area"),
  c("drug_roxatidine.Results_Area" = "Roxatidine acetate", "met_roxatidine.Results_Area" = "Deacetyl-roxatidine"),
  "Fig 2E: Roxatidine Acetate Metabolism Across Human/Mouse Microbiota"
)
p_bis_assay <- plot_assay_pair(
  c("drug_bisacodyl.Results_Area", "met_bisacodyl278.Results_Area", "met_bisacodyl320.Results_Area"),
  c("drug_bisacodyl.Results_Area" = "Bisacodyl", "met_bisacodyl278.Results_Area" = "Metabolite (278 m/z)", "met_bisacodyl320.Results_Area" = "Metabolite (320 m/z)"),
  "Fig 2F: Bisacodyl Hydrolysis Across Donors"
)
ggsave(file.path(res_dir, "04_fig2ef_microbiota_drug_assays.png"), p_rox / p_bis_assay + plot_layout(heights = c(1, 1.3)), width = 9.5, height = 7.5, dpi = 300)

cat(">>> [2/3] 复现 Figure 3 Bisacodyl 与 Roxatidine Transwell 通透性曲线...\n")
df_bisa <- parse_targeted_csv(file.path(in_dir, "AB012_drugs_mets_table.csv"), sample_col_out = "sample_name")
df_bisa$strain <- sapply(strsplit(df_bisa$sample_name, "_"), "[", 2)
df_bisa$condition <- sapply(strsplit(df_bisa$sample_name, "_"), "[", 3)
df_bisa$TPs <- sapply(strsplit(df_bisa$sample_name, "_"), "[", 4)
df_bisa$Rep <- sapply(strsplit(df_bisa$sample_name, "_"), "[", 6)
df_bisa$side <- sapply(strsplit(df_bisa$sample_name, "_"), "[", 7)
df_bisa$pool <- sapply(strsplit(df_bisa$sample_name, "_"), "[", 5)
df_bisa$Time <- suppressWarnings(as.numeric(gsub("T", "", df_bisa$TPs)))

prep_bisa <- df_bisa %>%
  dplyr::select(sample_name, condition, side, Time, Rep, strain, pool, ends_with("_Area")) %>%
  pivot_longer(cols = ends_with("_Area"), names_to = "drug", values_to = "Area") %>%
  filter(!is.na(Time), !is.na(Area))

p_rox_trans <- prep_bisa %>%
  filter((drug == "drug_roxatidine.Results_Area" & strain == "nobug" & pool == "pool2") |
           (drug == "met_roxatidine.Results_Area" & strain == "Bu" & condition == "TP8")) %>%
  ggplot(aes(x = Time, y = Area, color = side, group = side)) +
  geom_point(size = 1.1) +
  geom_smooth(method = "loess", formula = y ~ x, linewidth = 1.1, se = FALSE) +
  facet_grid(
    side ~ drug + strain,
    scales = "free_y",
    labeller = labeller(
      drug = c("drug_roxatidine.Results_Area" = "Roxatidine acetate", "met_roxatidine.Results_Area" = "Deacetyl-roxatidine"),
      strain = c("nobug" = "Sterile control", "Bu" = "B. uniformis"),
      side = c("api" = "Apical", "baso" = "Basolateral")
    )
  ) +
  scale_x_continuous(breaks = c(0, 120, 240)) +
  scale_y_continuous(labels = scales::scientific, n.breaks = 3) +
  scale_color_manual(values = c("api" = "#4A90E2", "baso" = "#E75480"), labels = c("api" = "Apical", "baso" = "Basolateral")) +
  labs(title = "Fig 3: Roxatidine Acetate Activation & Permeation", x = "Time (minutes)", y = "Peak Area Intensity") +
  theme_minimal() +
  cell_plot_theme +
  theme(legend.position = "top")

cat(">>> [3/3] 复现 Figure 3D Dexamethasone Transwell 顶侧/基底侧转运曲线...\n")
df_dex <- parse_targeted_csv(file.path(in_dir, "Data_AB012I_transwell_drugs.csv"), sample_col_out = "Sample_name")
df_dex$group <- sapply(strsplit(df_dex$Sample_name, "_"), "[", 2)
df_dex$TP <- sapply(strsplit(df_dex$Sample_name, "_"), "[", 3)
df_dex$rep <- sapply(strsplit(df_dex$Sample_name, "_"), "[", 6)
df_dex$side <- sapply(strsplit(df_dex$Sample_name, "_"), "[", 7)
df_dex$time <- suppressWarnings(as.numeric(gsub("T", "", sapply(strsplit(df_dex$Sample_name, "_"), "[", 4))))
df_dex$condition <- sapply(strsplit(df_dex$Sample_name, "_"), "[", 5)

prep_dex <- df_dex %>%
  filter(!grepl("water|IS", Sample_name), !is.na(time)) %>%
  dplyr::select(Sample_name, condition, TP, rep, time, group, side, ends_with("_Area")) %>%
  pivot_longer(cols = ends_with("_Area"), names_to = "drug", values_to = "Area") %>%
  filter(drug %in% c("drug_dexamethasone.Results_Area", "met_dexamethasone.Results_Area"), !is.na(Area))

p_dex_trans <- prep_dex %>%
  filter(TP == "TP16", condition == "drugs") %>%
  ggplot(aes(x = time, y = Area, color = side, linetype = group)) +
  geom_point(size = 0.8, alpha = 0.6) +
  geom_smooth(aes(group = interaction(side, group)), method = "loess", formula = y ~ x, linewidth = 1.0, se = FALSE) +
  facet_grid(
    side ~ drug,
    scales = "free_y",
    labeller = labeller(
      drug = c("drug_dexamethasone.Results_Area" = "Dexamethasone", "met_dexamethasone.Results_Area" = "Dexamethasone metabolite"),
      side = c("api" = "Apical", "baso" = "Basolateral")
    )
  ) +
  scale_color_manual(values = c("api" = "#4A90E2", "baso" = "#E75480"), labels = c("api" = "Apical", "baso" = "Basolateral")) +
  scale_x_continuous(breaks = c(0, 120, 240)) +
  scale_y_continuous(labels = scales::scientific, n.breaks = 3) +
  labs(title = "Fig 3D: Dexamethasone Bacterial Metabolism Alters Permeability", x = "Time (minutes)", y = "Peak Area Intensity", linetype = "Microbiota") +
  guides(color = guide_legend(nrow = 1), linetype = guide_legend(nrow = 2)) +
  theme_minimal() +
  cell_plot_theme +
  theme(legend.position = "top", legend.box = "vertical")

ggsave(file.path(res_dir, "04_fig3_bisacodyl_and_dexamethasone.png"), p_rox_trans | p_dex_trans, width = 12, height = 6.2, dpi = 300)

summary_stats <- prep_dex %>%
  filter(TP == "TP16", condition == "drugs") %>%
  group_by(drug, group, side, time) %>%
  summarise(mean_area = mean(Area, na.rm = TRUE), sd_area = sd(Area, na.rm = TRUE), .groups = "drop")
write.table(summary_stats, file.path(res_dir, "04_fig3_dexamethasone_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

cat("04_fig3_transwell_transport.R 运行成功。\n")
