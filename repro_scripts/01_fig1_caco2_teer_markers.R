# ==============================================================================
# 01_fig1_caco2_teer_markers.R
# 复现 ZimmermannLab/transwell_2026 Figure 1 (上皮细胞屏障完整性与基准通透性验证):
#
# 上游源码对应:
#   合并并重构上游 3 个独立脚本:
#   - Fig1_TEER_viability_script.R (Fig 1A-B 96 孔筛选细胞活力与 TEER)
#   - Fig1E_Cell_viability&TEER.R  (Fig 1E 24 孔 Transwell 恢复实验，移除冗余 xlsx 依赖)
#   - Fig1F_markers_transwell1.R   (Fig 1F 7 种参考化合物双室转运与 Peff 计算)
# 输入文件:
#   - input_folder/Caco2_viability_data.csv, input_folder/Caco2_TEER.csv
#   - input_folder/AB012I_cellviab_24well.txt, input_folder/AB012I_TEER_24well.xlsx
#   - input_folder/Data_AB012I_markers_transwell1.csv
# 输出产物:
#   - repro_results/01_fig1_viability_and_teer.png
#   - repro_results/01_fig1_ab012i_viability_teer.png
#   - repro_results/01_fig1_markers_and_peff.png & 01_fig1_peff_summary.tsv
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(data.table)
  library(RColorBrewer)
  library(ggforce)
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

strain_lab <- c(
  "An" = "A. naeslundii",
  "Ao" = "A. omnicolens",
  "Bt" = "B. thetaiotaomicron",
  "Bu" = "B. uniformis",
  "Cc" = "C. comes",
  "Cr" = "C. ramosum",
  "Cs" = "C. scindens",
  "Df" = "D. formicigenerans",
  "Ef" = "E. faecalis",
  "Gh" = "G. haemolysans",
  "GMM" = "media control"
)

cat(">>> [1/4] 绘制 Caco-2 11 株菌上清暴露细胞活力与 TEER (Fig 1A-B)...\n")
df_viab <- read.csv(file.path(in_dir, "Caco2_viability_data.csv"))
p_viab <- ggplot(df_viab, aes(x = condition, y = rel_viab, fill = bacterial_time)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(group = bacterial_time), position = position_jitterdodge(jitter.width = 0.05), size = 0.7) +
  ylab("Relative fluo 560/590") +
  labs(title = "Caco-2 Viability Across Gut Bacterial Supernatants", fill = "Bacterial incubation time") +
  scale_fill_manual(values = c("16_hour" = "#FF7F7F", "24_hour" = "#E75480")) +
  scale_x_discrete(labels = strain_lab) +
  theme_minimal() +
  cell_plot_theme +
  theme(
    legend.position = "top",
    legend.key.size = unit(0.35, "cm"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 1, face = "italic")
  )

# 剔除暴露前即处于低阻值状态的异常 C. ramosum (Cr) 孔 (< 1000 Ohm，保证单层细胞起始屏障完整)
df_teer <- read.csv(file.path(in_dir, "Caco2_TEER.csv")) %>%
  filter(!(condition == "Cr" & name == "T0" & TEER < 1000)) %>%
  mutate(
    interaction_factor = factor(
      paste(condition, name, exp_group, sep = "_"),
      levels = unlist(lapply(
        c("GMM", "An", "Ao", "Bt", "Bu", "Cc", "Cr", "Cs", "Df", "Ef", "Gh"),
        function(cond) {
          c(
            paste(cond, "T0", "group_16H", sep = "_"),
            paste(cond, "T0", "group_24H", sep = "_"),
            paste(cond, "T4_HBSS", "group_16H", sep = "_"),
            paste(cond, "T4_HBSS", "group_24H", sep = "_")
          )
        }
      ))
    ),
    condition_factor = factor(condition, levels = c("GMM", "An", "Ao", "Bt", "Bu", "Cc", "Cr", "Cs", "Df", "Ef", "Gh"))
  )

p_teer <- ggplot(df_teer, aes(x = interaction_factor, y = TEER, fill = name)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(group = name), position = position_jitterdodge(jitter.width = 0.05), size = 0.6) +
  scale_fill_manual(
    labels = c("T0" = "0 hours", "T4_HBSS" = "4 hours"),
    values = c("T0" = "#FFC0CB", "T4_HBSS" = "#5F9EA0")
  ) +
  facet_grid(~condition_factor, scales = "free_x", space = "free_x") +
  scale_x_discrete(labels = function(x) sub("group_", "", sub(".*_(group_16H|group_24H)$", "\\1", x))) +
  labs(
    title = "Caco-2 TEER Before & After 4h Supernatant Exposure",
    fill = "Exposure time",
    x = "Bacterial species & cultivation time",
    y = "TEER (Ohm)"
  ) +
  theme_minimal() +
  cell_plot_theme +
  theme(
    legend.position = "top",
    legend.key.size = unit(0.35, "cm"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 0.8)
  )

ggsave(file.path(res_dir, "01_fig1_viability_and_teer.png"), p_viab / p_teer, width = 10, height = 7.5, dpi = 300)

cat(">>> [2/4] 绘制 AB012I 24 孔 Transwell 细胞活力与 TEER 恢复 (Fig 1E, readxl 无 Java)...\n")
# Tecan 酶标仪导出的 .txt 含 39 行仪器元信息头，skip = 39 后方为标准逗号分隔表头
df_viab_i <- read.table(file.path(in_dir, "AB012I_cellviab_24well.txt"), skip = 39, sep = ",", header = TRUE) %>%
  na.omit() %>%
  mutate(gmm_value = mean(X560.590[group == "GMM"]), rel_viab = X560.590 / gmm_value)

p_viab_i <- ggplot(df_viab_i, aes(x = group, y = rel_viab, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 0.9) +
  ylab("Relative signal 560/590") +
  labs(title = "AB012I 24-well Caco-2 Viability") +
  scale_fill_manual(values = c("#DDA0DD", "#5F9EA0", "#E75480", "#F08080", "#FFD1DC", "#FFA07A")) +
  theme_minimal() +
  cell_plot_theme +
  theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"))

df_teer_i <- read_excel(file.path(in_dir, "AB012I_TEER_24well.xlsx"))
df_teer_i$measurement <- factor(df_teer_i$measurement, levels = c("before", "after_4H", "after_8H"))

p_teer_i <- df_teer_i %>%
  ggplot(aes(x = measurement, y = TEER, fill = measurement)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 0.7) +
  ylab("Electrical resistance (Ohm)") +
  labs(title = "AB012I TEER Exposure (4h) & Recovery (8h)") +
  facet_wrap(~group, nrow = 1) +
  scale_fill_manual(
    labels = c("before" = "0 h", "after_4H" = "4 h", "after_8H" = "8 h recovery"),
    values = c("before" = "#FFC0CB", "after_4H" = "#5F9EA0", "after_8H" = "orange")
  ) +
  scale_x_discrete(labels = c("before" = "0h", "after_4H" = "4h", "after_8H" = "8h")) +
  theme_minimal() +
  cell_plot_theme +
  theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(res_dir, "01_fig1_ab012i_viability_teer.png"), p_viab_i | p_teer_i + plot_layout(widths = c(1, 2)), width = 10.5, height = 4.2, dpi = 300)

cat(">>> [3/4] 解析 Transwell 标志物定性定量表与计算表观渗透系数 Peff (Fig 1F)...\n")
# Agilent MassHunter 导出 CSV 为双行表头，从 "Results" 列名动态反推各药物对应的 (_RT, _Area) 两列
raw_m <- read.csv(file.path(in_dir, "Data_AB012I_markers_transwell1.csv"))
target_drugs <- gsub(" ", "_", gsub(" Results", "", grep("Results", colnames(raw_m), fixed = TRUE, value = TRUE), fixed = TRUE), fixed = TRUE)
if (any(raw_m[1] == "!", na.rm = TRUE)) raw_m <- raw_m[, -c(1, 2)]
colnames(raw_m) <- as.character(c("Sample_name", as.vector(sapply(target_drugs, function(x) paste0(x, c("_RT", "_Area"))))))
merged_df <- raw_m %>%
  slice(-1) %>%
  mutate(across(everything(), ~ na_if(as.character(.), "")))
merged_df$group <- sapply(strsplit(merged_df$Sample_name, "_"), "[", 2)
merged_df$TP <- sapply(strsplit(merged_df$Sample_name, "_"), "[", 3)
merged_df$rep <- sapply(strsplit(merged_df$Sample_name, "_"), "[", 4)
merged_df$side <- sapply(strsplit(merged_df$Sample_name, "_"), "[", 5)
merged_df$exp <- sapply(strsplit(merged_df$Sample_name, "_"), "[", 1)

area_cols <- grep("_Area$", colnames(merged_df))
merged_df[, area_cols] <- lapply(merged_df[, area_cols], function(x) suppressWarnings(as.numeric(x)))
merged_df$time <- suppressWarnings(as.numeric(gsub("T", "", merged_df$TP)))
merged_df <- merged_df[!grepl("waterwash|IS|Water", merged_df$Sample_name), ]

tmp <- merged_df %>%
  dplyr::select(Sample_name, group, exp, time, rep, side, all_of(grep("marker.*_Area$", colnames(merged_df)))) %>%
  pivot_longer(cols = grep("marker.*_Area$", colnames(.)), names_to = "name", values_to = "value") %>%
  filter(!is.na(time), !is.na(value))

drug_labs <- c(
  "marker_antipyrine.Results_Area" = "Antipyrine",
  "marker_nadolol.Results_Area" = "Nadolol",
  "marker_terbutaline.Results_Area" = "Terbutaline",
  "marker_propanolol.Results_Area" = "Propranolol",
  "marker_haloperidol.Results_Area" = "Haloperidol",
  "marker_etoposide.Results_Area" = "Etoposide",
  "marker_warfarin.Results_Area" = "Warfarin"
)

p_gmm_markers <- tmp %>%
  filter(group == "GMM") %>%
  ggplot(aes(x = time, y = value, color = side)) +
  geom_point(size = 0.8, alpha = 0.8) +
  geom_smooth(se = FALSE, method = "loess", formula = y ~ x, linewidth = 0.8, aes(group = interaction(side, name))) +
  labs(title = "GMM Media Control: 7 Permeability Markers (Apical vs Basolateral)", x = "Time (minutes)", y = "Peak Area Intensity") +
  theme_minimal() +
  scale_color_manual(values = c("api" = "#4A90E2", "baso" = "#E75480"), labels = c("api" = "Apical", "baso" = "Basolateral")) +
  scale_x_continuous(breaks = c(0, 120, 240)) +
  facet_wrap(~name, scales = "free_y", ncol = 2, labeller = labeller(name = drug_labs)) +
  scale_y_continuous(labels = scales::scientific, n.breaks = 3) +
  cell_plot_theme +
  theme(legend.position = "top")

# 计算表观渗透系数 Peff = [V / (A * C0)] * (dC_baso / dt)：
#   volume = 100 uL (基底侧接收腔体积), surface = 0.11 cm^2 (96孔 Transwell 滤膜面积),
#   选取线性稳态区间 [30 min, 150 min] 估算基底侧累积斜率
calculate_peff <- function(df_long, marker_name, grp = "GMM", side_baso = "baso", side_donor = "api",
                           time_start = 30, time_end = 150, volume = 100, surface = 0.11) {
  baso_data <- df_long %>% filter(group == grp, side == side_baso, name == marker_name)
  donor_conc <- mean(df_long$value[df_long$group == grp & df_long$side == side_donor & df_long$name == marker_name & df_long$time == 0], na.rm = TRUE)
  dC_dt <- (mean(baso_data$value[baso_data$time == time_end], na.rm = TRUE) -
              mean(baso_data$value[baso_data$time == time_start], na.rm = TRUE)) / (time_end - time_start)
  (volume / (surface * donor_conc)) * dC_dt
}

# 7 种参照标志物的人体口服吸收率文献真值 FA% (用于检验体外 ln(Peff) 相关性)
fa_map <- c(
  "marker_antipyrine.Results_Area" = 98,
  "marker_etoposide.Results_Area" = 50,
  "marker_haloperidol.Results_Area" = 70,
  "marker_nadolol.Results_Area" = 35,
  "marker_propanolol.Results_Area" = 90,
  "marker_terbutaline.Results_Area" = 62,
  "marker_warfarin.Results_Area" = 98
)

peff_table <- tmp %>%
  distinct(name) %>%
  mutate(
    Marker = recode(name, !!!drug_labs),
    Peff = map_dbl(name, ~ calculate_peff(tmp, .x)),
    Human_FA_pct = unname(fa_map[name])
  ) %>%
  arrange(Peff)

p_peff <- ggplot(peff_table %>% filter(Peff > 0), aes(x = log(Peff), y = Human_FA_pct, label = Marker)) +
  geom_point(size = 2.8, color = "#5F9EA0") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#E75480", linewidth = 0.7) +
  geom_text(
    aes(vjust = ifelse(Marker == "Warfarin", 1.5, -0.7),
        hjust = ifelse(Marker %in% c("Antipyrine", "Warfarin", "Propranolol", "Haloperidol"), 0.85, 0.25)),
    size = 2.8
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.12, 0.12))) +
  scale_y_continuous(limits = c(15, 108)) +
  labs(title = "Measured log(Peff) vs Human Fraction Absorbed", x = "ln(Peff)", y = "Human Oral FA (%)") +
  theme_minimal() +
  cell_plot_theme

ggsave(file.path(res_dir, "01_fig1_markers_and_peff.png"), p_gmm_markers | p_peff, width = 11, height = 6.5, dpi = 300)
write.table(peff_table, file.path(res_dir, "01_fig1_peff_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

cat(">>> [4/4] Fig 1 验证指标:\n")
print(peff_table[, c("Marker", "Peff", "Human_FA_pct")])
cat("01_fig1_caco2_teer_markers.R 运行成功。\n")
