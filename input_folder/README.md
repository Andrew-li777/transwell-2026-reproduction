# `transwell_2026` 输入数据目录说明 (`input_folder/`)

本目录存放通过 [00_download_zenodo_data.R](../repro_scripts/00_download_zenodo_data.R) 从公开 Zenodo 归档 [**Record 22044327**](https://zenodo.org/records/22044327) 下载的 **26 个处理后原始输入文件**（合计约 `195.4 MB`），以及 Step 02 (`02_fig2a_targeted_drug_screen.R`) 运行时自动同步至此以兼容 Step 03 读取的 **`AB005_results_p02_10percent_241126.csv`**（共 27 个数据文件）。

> **Git 版本控制说明**：根据上级 `.gitignore` 规则，本目录下的 `.csv` / `.rds` / `.txt` / `.xlsx` 大体积输入文件会被自动忽略以保持仓库轻量化，仅本说明书 (`README.md`) 与 `.gitkeep` 纳入 Git 追踪。如目录为空，请运行 `Rscript reproduction/metabolomics/transwell_2026/repro_scripts/00_download_zenodo_data.R` 自动拉取。

---

## 1. 全量文件清单、规格与字段结构

### 1.1 Caco-2 细胞活力、TEER 跨上皮电阻与标志物表（支撑 Fig 1）

| 文件名 | 大小 | 行 × 列维度 | 核心列名 / 结构说明 | 消费脚本 |
|---|---:|---:|---|---|
| `Caco2_viability_data.csv` | 6.6 KB | `102 × 8` | `Well`, `value`, `experiment`, `condition` (`GMM` 与 10 株菌), `gmm_value`, `rel_viab`, `bacterial_time` (`16_hour`/`24_hour`) | `00_smoke_test.R`, `01` |
| `Caco2_TEER.csv` | 10.7 KB | `308 × 5` | `condition`, `TEER` ($\Omega$), `name` (`T0`/`T4_HBSS`), `exp_group` (`group_16H`/`group_24H`) | `00_smoke_test.R`, `01` |
| `AB012I_cellviab_24well.txt` | 1.3 KB | `87 × 2` (含 39 行仪器导出头) | 酶标仪原始导出文本，跳过前 39 行后含 `group` 与 `X560.590` 荧光信号 | `01_fig1_caco2_teer_markers.R` |
| `AB012I_TEER_24well.xlsx` | 10.3 KB | `72 × 4` | `well`, `TEER`, `group`, `measurement` (`before`/`after_4H`/`after_8H`) | `01_fig1_caco2_teer_markers.R` |
| `Data_AB012I_markers_transwell1.csv` | 202.7 KB | `439 × 85` | Agilent MassHunter 导出双行表头格式：`Sample_name` + 7 种参照标志物（`antipyrine`, `nadolol`, `terbutaline`, `propanolol`, `haloperidol`, `etoposide`, `warfarin`）等的 `_RT` 与 `_Area` 列 | `01_fig1_caco2_teer_markers.R` |
| `markers_AB012G1_AB012G5_AB012G6.csv` | 214.4 KB | `1,077 × 39` | 补充实验组 AB012G1/G5/G6 通透性标志物定量矩阵（原著 `Fig1_marker_analysis.R` 备用） | 原著参考对照 |
| `markers_AB012G2_AB012G3.csv` | 160.3 KB | `722 × 39` | 补充实验组 AB012G2/G3 通透性标志物定量矩阵（原著 `Fig1_marker_analysis.R` 备用） | 原著参考对照 |

### 1.2 靶向药物菌群代谢筛选表（支撑 Fig 2A）

| 文件名 | 大小 | 行 × 列维度 | 核心列名 / 结构说明 | 消费脚本 |
|---|---:|---:|---|---|
| `Fig2A_drugs_all_raw_with_metadata.rds` | 3.23 MB | `14,008 × 18` (`my_df`) | R `.RData` 归档（`load()` 载入 `my_df`），含 `file`, `column`, `rt`, `rtmin`, `rtmax`, `into`, `intb`, `drug`, `sample_name`, `pool`, `strain`, `TPs`, `well_location` 等 | `02_fig2a_targeted_drug_screen.R` |
| `NEW_pooling_scheme_info_with_IS.csv` | 56.7 KB | `487 × 35` | 481 种临床药物 + 6 种同位素内标（`IS_*`）在 `pool1–pool26` 中的二进制分配矩阵（`1`/`0`）及理论保留时间 `RT.NR`（min）和单同位素质量 `M.SD` | `02`, `03` |
| `AB005_results_p02_10percent_241126.csv` | 41.9 KB | `399 × 7` | **运行期派生**：由 `02` 生成并同步至本目录；列含 `Product.Name`, `strain`, `FC`, `FC_pval`, `padj_strain`, `TP`, `percentage` | `03_fig2be_untargeted_mzmine.R` |

### 1.3 11 株肠道菌非靶向 MZmine 特征定量矩阵（支撑 Fig 2B–E）

全部 `MZMINE_<strain>__iimn_gnps_quant.csv` 文件采用统一的 **MZmine IIMN-GNPS 特征表格式**：前 13 列为固定峰元数据（`row ID`, `row m/z`, `row retention time`, `row ion mobility`, `row ion mobility unit`, `row CCS`, `correlation group ID`, `annotation network number`, `best ion`, `auto MS2 verify`, `identified by n=`, `partners`, `neutral M mass`），第 14 列起为各 `AB005_<strain>_pool<X>_TP<Y>` 样本峰面积。

| 文件名 | 对应菌株 | 大小 (MB) | 特征行 × 总列数 | 消费脚本 |
|---|---|---:|---:|---|
| `MZMINE_Anaes__iimn_gnps_quant.csv` | *Actinomyces naeslundii* (`Anaeslundi`) | 12.81 MB | `13,069 × 119` | `03` |
| `MZMINE_Aomni__iimn_gnps_quant.csv` | *Actinomyces omnicolens* (`Aomnicolens`) | 17.39 MB | `17,725 × 119` | `03` |
| `MZMINE_Btheta__iimn_gnps_quant.csv` | *Bacteroides thetaiotaomicron* (`Btheta`) | 11.75 MB | `13,055 × 110` | `03` |
| `MZMINE_Buni__iimn_gnps_quant.csv` | *Bacteroides uniformis* (`Buni`) | 15.65 MB | `17,192 × 126` | `03` |
| `MZMINE_Ccomes__iimn_gnps_quant.csv` | *Coprococcus comes* (`Ccomes`) | 25.03 MB | `25,011 × 127` | `03` |
| `MZMINE_Cramosum__iimn_gnps_quant.csv` | *Clostridium ramosum* (`Cramosum`) | 20.92 MB | `21,124 × 125` | `03` |
| `MZMINE_Cscindens__iimn_gnps_quant.csv` | *Clostridium scindens* (`Cscindens`) | 12.82 MB | `13,045 × 119` | `03` |
| `MZMINE_Dformi__iimn_gnps_quant.csv` | *Dorea formicigenerans* (`Dformi`) | 12.45 MB | `12,587 × 119` | `03` |
| `MZMINE_Ef2__iimn_gnps_quant.csv` | *Enterococcus faecalis* (`Ef2`) | 11.06 MB | `11,343 × 119` | `03` |
| `MZMINE_Ghaemo__iimn_gnps_quant.csv` | *Gemella haemolysans* (`Ghaemolysans`) | 15.03 MB | `15,228 × 119` | `03` |
| `MZMINE_nobug__iimn_gnps_quant.csv` | 无菌培养基对照 (`nobug`) | 15.58 MB | `17,219 × 118` | `03` |

### 1.4 粪便验证、Transwell 药物转运与 Fig 4 双室非靶向表（支撑 Fig 2E–F, Fig 3, Fig 4）

| 文件名 | 大小 | 行 × 列维度 | 核心列名 / 结构说明 | 消费脚本 |
|---|---:|---:|---|---|
| `Data_AB012I_drugassay3.csv` | 42.4 KB | `96 × 85` | 人源供体 (`MB002`/`MB003`/`MB005`) 与小鼠 (`mouse`) 粪便菌群体外代谢 Roxatidine / Nicergoline / Bisacodyl 靶向矩阵（1 个前置元数据列） | `04` |
| `AB012_drugs_mets_table.csv` | 239.7 KB | `883 × 39` | `B. uniformis` (`Bu`) 与 `nobug` 孵育液中 Roxatidine acetate / Deacetyl-roxatidine / Bisacodyl 的 Apical (`api`) 与 Basolateral (`baso`) 转运表（5 个前置元数据列，样本名在第 1 行标记 `"Name"`） | `04` |
| `Data_AB012I_transwell_drugs.csv` | 713.9 KB | `1,599 × 85` | Dexamethasone 及其菌源产物在 `TP16` 下的双室转运定量表 | `04` |
| `mzmine_dmso_api_baso__quant.csv` | 19.63 MB | `35,434 × 102` | **Fig 4 核心矩阵**：DMSO 处理下 9 株菌 + `nobug` 在 `api`（40 列）与 `baso`（40 列）的 35,434 个正离子特征定量表 | `05` |
| `FT7682_ion.csv` | 168.2 KB | `1,115 × 12` | Vanillic acid (`FT7682`) 预先从 `MTBLS12528` 提取的 EIC 积分坐标表（含 `sample_name`, `column`, `rt`, `into`, `intb`, `sn`, `ion` 等） | `05` |
| `FT5225_ion.csv` | 268.5 KB | `1,778 × 12` | 5-Hydroxynorvaline (`FT5225`) 预先提取的 EIC 积分坐标表 | `05` |

---

## 2. 样本列名命名规范速查

- **AB005（非靶向/靶向筛选）**：`AB005_<strain>_pool<N>_TP<0|2|8>`（例如 `AB005_Btheta_pool1_TP8`）
- **AB012 / AB012I（Transwell 转运）**：`<exp>_<strain/group>_<cond/TP>_<timepoint>_<rep/pool>_<side>`，其中 `<side>` 为 `api`（顶侧/肠腔侧）或 `baso`（基底侧/血浆吸收侧）。
