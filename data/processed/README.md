# `transwell_2026` 中间处理产物目录 (`data/processed/`)

本目录存放复现流水线在各阶段动态生成、缓存并供下游脚本接力的中间数据对象，用以修复上游仓库断链（如原著注释未保存的 `TP8_ttest_results.RDS`）并加速二次复跑。

> **Git 策略**：为避免将高达 `2.27 GB` 的中间 `.RDS` 与大体积 `.csv` 提交至 Git，本目录中除 `README.md` 外的所有数据产物均受 `.gitignore` 屏蔽，均可通过运行 `repro_scripts/run_all.R --mode full` 一键全量再生。

---

## 中间文件清单与流转关系

| 文件名 | 格式 | 实测大小 | 维度 / 内部结构 | 生产脚本 $\rightarrow$ 消费脚本 | 工程说明 |
|---|---|---:|---|---|---|
| `02_fig2a_norm_df.rds` | R `.rds` | 248.1 KB | `data.frame`: `9,984 × 9` (`sample_name`, `Product.Name`, `rt`, `intb`, `pool`, `strain`, `TPs`, `well_location`, `correctedArea`) | `02_fig2a_targeted_drug_screen.R` $\rightarrow$ 存档/延伸检验 | 内标中位数比例校正（已剔除 `IS_sulfamethoxazole`）后的 `AB005` 靶向药物峰面积全表 |
| `TP8_ttest_results.RDS` | R `.RDS` | **2.27 GB** (`2,439,559,313 B`) | 嵌套 `list`：11 株菌 × ~449 种药物 × `11,343–25,011` 特征矩阵（列含 `pvalue`, `statistic`, `padj`, `FoldChange`） | `03_fig2be_untargeted_mzmine.R`（自产自缓存） | **关键断链修复**：原著 `Fig2B_E_untar_analysis.Rmd` 注释了保存代码却直接读取；此处通过 `matrixStats::colVars` 向量化 Welch t 检验现场计算并缓存 |
| `240607_mzmine_candidates.csv` | CSV | 265.9 KB | `1,248 × 15` (`value`, `strain`, `drug`, `padj`, `FoldChange`, `row m/z`, `row retention time`, `deltaMZ`, `deltaRT`, `M.SD`, `RT.NR`, `pool_percentage`, `n_pools_det`, `combo`, `category`) | `03_fig2be_untargeted_mzmine.R` $\rightarrow$ 归档 | 经 4 级漏斗筛选并联合 `AB005` 降解清单标记后的 `123` 个 `Drug_Metabolite` + `1,125` 个 `Physiological_Metabolite` 候选全表（同 `repro_results/03_fig2be_mzmine_candidates.tsv` 内容） |
| `AB007_api_ttest_candidates.csv` | CSV | 1.57 MB | `10,394 × 11` (`strain`, `ion`, `pval`, `statistic`, `padj`, `FC`, `row ID`, `row m/z`, `row retention time`, `neutral M mass`, `massbank_ms1_match`) | `05_fig4_api_baso_eic_and_stats.R` | DMSO 顶侧 (`api`) `FC > 5 & padj < 0.05` 菌源显著离子表（含 `4,124` 唯一特征，其中 `945` 个获本地 MassBank MS1 注释） |
| `AB007_baso_ttest_candidates.csv` | CSV | 543.8 KB | `3,178 × 11`（同上 11 列） | `05_fig4_api_baso_eic_and_stats.R` | DMSO 基底侧 (`baso`) `FC > 3 & padj < 0.05` 菌源显著离子表（含 `1,821` 唯一特征，其中 `789` 个获本地 MassBank MS1 注释） |

---

## 关联文件说明

- 上游跨步骤依赖表 `AB005_results_p02_10percent_241126.csv`（`399 × 7`）由 `02_fig2a_targeted_drug_screen.R` 输出至 `../input_folder/`（为兼容原著默认读入路径），并同步备份于 `../repro_results/02_fig2a_degraded_drugs_results.tsv`。
