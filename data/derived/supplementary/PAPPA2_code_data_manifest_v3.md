# PAPPA2 论文代码与数据可用性清单 v3

日期：2026-09-30
对应稿件：`submission/manuscript/PAPPA2_manuscript_draft_v6_full_data.md`

## 一、数据来源与登录号

| 数据 | 角色 | 登录号或来源 | 备注 |
|---|---|---|---|
| KPMP snRNA-seq | 发现集 | GEO 超级系列 `GSE183279`，snRNA 子系列 `GSE183277` | KPMP 联盟资源 |
| KPMP Visium | 空间跨模态 | GEO `GSE183456` | 同属 KPMP 超级系列 |
| GSE137570 | 外部公开 eGFR | GEO `GSE137570` | Cohort 1，24 例公开 eGFR |
| Nephroseq v5 | 外部 GFR 关联 | https://www.nephroseq.org | 预计算分析 |
| GSE115098 | 候选外部资源，未使用 | GEO `GSE115098` | 逐样本 GFR 未公开，未纳入分析 |

注意：

- `GSE183456` 不是外部队列，属于 KPMP 同联盟空间资源；
- `GSE115098` 的表达数据公开，但逐样本 GFR 不在公开补充材料中，
  因此本版本未把它写成外部验证。

## 二、主要分析脚本

### 1. 完整 TAL 伪批量与分析单元比较

| 脚本 | 作用 |
|---|---|
| `repo/analysis/extract_full_kpmp_tal_focus_expression.py` | 从完整 KPMP h5ad 提取四个 TAL 亚类和五个焦点基因 |
| `repo/analysis/prepare_full_tal_milo_metadata.py` | 从完整 h5ad 提取 TAL metadata 和 UMAP |
| `repo/analysis/analyze_full_tal_donor_pseudobulk.py` | 完整 TAL 供者伪批量、年龄性别调整和旧子集方向比较 |
| `repo/analysis/analyze_full_tal_pseudobulk_robustness.py` | 留一供者、bootstrap 和皮质/髓质分层 |
| `repo/analysis/analyze_full_tal_nhood_expression.py` | 完整邻域表达关联 |
| `repo/analysis/sensitivity_full_tal_nhood_expression.py` | 完整邻域 bootstrap、LODO 和纯度偏相关 |
| `repo/analysis/analyze_full_tal_within_donor_cell_rank.py` | 完整供者内细胞秩相关 |
| `repo/analysis/summarize_full_tal_analysis_units.py` | 完整四层分析单元一致性表 |

### 2. TAL 亚类邻域分析

| 脚本 | 作用 |
|---|---|
| `scripts/run_milo_tal_kpmp.R` | TAL 亚类 Milo 邻域丰度分析 |
| `scripts/run_pyaugur_tal_kpmp.py` | TAL 亚类 Augur 可分性分析 |
| `scripts/summarize_tal_milo_augur.py` | 汇总 Milo 和 Augur 结果 |
| `scripts/prepare_milo_tal_metadata.py` | 准备 Milo 元数据 |

### 3. 邻域层面表达关联与稳健性

| 脚本 | 作用 |
|---|---|
| `work/analyze_tal_nhood_focus_expression.py` | 邻域平均表达与 Milo logFC 关联 |
| `work/sensitivity_tal_nhood_expression.py` | 供者分层重采样、留一供者、纯度偏相关 |
| `work/validate_milo_seed_export.py` | 固定种子 Milo 结果校验 |

### 4. 供者内细胞层面关联

| 脚本 | 作用 |
|---|---|
| `work/analyze_tal_within_donor_cell_rank.py` | 供者内细胞秩相关与置换检验 |
| `work/analyze_tal_within_donor_nhood_effect.py` | 供者内邻域效应辅助分析 |

### 5. 临床关联

| 脚本 | 作用 |
|---|---|
| `scripts/analyze_pappa2_clinical_egfr.py` | KPMP eGFR 初步关联 |
| `work/analyze_pappa2_egfr_adjusted.py` | KPMP eGFR 有序 logistic 模型 |
| `work/analyze_egfr_gene_specificity.py` | 五基因 eGFR 模型比较 |
| `work/analyze_gse137570_pappa2_egfr.py` | GSE137570 Cohort 1 外部 eGFR 检验 |
| `work/query_nephroseq_gfr_genes.py` | Nephroseq GFR 关联查询 |

### 6. 空间分析

| 脚本 | 作用 |
|---|---|
| `scripts/analyze_pappa2_branch_signatures_spatial.py` | GSE183456 分支签名空间共现 |
| `scripts/summarize_pappa2_23sample_spatial.py` | 23 切片空间结果汇总 |
| `scripts/analyze_pappa2_branch_pathology_regions.py` | 病理区域富集分析 |
| `scripts/analyze_pappa2_directional_spatial.py` | 空间方向性分析 |

### 7. 论文与图版

| 脚本 | 作用 |
|---|---|
| `work/build_manuscript_evidence_matrix.py` | 证据矩阵生成 |
| `repo/manuscript/check_submission_consistency.py` | 投稿包图、表、补充文件和 PDF 一致性检查 |
| `work/make_figures_v5_full.py` | 完整数据 v5 八图版生成 |
| `work/make_figures_v2.py` | 图 1-8 生成 |
| `work/check_figures.py` | 图版非空白检查 |
| `work/search_refs.py` | OpenAlex/CrossRef 文献检索 |
| `work/verify_refs.py` | Europe PMC DOI/PMID 核对 |

## 三、随机种子与关键参数

| 参数 | 值 |
|---|---|
| Milo 固定种子 `MILO_SEED` | `20260927` |
| 供者分层重采样次数 | 500 |
| 供者内细胞置换次数 | 1,000 |
| 完整 TAL donor bootstrap 次数 | 1,000 |
| 完整 TAL donor bootstrap 种子 | `20260930` |
| Augur 子采样次数 | 20 |
| Augur 每子采样细胞数 | 30 |
| Augur 交叉验证折数 | 3 |
| 供者内细胞阈值 | 3、5、10 |
| 主阈值 | 5 cells per donor |
| Milo 阈值 | `SpatialFDR < 0.05` |
| BH 校正范围 | 60 个基因 × 亚类 × 比较检验 |
| 完整 TAL 细胞数 | 106,851 |
| 旧邻域分析子集细胞数 | 4,446 |
| 图版抖动随机种子 | 1 |

## 四、关键结果文件

### 发现集与 TAL 亚类

```text
outputs\KPMP_full_TAL_focus_expression.csv.gz
outputs\KPMP_full_TAL_focus_expression_manifest.csv
outputs\KPMP_full_TAL_focus_gene_pseudobulk.csv
outputs\KPMP_full_vs_subset_TAL_focus_pseudobulk_comparison.csv
outputs\KPMP_TAL_analysis_unit_concordance.csv
outputs\KPMP_TAL_analysis_unit_concordance_summary.md
outputs\KPMP_full_TAL_focus_pseudobulk_lodo.csv
outputs\KPMP_full_TAL_focus_pseudobulk_bootstrap.csv
outputs\KPMP_full_TAL_focus_pseudobulk_region_strata.csv
outputs\KPMP_milo_TAL_subclass_summary.csv
outputs\KPMP_Augur_TAL_full_AUC_summary.csv
outputs\KPMP_TAL_focus_gene_pseudobulk.csv
```

### 邻域层面与供者内层面

```text
outputs\KPMP_TAL_seed20260927_nhood_focus_robustness.csv
outputs\KPMP_TAL_seed20260927_within_donor_cell_rank.csv
outputs\KPMP_TAL_seed20260927_within_donor_cell_rank_threshold_sensitivity.csv
```

### 临床关联

```text
outputs\KPMP_PAPPA2_TAL_egfr_models.csv
outputs\KPMP_PAPPA2_TAL_egfr_gene_specificity.csv
outputs\KPMP_PAPPA2_TAL_egfr_linear_sensitivity.csv
outputs\GSE137570_PAPPA2_eGFR_validation.csv
outputs\nephroseq_GFR_gene_summary.csv
outputs\nephroseq_GFR_significant_analyses.csv
```

### 空间

```text
outputs\KPMP_Visium_PAPPA2_branch_signature_sample_associations.csv
outputs\KPMP_Visium_PAPPA2_branch_signature_group_summary.csv
outputs\PAPPA2_23sample_spatial_validation_summary.md
```

### 图版

```text
outputs\Fig1_study_design_v5.png
outputs\Fig2_TAL_full_neighborhood_heterogeneity_v5.png
outputs\Fig3_analysis_unit_framework_v5.png
outputs\Fig4_full_vs_subset_coverage_v5.png
outputs\Fig5_PAPPA2_case_analysis_v5.png
outputs\Fig6_aTAL2_frTAL_unit_conflict_v5.png
outputs\Fig7_internal_robustness_v5.png
outputs\Fig8_TAL_state_model_v5.png
outputs\Fig1_study_design.png
outputs\Fig2_TAL_neighborhood_remodeling.png
outputs\Fig3_PAPPA2_pseudobulk.png
outputs\Fig4_two_level_association.png
outputs\Fig5_KPMP_eGFR_specificity.png
outputs\Fig6_KPMP_spatial_crossmodality.png
outputs\Fig7_external_GFR_associations.png
outputs\Fig8_TAL_state_model.png
```

一致性报告：

```text
outputs\PAPPA2_submission_consistency_report_v1.md
```

## 五、软件环境

已确认的本机版本：

```text
Python        3.12.13
numpy         2.3.5
pandas        3.0.1
scipy         1.18.0
statsmodels   0.14.6
scikit-learn  1.9.0
matplotlib    3.11.0
scanpy        1.12.3
anndata       0.13.2
pdfplumber    可用
R             4.6.1
```

已确认：

```text
R             4.6.1
miloR         2.8.1
SingleCellExperiment 1.34.0
edgeR         4.10.5
limma         3.68.5
BiocParallel  1.46.0
```

仍需补齐：

```text
Augur 包精确版本（Python 和 R）
```

Milo 运行日志：

```text
work\milo_tal_seed_20260927.log
```

完整 TAL Milo 输出前缀：

```text
KPMP_milo_full_TAL_final
```

完整数据主分析链路支持以下环境变量：

```text
PAPPA2_PROJECT_ROOT
PAPPA2_DATA_ROOT
PAPPA2_OUTPUT_ROOT
PAPPA2_WORK_ROOT
```

## 六、数据与代码可用性声明

```text
All data analyzed in this study are publicly available. KPMP snRNA-seq
and Visium data are available through the GEO superseries GSE183279
(snRNA-seq component GSE183277; Visium component GSE183456). GSE137570
was obtained from GEO. Nephroseq v5 was queried through its public
interface. Sample-level GFR annotations for the GSE115098 tubule cohort
are not publicly available, and that resource was therefore not used for
external validation. Analysis scripts and result files are stored in the
public repository at
https://github.com/ZhiPeng1022/pappa2-tal-analysis-unit-sensitivity and
are archived through the Zenodo concept DOI
https://doi.org/10.5281/zenodo.23196128.
```

Mendelian randomization was not performed because the public protein
resources examined did not contain a usable PAPPA2 cis-pQTL. This is a
data-availability limitation and is not reported as a negative causal
test.

The complete TAL pseudobulk analysis used 106,851 cells from four target
TAL subtypes. Milo, neighborhood association, and within-donor cell-rank
analyses also used the complete TAL set. The 4,446-cell derivative subset
was retained only as a coverage and direction comparator.

## 七、正式归档状态

已完成：

1. 已建立公开 GitHub 仓库并整理脚本、派生表和生成图目录；
2. 已记录 Python、R、Milo、Augur、scipy、statsmodels 等运行环境和种子；
3. 已生成 Zenodo 概念 DOI `10.5281/zenodo.23196128`；
4. 大型派生补充表以 `.csv.gz` 和补充表清单形式归档；
5. 核心脚本已改用 `PAPPA2_*` 环境变量，不再包含个人电脑路径；
6. 已在正文和补充材料中写明公开数据来源、伦理和数据使用边界。

后续可选优化：如果目标期刊要求，可将完整 TAL 分析链合并为单一
入口脚本，减少复现步骤。
