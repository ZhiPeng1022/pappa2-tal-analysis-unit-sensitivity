# PAPPA2 全量数据再分析摘要 v1

日期：2026-09-30

## 方法与输入

- 完整 TAL：106,851 个细胞，223 位供者；
- 使用完整 h5ad 中已有 UMAP；
- Milo 使用 miloR 2.8.1，固定种子 20260927；
- 归一化方法：logMS；
- AKI 比较剔除 1 位没有邻域贡献的供者 163-5；
- 全量 neighborhood membership 已导出；
- 全量供者内秩相关使用 1,000 次置换。

## 全量 Milo 汇总

| comparison | subclass | total_annotated_nhoods | significant_nhoods | positive_nhoods | negative_nhoods | median_logFC | median_subclass_fraction | median_nhood_purity |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| AKI_vs_Reference | C-TAL-A | 3077 | 1747 | 151 | 1596 | -1.074 | 0.9722 | 0.9722 |
| AKI_vs_Reference | C-TAL-B | 2158 | 1002 | 158 | 844 | -0.8777 | 0.9722 | 0.9722 |
| AKI_vs_Reference | aTAL2 | 1223 | 842 | 822 | 20 | 1.364 | 0.6857 | 0.6857 |
| AKI_vs_Reference | frTAL | 638 | 246 | 168 | 78 | 0.3209 | 0.9487 | 0.9487 |
| CKD_vs_Reference | C-TAL-A | 3735 | 1172 | 290 | 882 | -0.4245 | 0.9706 | 0.9706 |
| CKD_vs_Reference | C-TAL-B | 2672 | 342 | 102 | 240 | -0.2213 | 0.973 | 0.973 |
| CKD_vs_Reference | aTAL2 | 782 | 282 | 34 | 248 | -0.4056 | 0.6216 | 0.6216 |
| CKD_vs_Reference | frTAL | 1268 | 982 | 956 | 26 | 1.138 | 0.971 | 0.971 |
| all_disease | C-TAL-A | 4012 | 1797 | 333 | 1464 | -0.4963 | 0.9706 | 0.9706 |
| all_disease | C-TAL-B | 2805 | 715 | 202 | 513 | -0.2942 | 0.9714 | 0.9714 |
| all_disease | aTAL2 | 1465 | 632 | 548 | 84 | 0.4371 | 0.6774 | 0.6774 |
| all_disease | frTAL | 1332 | 1023 | 996 | 27 | 0.9508 | 0.9706 | 0.9706 |

## PAPPA2 四层结果

| comparison | subclass | full_pb_mean_diff | full_pb_FDR | nhood_rho | nhood_direction | within_donor_rho_reference | within_donor_rho_disease | within_donor_disease_minus_reference | within_donor_difference_FDR | direction_conflict |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| AKI_vs_Reference | C-TAL-A | 0.8001 | 7.524e-12 | 0.5261 | 1 | 0.143 | 0.1677 | 0.02472 | 0.5278 | False |
| AKI_vs_Reference | C-TAL-B | 0.3887 | 3.741e-09 | 0.4722 | 1 | 0.05966 | 0.1057 | 0.04607 | 0.3885 | False |
| AKI_vs_Reference | aTAL2 | 0.2648 | 0.0001119 | 0.02879 | 1 | 0.1114 | -0.1742 | -0.2857 | 0.008563 | True |
| AKI_vs_Reference | frTAL | 0.139 | 0.08009 | 0.3324 | 1 | 0.1964 | 0.07891 | -0.1175 | 0.1456 | True |
| CKD_vs_Reference | C-TAL-A | 0.9387 | 6.186e-16 | 0.6885 | 1 | 0.1996 | 0.2066 | 0.007038 | 0.8203 | False |
| CKD_vs_Reference | C-TAL-B | 0.5181 | 2.558e-11 | 0.5354 | 1 | 0.06668 | 0.1596 | 0.09295 | 0.008563 | False |
| CKD_vs_Reference | aTAL2 | 0.4837 | 3.136e-09 | 0.4921 | 1 | 0.25 | 0.1613 | -0.08865 | 0.224 | True |
| CKD_vs_Reference | frTAL | 0.4764 | 9.869e-12 | 0.6077 | 1 | 0.2682 | 0.2125 | -0.0557 | 0.224 | True |
| all_disease | C-TAL-A | 0.894 | 1.595e-19 | 0.729 | 1 | 0.2005 | 0.2117 | 0.01116 | 0.6767 | False |
| all_disease | C-TAL-B | 0.4768 | 4.914e-14 | 0.6021 | 1 | 0.05277 | 0.1534 | 0.1006 | 0.01332 | False |
| all_disease | aTAL2 | 0.4114 | 9.219e-10 | 0.1453 | 1 | 0.1825 | -0.02375 | -0.2063 | 0.008563 | True |
| all_disease | frTAL | 0.3704 | 2.586e-09 | 0.5251 | 1 | 0.2707 | 0.181 | -0.08969 | 0.01798 | True |

## PAPPA2 完整邻域稳健性

| comparison | subclass | rho_full | lodo_same_sign_fraction | bootstrap_ci_lower | bootstrap_ci_upper | bootstrap_fraction_positive | partial_rho_controlling_purity |
| --- | --- | --- | --- | --- | --- | --- | --- |
| all_disease | C-TAL-A | 0.729 | 1 | 0.6651 | 0.7134 | 1 | 0.73 |
| all_disease | C-TAL-B | 0.6021 | 1 | 0.4748 | 0.5652 | 1 | 0.6024 |
| all_disease | aTAL2 | 0.1453 | 1 | 0.004537 | 0.3415 | 0.978 | 0.153 |
| all_disease | frTAL | 0.5251 | 1 | 0.4576 | 0.5347 | 1 | 0.5498 |
| AKI_vs_Reference | C-TAL-A | 0.5261 | 1 | 0.4576 | 0.5336 | 1 | 0.5603 |
| AKI_vs_Reference | C-TAL-B | 0.4722 | 1 | 0.318 | 0.4659 | 1 | 0.4778 |
| AKI_vs_Reference | aTAL2 | 0.02879 | 0.9917 | -0.08628 | 0.2551 | 0.632 | 0.04266 |
| AKI_vs_Reference | frTAL | 0.3324 | 1 | 0.2476 | 0.4146 | 1 | 0.4894 |
| CKD_vs_Reference | C-TAL-A | 0.6885 | 1 | 0.6138 | 0.6804 | 1 | 0.6779 |
| CKD_vs_Reference | C-TAL-B | 0.5354 | 1 | 0.4093 | 0.5104 | 1 | 0.5341 |
| CKD_vs_Reference | aTAL2 | 0.4921 | 1 | 0.3412 | 0.6076 | 1 | 0.5039 |
| CKD_vs_Reference | frTAL | 0.6077 | 1 | 0.529 | 0.6058 | 1 | 0.6092 |

## PAPPA2 方向汇总

- 供者伪批量正向：12/12；
- 供者伪批量 FDR < 0.05：11/12；
- 邻域 bootstrap 区间完全为正：11/12。

## 关键解释

1. 完整数据中 PAPPA2 的供者伪批量和邻域表达关联仍为正向；
2. aTAL2 的供者内秩相关在 all_disease 和 AKI 中相对负向，并在疾病组与 Reference 组差异检验中通过 FDR；
3. CKD 的 aTAL2 方向差为负，但差异检验未通过 FDR；
4. 因此应写成“跨分析单元冲突在完整数据中仍存在，但显著程度依疾病比较而异”，不能写成三个比较全部显著；
5. 全量 Milo 与旧子集 Milo 不能直接按显著数量比较，应分别报告输入细胞数、供者数和效应量。
