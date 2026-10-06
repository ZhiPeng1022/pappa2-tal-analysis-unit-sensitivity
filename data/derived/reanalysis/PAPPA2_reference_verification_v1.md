# PAPPA2 论文参考文献核对报告 v1

日期：2026-09-29

## 一、核对方法

1. 用 Europe PMC REST API 按 DOI 检索目标文献；
2. 核对标题、期刊、年份、卷、期、页码和 PMID；
3. 用 CrossRef DOI 解析交叉核对；
4. 网络资源类引用核对网址和访问日期。

## 二、核对结果

共 16 条参考文献：

- 15 条期刊论文，全部通过 DOI 和 PMID 核对；
- 1 条网络资源（Nephroseq v5），按数据库网址引用；
- 未使用任何无法解析的 DOI 或 PMID。

## 三、逐条状态

| 编号 | 第一作者或资源 | 期刊 | 年份 | DOI | PMID | 状态 |
|---:|---|---|---:|---|---|---|
| 1 | Bikbov B | Lancet | 2020 | 10.1016/S0140-6736(20)30045-3 | 32061315 | 已核对 |
| 2 | Hoste EAJ | Nat Rev Nephrol | 2018 | 10.1038/s41581-018-0052-0 | 30135570 | 已核对 |
| 3 | Ferenbach DA | Nat Rev Nephrol | 2015 | 10.1038/nrneph.2015.3 | 25643664 | 已核对 |
| 4 | Bonventre JV | J Clin Invest | 2011 | 10.1172/JCI45161 | 22045571 | 已核对 |
| 5 | Lake BB | Nature | 2023 | 10.1038/s41586-023-05769-3 | 37468583 | 已核对 |
| 6 | Knepper MA | J Am Soc Nephrol | 1999 | 10.1681/ASN.V103628 | 10073614 | 已核对 |
| 7 | Brezis M | J Clin Invest | 1984 | 10.1172/JCI111189 | 6690477 | 已核对 |
| 8 | Oxvig C | J Cell Commun Signal | 2015 | 10.1007/s12079-015-0259-9 | 25617049 | 已核对 |
| 9 | Dauber A | EMBO Mol Med | 2016 | 10.15252/emmm.201506106 | 26902202 | 已核对 |
| 10 | Conover CA | Endocrinology | 2011 | 10.1210/en.2011-0036 | 21586553 | 已核对 |
| 11 | Asghari M | Sci Adv | 2025 | 10.1126/sciadv.adv8918 | 40815665 | 已核对 |
| 12 | Kennedy C | Kidney360 | 2024 | 10.34067/KID.0000000602 | 39450948 | 已核对 |
| 13 | Nephroseq v5 | 网络资源 | - | - | - | 网址核对 |
| 14 | Dann E | Nat Biotechnol | 2022 | 10.1038/s41587-021-01033-z | 34594043 | 已核对 |
| 15 | Skinnider MA | Nat Biotechnol | 2021 | 10.1038/s41587-020-0605-1 | 32690972 | 已核对 |
| 16 | Cippa PE | JCI Insight | 2018 | 10.1172/jci.insight.123151 | 30429361 | 已核对 |

## 四、需要作者最终确认的事项

1. 目标期刊确定后，把 Vancouver 编号格式改成该期刊要求的格式；
2. 如果期刊要求，补充 PMCID 或 URL；
3. 核对作者单位、基金和致谢信息；
4. 核对 Nephroseq 的访问日期是否需要改为最终投稿日期；
5. 确认第 11 条 Asghari 等论文是否作为 GSE183456 空间资源的主要
   引用，或改用 KPMP 主图谱论文第 5 条；
6. 确认第 12 条 Kennedy 等论文是否为 GSE137570 的正式来源引用。

## 五、对应稿件

```text
outputs\PAPPA2_manuscript_draft_v2.md
```

## 六、复现脚本

```text
work\verify_refs.py
```
