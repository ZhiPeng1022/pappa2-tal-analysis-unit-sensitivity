from pathlib import Path
import csv


PROJECT_ROOT = Path(
    r"C:\Users\Elsa\Documents\Codex\2026-09-25\9-2"
)
OUTPUT_ROOT = PROJECT_ROOT / "outputs"

MARKDOWN_OUTPUT = (
    OUTPUT_ROOT / "PAPPA2_manuscript_evidence_matrix_v1.md"
)
CSV_OUTPUT = (
    OUTPUT_ROOT / "PAPPA2_manuscript_evidence_matrix_v1.csv"
)

CLAIMS = [
    {
        "claim_id": "E01",
        "paper_position": "引言、方法总览",
        "claim": (
            "供者内邻域关联框架可以整合 TAL 亚类、疾病相关邻域和"
            "外部 GFR 证据。"
        ),
        "evidence_level": "方法整合",
        "status": "可作为研究框架",
        "key_evidence": (
            "整合 4,446 个 TAL 细胞、Milo、Augur、供者内置换、"
            "GSE183456、GSE137570 和 Nephroseq。"
        ),
        "supporting_files": (
            "outputs/PAPPA2_CAS2_frozen_analysis_protocol_v1.md; "
            "outputs/PAPPA2_CAS2_decisive_validation_summary.md"
        ),
        "allowed_wording": (
            "我们提出并应用了一个供者内、亚类分辨的整合分析框架。"
        ),
        "forbidden_wording": (
            "不能声称该框架已经优于所有现有邻域方法。"
        ),
        "reviewer_risk": (
            "如果投稿定位为方法学论文，审稿人会要求与其他方法比较。"
        ),
        "next_verification": (
            "若走方法学路线，补一个与标准 Milo 或简单供者混合模型的"
            "基准比较。"
        ),
    },
    {
        "claim_id": "E02",
        "paper_position": "结果 1",
        "claim": (
            "TAL 疾病相关邻域变化具有明显亚类异质性。"
        ),
        "evidence_level": "B",
        "status": "主结果，可写",
        "key_evidence": (
            "all_disease：C-TAL-A 62/128显著；C-TAL-B 32/131；"
            "aTAL2 58/68；frTAL 67/84。CKD 中 aTAL2 median logFC "
            "-1.12，而 frTAL +1.24。"
        ),
        "supporting_files": (
            "outputs/KPMP_milo_TAL_subclass_summary.csv; "
            "outputs/PAPPA2_TAL_milo_augur_summary.md"
        ),
        "allowed_wording": (
            "TAL 亚类的疾病相关邻域状态和方向存在明显差异。"
        ),
        "forbidden_wording": (
            "不能写四个 TAL 亚类整体扩增或方向一致。"
        ),
        "reviewer_risk": (
            "Milo 固定种子敏感性提示 AKI 邻域数不稳定。"
        ),
        "next_verification": (
            "正式结果继续使用原 Milo 文件，固定种子只作为敏感性。"
        ),
    },
    {
        "claim_id": "E03",
        "paper_position": "结果 1、讨论",
        "claim": (
            "aTAL2 与其余 TAL 亚类呈现候选状态解耦。"
        ),
        "evidence_level": "D",
        "status": "探索性结果",
        "key_evidence": (
            "供者内 PAPPA2 关联在 aTAL2 三个比较中均为负；"
            "all_disease FDR=0.0466，CKD 和 AKI 未通过 FDR。"
        ),
        "supporting_files": (
            "outputs/KPMP_TAL_seed20260927_within_donor_cell_rank.csv; "
            "outputs/KPMP_TAL_seed20260927_within_donor_cell_rank_"
            "threshold_sensitivity.csv"
        ),
        "allowed_wording": (
            "aTAL2 提示与其余亚类相反的邻域状态，值得进一步验证。"
        ),
        "forbidden_wording": (
            "不能写 aTAL2 在所有条件下稳定反向，或写成独立机制。"
        ),
        "reviewer_risk": (
            "CKD aTAL2 有效供者数和邻域数较少。"
        ),
        "next_verification": (
            "用供者层级模型或独立单细胞队列验证 aTAL2 反向。"
        ),
    },
    {
        "claim_id": "E04",
        "paper_position": "结果 2",
        "claim": (
            "PAPPA2 在四种 TAL 亚类的疾病状态中均升高。"
        ),
        "evidence_level": "B",
        "status": "主结果，可写",
        "key_evidence": (
            "12 个基因 x 亚类 x 比较组合全部 FDR<0.05 且方向为正；"
            "例如 all_disease C-TAL-A +1.434，frTAL +0.954。"
        ),
        "supporting_files": (
            "outputs/KPMP_TAL_focus_gene_pseudobulk.csv"
        ),
        "allowed_wording": (
            "PAPPA2 是 TAL 亚类中最一致的疾病升高候选基因。"
        ),
        "forbidden_wording": (
            "不能由表达升高推断 PAPPA2 导致疾病。"
        ),
        "reviewer_risk": (
            "伪批量结果仍是关联分析，并受供者组成影响。"
        ),
        "next_verification": (
            "报告供者数、效应方向和 FDR，不把表达变化写成功能效应。"
        ),
    },
    {
        "claim_id": "E05",
        "paper_position": "结果 2",
        "claim": (
            "PAPPA2 与疾病相关 Milo 邻域暴露存在供者内关联，"
            "C-TAL-A、C-TAL-B 和 frTAL 最稳定。"
        ),
        "evidence_level": "B/C",
        "status": "主结果，可写",
        "key_evidence": (
            "12 个组中 11 个 Spearman 为正，中位数 0.711；"
            "10/12 的供者重采样区间不跨 0；9/9 非 aTAL2 "
            "在 3、5、10 细胞阈值下均为正。"
        ),
        "supporting_files": (
            "outputs/PAPPA2_TAL_within_donor_rank_summary.md; "
            "outputs/KPMP_TAL_seed20260927_within_donor_cell_rank.csv"
        ),
        "allowed_wording": (
            "PAPPA2 表达与疾病相关邻域状态在同供者内部呈稳定关联。"
        ),
        "forbidden_wording": (
            "不能写邻域是独立样本、因果关系或直接细胞通讯。"
        ),
        "reviewer_risk": (
            "邻域高度重叠，标准 p 值不能按独立样本解释。"
        ),
        "next_verification": (
            "保留置换检验和阈值敏感性，强调方向稳定性。"
        ),
    },
    {
        "claim_id": "E06",
        "paper_position": "结果 3",
        "claim": (
            "KPMP 供者中 PAPPA2 阳性比例与较低 eGFR 独立相关。"
        ),
        "evidence_level": "B",
        "status": "内部临床主结果",
        "key_evidence": (
            "主模型每增加 10 个百分点 OR=0.819（95% CI "
            "0.682-0.985，p=0.0335）；疾病合并 OR=0.789；"
            "CKD OR=0.733。AKI 无关联。"
        ),
        "supporting_files": (
            "outputs/KPMP_PAPPA2_TAL_egfr_models.csv; "
            "outputs/KPMP_PAPPA2_TAL_egfr_gene_specificity.csv"
        ),
        "allowed_wording": (
            "KPMP 内部模型支持 PAPPA2 与较低 eGFR 相关。"
        ),
        "forbidden_wording": (
            "不能写成独立临床验证，或连续 eGFR 的精确剂量反应。"
        ),
        "reviewer_risk": (
            "eGFR 是分箱变量，AKI 样本量不足，模型包含多重比较。"
        ),
        "next_verification": (
            "在讨论中报告 bin 限制，并明确 CKD 是主要信号来源。"
        ),
    },
    {
        "claim_id": "E07",
        "paper_position": "结果 3",
        "claim": (
            "PAPPA2 的 GFR 和邻域关联方向具有候选基因特异性。"
        ),
        "evidence_level": "B/C",
        "status": "支持性主结果",
        "key_evidence": (
            "KPMP 供者内秩相关中 PAPPA2 有 7 个 FDR<0.05，"
            "CXCL12 2 个，其他基因更少；Nephroseq 中 PAPPA2 "
            "tubulointerstitium 负向结果最多（7 个）。"
        ),
        "supporting_files": (
            "outputs/KPMP_TAL_seed20260927_within_donor_cell_rank.csv; "
            "outputs/nephroseq_GFR_gene_summary.csv"
        ),
        "allowed_wording": (
            "在本次五个候选基因中，PAPPA2 的方向一致性和"
            "tubulointerstitium 负向结果最多。"
        ),
        "forbidden_wording": (
            "不能写 PAPPA2 是唯一相关基因或效应量最大基因。"
        ),
        "reviewer_risk": (
            "HES1 和 EPHB2 也有显著结果，CXCL12 方向相反。"
        ),
        "next_verification": (
            "使用随机匹配基因或负对照进一步验证特异性。"
        ),
    },
    {
        "claim_id": "E08",
        "paper_position": "结果 4",
        "claim": (
            "GSE183456 跨队列空间复现 PAPPA2 分支签名的共现方向。"
        ),
        "evidence_level": "C",
        "status": "支持性结果",
        "key_evidence": (
            "23 个切片中 20 个共现 OR>1；DKD 11/11；供者去重后"
            "10/10；但 DKD vs Reference p=0.078，Reference "
            "本身也有 5/6 共现。"
        ),
        "supporting_files": (
            "outputs/PAPPA2_23sample_spatial_validation_summary.md; "
            "outputs/KPMP_Visium_PAPPA2_branch_signature_"
            "sample_associations.csv"
        ),
        "allowed_wording": (
            "空间方向在跨队列 Visium 数据中保留。"
        ),
        "forbidden_wording": (
            "不能写疾病特异验证或完全独立验证。"
        ),
        "reviewer_risk": (
            "GSE183456 已被前期探索使用，且 Reference 也有共现。"
        ),
        "next_verification": (
            "明确写 external cross-cohort spatial replication。"
        ),
    },
    {
        "claim_id": "E09",
        "paper_position": "补充结果",
        "claim": (
            "GSE137570 Cohort 1 的 PAPPA2 与较低 eGFR 方向一致，"
            "但样本量不足。"
        ),
        "evidence_level": "D",
        "status": "有限支持",
        "key_evidence": (
            "24 个供者 rho=-0.365，p=0.079；调整后 p=0.351；"
            "EPHB2 未调整 rho=-0.646，p=0.0006。"
        ),
        "supporting_files": (
            "outputs/GSE137570_PAPPA2_eGFR_validation.csv; "
            "outputs/PAPPA2_external_cohort_search_v1.md"
        ),
        "allowed_wording": (
            "方向一致但 underpowered 的公开外部队列支持。"
        ),
        "forbidden_wording": (
            "不能写显著性验证或独立临床验证。"
        ),
        "reviewer_risk": (
            "只有 24 例公开 eGFR，调整模型置信区间很宽。"
        ),
        "next_verification": (
            "向作者申请 Cohort 2 eGFR，纳入 41 例分析。"
        ),
    },
    {
        "claim_id": "E10",
        "paper_position": "结果 5",
        "claim": (
            "Nephroseq 多数据集显示 PAPPA2 与较低 GFR 的负向关联。"
        ),
        "evidence_level": "A-",
        "status": "外部临床支持主结果",
        "key_evidence": (
            "50 个 GFR 分析中 10 个显著，9 个负向，7 个来自"
            "tubulointerstitium；主要结果包括 Sampson TubInt "
            "r=-0.558（n=49）、Ju CKD TubInt 2 r=-0.530（n=42）、"
            "ERCB TubInt r=-0.712（n=17）。"
        ),
        "supporting_files": (
            "outputs/PAPPA2_Nephroseq_GFR_validation_summary.md; "
            "outputs/nephroseq_GFR_significant_analyses.csv"
        ),
        "allowed_wording": (
            "Nephroseq 的多个公开肾脏数据集中，PAPPA2 与较低 GFR "
            "呈一致负向关联。"
        ),
        "forbidden_wording": (
            "不能把 10 个分析写成 10 个独立队列，也不能写成因果。"
        ),
        "reviewer_risk": (
            "Nephroseq 是预计算分析，GFR 公式不同，部分数据集"
            "样本较小或同源。"
        ),
        "next_verification": (
            "核对原始数据集和引用，按数据集而非分析行数报告。"
        ),
    },
    {
        "claim_id": "E11",
        "paper_position": "结果 5",
        "claim": (
            "在测试的五个候选基因中，PAPPA2 的 tubulointerstitium "
            "负向 GFR 结果最多。"
        ),
        "evidence_level": "A-",
        "status": "支持性结果",
        "key_evidence": (
            "Nephroseq tubulointerstitium 负向显著分析：PAPPA2 7、"
            "HES1 6、EPHB2 3、NOTCH3 2、CXCL12 1。"
        ),
        "supporting_files": (
            "outputs/nephroseq_GFR_gene_summary.csv; "
            "outputs/nephroseq_GFR_significant_analyses.csv"
        ),
        "allowed_wording": (
            "PAPPA2 在本次对照基因中显示出最多的 tubulointerstitium "
            "负向 GFR 证据。"
        ),
        "forbidden_wording": (
            "不能写 PAPPA2 是唯一或经验证的唯一 GFR 相关基因。"
        ),
        "reviewer_risk": (
            "HES1 的结果数量接近，且有 2 个正向分析。"
        ),
        "next_verification": (
            "按数据集去重后重新统计，并增加随机基因对照。"
        ),
    },
    {
        "claim_id": "E12",
        "paper_position": "讨论",
        "claim": (
            "JAG1-NOTCH3-HES1 在 TAL 亚类中没有被统一复现。"
        ),
        "evidence_level": "E",
        "status": "阴性结果",
        "key_evidence": (
            "JAG1、HES4 未显著，NOTCH3 只在 CKD frTAL 中显著，"
            "HES1 方向不一致；Nephroseq 中 HES1 也正负混合。"
        ),
        "supporting_files": (
            "outputs/KPMP_TAL_focus_gene_pseudobulk.csv; "
            "outputs/PAPPA2_SCENIC_summary.md"
        ),
        "allowed_wording": (
            "Notch 相关结果未形成统一 TAL 亚类模式。"
        ),
        "forbidden_wording": (
            "不能写 JAG1-NOTCH3-HES1 已成立或已复现。"
        ),
        "reviewer_risk": (
            "把零散显著项组合成机制链会被审稿人质疑。"
        ),
        "next_verification": (
            "保持为阴性或探索性结果，不强行拼接通路。"
        ),
    },
    {
        "claim_id": "E13",
        "paper_position": "讨论",
        "claim": (
            "EPHB2-STAT6 轴没有成立。"
        ),
        "evidence_level": "E",
        "status": "阴性结果",
        "key_evidence": (
            "STAT6 在供者伪批量中无显著变化；EPHB2 在不同亚类"
            "方向不一致，Nephroseq 中虽有负向 GFR 结果，"
            "但没有同步 STAT6 证据。"
        ),
        "supporting_files": (
            "outputs/KPMP_TAL_focus_gene_pseudobulk.csv; "
            "outputs/PAPPA2_directional_communication_summary.md"
        ),
        "allowed_wording": (
            "EPHB2 与 STAT6 没有形成一致 TAL 亚类证据链。"
        ),
        "forbidden_wording": (
            "不能写 EPHB2-STAT6 轴已成立。"
        ),
        "reviewer_risk": (
            "只有配体或受体一侧变化不能证明通路激活。"
        ),
        "next_verification": (
            "不继续扩展；若需要，必须依赖功能实验。"
        ),
    },
    {
        "claim_id": "E14",
        "paper_position": "讨论、限制",
        "claim": (
            "当前结果支持关联，不支持 PAPPA2 的因果作用。"
        ),
        "evidence_level": "E",
        "status": "限制边界",
        "key_evidence": (
            "所有结果来自观察性组学数据和预计算关联；"
            "没有 PAPPA2 扰动、蛋白验证或功能实验。"
        ),
        "supporting_files": (
            "outputs/PAPPA2_CAS2_decisive_validation_summary.md"
        ),
        "allowed_wording": (
            "PAPPA2 与疾病相关 TAL 状态和较低 GFR 相关。"
        ),
        "forbidden_wording": (
            "不能写导致、驱动、介导或直接调控。"
        ),
        "reviewer_risk": (
            "因果措辞会显著降低可信度。"
        ),
        "next_verification": (
            "论文中统一使用 associated with、linked to。"
        ),
    },
    {
        "claim_id": "E15",
        "paper_position": "讨论、限制",
        "claim": (
            "完全独立的外部临床验证尚未完成。"
        ),
        "evidence_level": "E",
        "status": "剩余限制",
        "key_evidence": (
            "GSE183456已被前期使用；GSE137570只有24例公开eGFR；"
            "Nephroseq是预计算外部支持，不等同于统一模型验证。"
        ),
        "supporting_files": (
            "outputs/PAPPA2_external_cohort_search_v1.md; "
            "outputs/PAPPA2_CAS2_decisive_validation_summary.md"
        ),
        "allowed_wording": (
            "结果得到跨队列空间方向和多数据集临床关联支持。"
        ),
        "forbidden_wording": (
            "不能写完全独立外部验证。"
        ),
        "reviewer_risk": (
            "这是二区审稿时最可能被追问的限制。"
        ),
        "next_verification": (
            "申请GSE137570 Cohort 2或432例微分离小管队列。"
        ),
    },
    {
        "claim_id": "E16",
        "paper_position": "数据与方法",
        "claim": (
            "主要分析流程和关键结果具有可复现文件。"
        ),
        "evidence_level": "A",
        "status": "可写",
        "key_evidence": (
            "分析脚本、冻结协议、置换结果、空间哈希核验和"
            "Nephroseq查询脚本均已保存。"
        ),
        "supporting_files": (
            "outputs/PAPPA2_CAS2_frozen_analysis_protocol_v1.md; "
            "work/query_nephroseq_gfr_genes.py; "
            "work/sensitivity_tal_nhood_expression.py"
        ),
        "allowed_wording": (
            "研究和结果文件均可复现。"
        ),
        "forbidden_wording": (
            "不能省略GSE183456前期使用等数据角色说明。"
        ),
        "reviewer_risk": (
            "数据版本、阈值和脚本哈希需要进入补充材料。"
        ),
        "next_verification": (
            "整理代码清单、版本号和随机种子表。"
        ),
    },
    {
        "claim_id": "E17",
        "paper_position": "摘要、结论",
        "claim": (
            "供者内邻域分析识别出与肾功能相关的 PAPPA2 相关 TAL "
            "状态。"
        ),
        "evidence_level": "条件性综合",
        "status": "可作条件性主主张",
        "key_evidence": (
            "KPMP供者内关联、KPMP内部eGFR、GSE183456空间方向、"
            "GSE137570有限支持、Nephroseq多数据集GFR共同支持。"
        ),
        "supporting_files": (
            "outputs/PAPPA2_Nephroseq_GFR_validation_summary.md; "
            "outputs/PAPPA2_TAL_within_donor_rank_summary.md; "
            "outputs/PAPPA2_CAS2_decisive_validation_summary.md"
        ),
        "allowed_wording": (
            "PAPPA2 相关 TAL 邻域状态与肾功能降低相关。"
        ),
        "forbidden_wording": (
            "不能写 PAPPA2 是疾病驱动因子或候选治疗靶点。"
        ),
        "reviewer_risk": (
            "主主张仍需依赖方法整合和多层关联证据。"
        ),
        "next_verification": (
            "所有图表围绕这一关联主线，不扩展到因果机制。"
        ),
    },
]


def write_markdown():
    lines = [
        "# PAPPA2 论文证据矩阵 v1",
        "",
        "日期：2026-09-28",
        "",
        "## 一、证据等级",
        "",
        "- `A`：多数据集外部队列或公开数据库支持，"
        "但仍需核对数据集重叠。",
        "- `A-`：多数据集外部关联支持，但为预计算结果，"
        "不能统一调整协变量。",
        "- `B`：单一主队列中的调整后人类数据结果。",
        "- `C`：不同队列或不同分析层之间方向复现。",
        "- `D`：探索性或敏感性分析，不能单独承担论文主结论。",
        "- `E`：机制或因果证据不足。",
        "",
        "## 二、核心主张",
        "",
        "| ID | 论文位置 | 主张 | 等级 | 状态 | 关键证据 |",
        "|---|---|---|---|---|---|",
    ]
    for claim in CLAIMS:
        lines.append(
            "| {claim_id} | {paper_position} | {claim} | "
            "{evidence_level} | {status} | {key_evidence} |".format(
                **claim
            )
        )

    lines.extend(
        [
            "",
            "## 三、允许与禁止表述",
            "",
            "| ID | 可写 | 不可写 | 审稿风险 | 下一步 |",
            "|---|---|---|---|---|",
        ]
    )
    for claim in CLAIMS:
        lines.append(
            "| {claim_id} | {allowed_wording} | "
            "{forbidden_wording} | {reviewer_risk} | "
            "{next_verification} |".format(**claim)
        )

    lines.extend(
        [
            "",
            "## 四、论文主结论建议",
            "",
            "优先使用：",
            "",
            "```text",
            "供者内邻域分析识别出一个与较低肾功能相关的",
            "PAPPA2 相关 TAL 状态。",
            "```",
            "",
            "不要把主结论写成：",
            "",
            "```text",
            "PAPPA2 驱动 TAL 损伤或肾病进展。",
            "```",
            "",
            "## 五、结构化文件",
            "",
            "```text",
            "outputs\\PAPPA2_manuscript_evidence_matrix_v1.csv",
            "```",
            "",
            "## 六、主要原始证据",
            "",
            "```text",
            "outputs\\KPMP_milo_TAL_subclass_summary.csv",
            "outputs\\KPMP_Augur_TAL_full_AUC_summary.csv",
            "outputs\\KPMP_TAL_focus_gene_pseudobulk.csv",
            "outputs\\KPMP_TAL_seed20260927_within_donor_cell_rank.csv",
            "outputs\\KPMP_PAPPA2_TAL_egfr_models.csv",
            "outputs\\PAPPA2_23sample_spatial_validation_summary.md",
            "outputs\\GSE137570_PAPPA2_eGFR_validation.csv",
            "outputs\\PAPPA2_Nephroseq_GFR_validation_summary.md",
            "outputs\\nephroseq_GFR_gene_summary.csv",
            "outputs\\nephroseq_GFR_significant_analyses.csv",
            "```",
        ]
    )
    MARKDOWN_OUTPUT.write_text(
        "\n".join(lines),
        encoding="utf-8",
    )


def write_csv():
    columns = [
        "claim_id",
        "paper_position",
        "claim",
        "evidence_level",
        "status",
        "key_evidence",
        "supporting_files",
        "allowed_wording",
        "forbidden_wording",
        "reviewer_risk",
        "next_verification",
    ]
    with CSV_OUTPUT.open(
        "w",
        encoding="utf-8-sig",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
        )
        writer.writeheader()
        writer.writerows(CLAIMS)


def main():
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    write_markdown()
    write_csv()
    print("Saved:", MARKDOWN_OUTPUT)
    print("Saved:", CSV_OUTPUT)
    print("Claims:", len(CLAIMS))


if __name__ == "__main__":
    main()
