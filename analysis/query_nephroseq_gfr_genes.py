import os
from pathlib import Path
import json
import time
from urllib.request import Request, urlopen

import pandas as pd


PROJECT_ROOT = Path(
    os.environ.get(
        "PAPPA2_PROJECT_ROOT",
        Path(__file__).resolve().parents[2],
    )
)
WORK_ROOT = Path(
    os.environ.get(
        "PAPPA2_WORK_ROOT",
        PROJECT_ROOT / "work",
    )
)
SUMMARY_OUTPUT = (
    WORK_ROOT / "nephroseq_GFR_gene_summary.csv"
)
DETAIL_OUTPUT = (
    WORK_ROOT / "nephroseq_GFR_significant_analyses.csv"
)

API_ROOT = "https://nephroseq.org/api"
GENES = [
    "PAPPA2",
    "CXCL12",
    "NOTCH3",
    "HES1",
    "EPHB2",
]
THRESHOLDS = {
    "pValue": 0.05,
    "rValue": 0.5,
    "foldChange": 1.5,
}
GFR_FILTER_ID = 77


def get_json(url):
    request = Request(
        url,
        headers={
            "Accept": "application/json",
            "User-Agent": "Codex nephroseq validation query",
        },
    )
    with urlopen(request, timeout=120) as response:
        return json.load(response)


def post_json(url, payload):
    request = Request(
        url,
        data=json.dumps(payload).encode("utf-8"),
        headers={
            "Accept": "application/json",
            "Content-Type": "application/json",
            "Origin": "https://nephroseq.org",
            "Referer": "https://nephroseq.org/",
            "User-Agent": "Codex nephroseq validation query",
        },
        method="POST",
    )
    with urlopen(request, timeout=120) as response:
        return json.load(response)


def resolve_gene(gene):
    request_id = int(time.time_ns() % 900) + 100
    result = get_json(
        f"{API_ROOT}/gene/{gene}/requestId/{request_id}"
    )
    genes = [
        item
        for item in result.get("genes", [])
        if item.get("symbol") == gene
    ]
    if len(genes) != 1:
        raise ValueError(
            f"Expected one Nephroseq gene for {gene}, got {genes}"
        )
    return int(genes[0]["id"])


def extract_analysis_rows(
    gene,
    analysis_response,
):
    rows = []
    for analysis in analysis_response.get(
        "analysisDisplayModels",
        [],
    ):
        dataset = analysis["datasetDisplayModel"]
        rows.append(
            {
                "gene": gene,
                "analysis_id": analysis["analysisId"],
                "analysis_name": analysis["analysisName"],
                "analysis_type": analysis["typeOfAnalysis"],
                "analysis_synopsis": analysis[
                    "analysisSynopsis"
                ],
                "p_value": analysis["pValue"],
                "r_value": analysis["rValue"],
                "fold_change": analysis["foldChange"],
                "reporter": analysis["reporterName"],
                "dataset_id": dataset["datasetId"],
                "dataset_name": dataset["name"],
                "organism": dataset["organism"],
                "tissue_types": dataset["tissueTypes"],
                "experiment_type": dataset[
                    "experimentType"
                ],
                "platforms": "; ".join(
                    dataset["platforms"]
                ),
                "dataset_title": dataset["title"],
                "pubmed_id": dataset["pubmedId"],
                "publication_date": dataset[
                    "publicationDate"
                ],
                "analysis_samples": analysis[
                    "analysisSampleCount"
                ],
                "dataset_samples": analysis[
                    "datasetSampleCount"
                ],
            }
        )
    return rows


def main():
    summary_rows = []
    detail_rows = []
    for gene in GENES:
        gene_id = resolve_gene(gene)
        summary = post_json(
            (
                f"{API_ROOT}/filters/geneSummary/"
                f"gene/{gene_id}"
            ),
            THRESHOLDS,
        )
        gfr_summary = next(
            (
                row
                for row in summary
                if row["filter"]["id"] == GFR_FILTER_ID
            ),
            None,
        )
        if gfr_summary is None:
            raise ValueError(
                f"No GFR summary returned for {gene}"
            )
        summary_rows.append(
            {
                "gene": gene,
                "nephroseq_gene_id": gene_id,
                "significant_gfr_analyses": (
                    gfr_summary["significantAnalyses"]
                ),
                "unique_gfr_analyses": (
                    gfr_summary["uniqueAnalyses"]
                ),
            }
        )

        request_id = int(
            (time.time_ns() + gene_id) % 900
        ) + 100
        details = post_json(
            (
                f"{API_ROOT}/analysis/filters/"
                f"{GFR_FILTER_ID}/gene/{gene_id}/"
                f"requestId/{request_id}"
            ),
            THRESHOLDS,
        )
        detail_rows.extend(
            extract_analysis_rows(
                gene,
                details,
            )
        )

    summary_table = pd.DataFrame(summary_rows)
    detail_table = pd.DataFrame(detail_rows)
    summary_table["significant_negative"] = (
        summary_table["gene"].map(
            detail_table.assign(
                negative=detail_table["r_value"] < 0
            )
            .groupby("gene")["negative"]
            .sum()
        )
    )
    summary_table["significant_positive"] = (
        summary_table["gene"].map(
            detail_table.assign(
                positive=detail_table["r_value"] > 0
            )
            .groupby("gene")["positive"]
            .sum()
        )
    )
    summary_table["significant_tubulointerstitium"] = (
        summary_table["gene"].map(
            detail_table.assign(
                tubulointerstitium=detail_table[
                    "tissue_types"
                ]
                .astype(str)
                .str.contains(
                    "Tubulointerstitium",
                    case=False,
                    na=False,
                )
            )
            .groupby("gene")["tubulointerstitium"]
            .sum()
        )
    )
    summary_table = summary_table.fillna(0)

    summary_table.to_csv(
        SUMMARY_OUTPUT,
        index=False,
        encoding="utf-8-sig",
    )
    detail_table.to_csv(
        DETAIL_OUTPUT,
        index=False,
        encoding="utf-8-sig",
    )
    print("Saved:", SUMMARY_OUTPUT)
    print("Saved:", DETAIL_OUTPUT)
    print("\nGFR summary:")
    print(summary_table.to_string(index=False))
    print("\nPAPPA2 significant analyses:")
    print(
        detail_table[
            detail_table["gene"] == "PAPPA2"
        ]
        .sort_values("p_value")[
            [
                "dataset_name",
                "tissue_types",
                "analysis_synopsis",
                "analysis_name",
                "p_value",
                "r_value",
                "analysis_samples",
            ]
        ]
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
