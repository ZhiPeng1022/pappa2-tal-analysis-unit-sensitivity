"""Final broad scan of all human kidney bulk GEO candidates."""

from __future__ import annotations

import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd

from audit_bulk_clinical_headers import (
    OUT,
    download_sample_metadata,
    is_meaningful_characteristic,
    parse_sample_records,
    summarize_dataset,
)


MAX_WORKERS = 6


def audit_one(accession: str, title: str) -> tuple[str, dict, list, list, str]:
    path, status = download_sample_metadata(accession)
    records = parse_sample_records(path.read_text(encoding="utf-8"))
    summary, characteristics, samples = summarize_dataset(
        accession,
        title,
        records,
    )
    return accession, summary, characteristics, samples, status


def main() -> None:
    raw = pd.read_csv(
        OUT / "geo_bulk_kidney_search_raw.csv",
        dtype={"accession": str},
    )
    kidney_mask = (
        raw["title"].str.contains(
            r"kidney|renal|nephro|allograft",
            case=False,
            regex=True,
            na=False,
        )
        | raw["summary"].str.contains(
            r"kidney|renal|nephro|allograft",
            case=False,
            regex=True,
            na=False,
        )
    )
    candidates = raw[
        raw["taxon"].eq("Homo sapiens")
        & raw["n_samples"].ge(10)
        & kidney_mask
    ].drop_duplicates("accession")
    candidates = candidates.sort_values(
        ["candidate_score", "n_samples"],
        ascending=[False, False],
    )
    print(
        f"START broad_candidates={len(candidates)} workers={MAX_WORKERS}",
        flush=True,
    )

    summary_rows = []
    characteristic_rows = []
    sample_rows = []
    fetch_status = {}
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {
            executor.submit(
                audit_one,
                row["accession"],
                row["title"],
            ): row["accession"]
            for _, row in candidates.iterrows()
        }
        completed = 0
        for future in as_completed(futures):
            accession = futures[future]
            completed += 1
            try:
                (
                    accession,
                    summary,
                    characteristics,
                    samples,
                    status,
                ) = future.result()
                summary_rows.append(summary)
                characteristic_rows.extend(characteristics)
                sample_rows.extend(samples)
                fetch_status[accession] = status
            except Exception as error:
                fetch_status[accession] = f"ERROR: {error}"
                summary_rows.append(
                    {
                        "accession": accession,
                        "title": "",
                        "n_samples": 0,
                        "error": str(error),
                    }
                )
            if completed % 20 == 0 or completed == len(futures):
                print(
                    f"PROGRESS {completed}/{len(futures)}",
                    flush=True,
                )
            time.sleep(0.05)

    summary = pd.DataFrame(summary_rows)
    characteristics = pd.DataFrame(characteristic_rows)
    samples = pd.DataFrame(sample_rows)
    summary["clinical_priority_score"] = (
        summary["gfr_numeric"].fillna(False).astype(int) * 10
        + summary["fibrosis_keyword"].fillna(False).astype(int) * 4
        + summary["progression_keyword"].fillna(False).astype(int) * 4
        + summary["creatinine_keyword"].fillna(False).astype(int) * 3
        + summary["age_key"].fillna(False).astype(int) * 2
        + summary["sex_key"].fillna(False).astype(int) * 2
        + summary["disease_key"].fillna(False).astype(int)
    )
    summary = summary.sort_values(
        ["clinical_priority_score", "n_samples"],
        ascending=[False, False],
    )
    summary.to_csv(
        OUT / "geo_broad_bulk_kidney_audit_summary.csv",
        index=False,
    )
    characteristics.to_csv(
        OUT / "geo_broad_bulk_kidney_characteristics.csv",
        index=False,
    )
    characteristics[
        characteristics["characteristic"].map(
            is_meaningful_characteristic
        )
    ].to_csv(
        OUT / "geo_broad_bulk_kidney_clinical_characteristics_curated.csv",
        index=False,
    )
    samples.to_csv(
        OUT / "geo_broad_bulk_kidney_sample_characteristics_long.csv",
        index=False,
    )
    pd.DataFrame(
        [
            {"accession": accession, "status": status}
            for accession, status in sorted(fetch_status.items())
        ]
    ).to_csv(
        OUT / "geo_broad_bulk_kidney_fetch_status.csv",
        index=False,
    )
    top = summary[
        summary["gfr_numeric"].fillna(False)
    ].head(30)
    print("NUMERIC GFR CANDIDATES", flush=True)
    print(
        top[
            [
                "accession",
                "n_samples",
                "age_key",
                "sex_key",
                "fibrosis_keyword",
                "progression_keyword",
                "title",
            ]
        ].to_string(index=False),
        flush=True,
    )


if __name__ == "__main__":
    main()
