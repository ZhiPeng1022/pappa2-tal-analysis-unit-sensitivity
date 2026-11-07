"""Audit sample metadata for all kidney-related GEO bulk candidates."""

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


MAX_WORKERS = 4


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
    candidates = pd.read_csv(
        OUT / "geo_bulk_kidney_clinical_candidates.csv",
        dtype={"accession": str},
    )
    kidney_mask = (
        candidates["title"].str.contains(
            r"kidney|renal|nephro|allograft",
            case=False,
            regex=True,
            na=False,
        )
        | candidates["summary"].str.contains(
            r"kidney|renal|nephro|allograft",
            case=False,
            regex=True,
            na=False,
        )
    )
    candidates = candidates[
        kidney_mask & candidates["n_samples"].ge(20)
    ].copy()
    candidates = candidates.sort_values(
        ["candidate_score", "n_samples"],
        ascending=[False, False],
    )
    metadata = candidates.set_index("accession")
    print(
        f"START candidates={len(candidates)} workers={MAX_WORKERS}",
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
                summary.update(
                    {
                        "candidate_score": metadata.loc[
                            accession,
                            "candidate_score",
                        ],
                        "query_hits": metadata.loc[
                            accession,
                            "query_hits",
                        ],
                    }
                )
                summary_rows.append(summary)
                characteristic_rows.extend(characteristics)
                sample_rows.extend(samples)
                fetch_status[accession] = status
            except Exception as error:
                fetch_status[accession] = f"ERROR: {error}"
                summary_rows.append(
                    {
                        "accession": accession,
                        "title": metadata.loc[accession, "title"],
                        "n_samples": int(
                            metadata.loc[accession, "n_samples"]
                        ),
                        "candidate_score": metadata.loc[
                            accession,
                            "candidate_score",
                        ],
                        "query_hits": metadata.loc[
                            accession,
                            "query_hits",
                        ],
                        "error": str(error),
                    }
                )
            if completed % 10 == 0 or completed == len(futures):
                print(
                    f"PROGRESS {completed}/{len(futures)}",
                    flush=True,
                )
            time.sleep(0.1)

    summary_frame = pd.DataFrame(summary_rows)
    summary_frame["clinical_priority_score"] = (
        summary_frame["gfr_numeric"].fillna(False).astype(int) * 10
        + summary_frame["fibrosis_keyword"].fillna(False).astype(int) * 4
        + summary_frame["progression_keyword"].fillna(False).astype(int) * 4
        + summary_frame["creatinine_keyword"].fillna(False).astype(int) * 3
        + summary_frame["age_key"].fillna(False).astype(int) * 2
        + summary_frame["sex_key"].fillna(False).astype(int) * 2
        + summary_frame["disease_key"].fillna(False).astype(int)
    )
    summary_frame = summary_frame.sort_values(
        ["clinical_priority_score", "n_samples"],
        ascending=[False, False],
    )
    summary_frame.to_csv(
        OUT / "geo_all_kidney_sample_audit_summary.csv",
        index=False,
    )
    characteristic_frame = pd.DataFrame(characteristic_rows)
    characteristic_frame.to_csv(
        OUT / "geo_all_kidney_characteristics_summary.csv",
        index=False,
    )
    characteristic_frame[
        characteristic_frame["characteristic"].map(
            is_meaningful_characteristic
        )
    ].to_csv(
        OUT / "geo_all_kidney_clinical_characteristics_curated.csv",
        index=False,
    )
    pd.DataFrame(sample_rows).to_csv(
        OUT / "geo_all_kidney_characteristics_long.csv",
        index=False,
    )
    pd.DataFrame(
        [
            {"accession": accession, "status": status}
            for accession, status in sorted(fetch_status.items())
        ]
    ).to_csv(
        OUT / "geo_all_kidney_fetch_status.csv",
        index=False,
    )
    print("TOP CANDIDATES", flush=True)
    print(
        summary_frame[
            [
                "accession",
                "n_samples",
                "clinical_priority_score",
                "gfr_numeric",
                "fibrosis_keyword",
                "progression_keyword",
                "creatinine_keyword",
                "age_key",
                "sex_key",
                "title",
            ]
        ]
        .head(40)
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
