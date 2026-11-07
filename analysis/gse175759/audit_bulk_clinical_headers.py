"""Audit sample-level clinical annotations through the GEO text interface."""

from __future__ import annotations

import json
import re
import time
from pathlib import Path

import pandas as pd
import requests


ROOT = Path(r"C:\Users\Elsa\Documents\Codex\2026-10-15\ji")
OUT = ROOT / "work" / "bulk_clinical_cohort_audit"
METADATA_DIR = OUT / "sample_metadata"
GEO_TEXT = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"

SHORTLIST = [
    "GSE104954",
    "GSE99325",
    "GSE57387",
    "GSE22459",
    "GSE126805",
    "GSE222889",
    "GSE21374",
    "GSE25902",
    "GSE273763",
    "GSE96804",
    "GSE30122",
    "GSE66494",
    "GSE115857",
    "GSE116626",
    "GSE175759",
    "GSE127797",
    "GSE32591",
    "GSE290167",
    "GSE147451",
    "GSE295373",
    "GSE248752",
    "GSE261240",
    "GSE261892",
    "GSE298999",
    "GSE181757",
    "GSE93798",
    "GSE139061",
    "GSE163603",
    "GSE47184",
    "GSE30528",
    "GSE30529",
    "GSE47109",
    "GSE37460",
    "GSE35487",
    "GSE35489",
    "GSE30718",
    "GSE90861",
    "GSE95026",
    "GSE82337",
    "GSE37838",
    "GSE83486",
    "GSE19130",
    "GSE26578",
    "GSE11166",
    "GSE61739",
    "GSE217427",
]

MEANINGFUL_CHARACTERISTIC = re.compile(
    r"age|sex|gender|female|male|bmi|race|ethnicity|"
    r"diagnos|disease|condition|group|subclass|"
    r"egfr|gfr|creatinine|proteinuria|albuminuria|"
    r"fibrosis|ifta|cadi|banff|histolog|glomerulosclerosis|"
    r"kdri|remuzzi|schold|mapi|irish|"
    r"(?:^|[\s_])(?:ci|ct|cv|ah|i)(?:$|[\s_])|"
    r"progress|outcome|fail|graft loss|survival|follow.?up|"
    r"time (?:from|to|point)|months|years|duration|"
    r"treatment|immunosuppression|felzartamab|"
    r"patient|donor|subject|cohort|protocol|batch|center|"
    r"tissue|storage|source|specimen|outlier|quality",
    re.IGNORECASE,
)


def is_meaningful_characteristic(name: str) -> bool:
    return bool(MEANINGFUL_CHARACTERISTIC.search(name))


def load_title_map() -> dict[str, str]:
    path = OUT / "geo_bulk_kidney_search_raw.csv"
    if not path.exists():
        return {}
    frame = pd.read_csv(path, usecols=["accession", "title"])
    return dict(zip(frame["accession"], frame["title"]))


def download_sample_metadata(accession: str) -> tuple[Path, str]:
    METADATA_DIR.mkdir(parents=True, exist_ok=True)
    path = METADATA_DIR / f"{accession}_samples.txt"
    if path.exists() and path.stat().st_size > 0:
        return path, "cached"
    response = requests.get(
        GEO_TEXT,
        params={
            "acc": accession,
            "targ": "gsm",
            "form": "text",
            "view": "brief",
        },
        headers={"User-Agent": "PAPPA2 clinical-cohort audit"},
        timeout=120,
    )
    response.raise_for_status()
    if "^SAMPLE = " not in response.text:
        raise ValueError(f"No sample records returned for {accession}")
    path.write_text(response.text, encoding="utf-8")
    return path, f"downloaded {len(response.content)} bytes"


def parse_sample_records(text: str) -> list[dict]:
    records = []
    for chunk in text.split("^SAMPLE = ")[1:]:
        lines = chunk.splitlines()
        gsm = lines[0].strip() if lines else ""
        record = {
            "gsm": gsm,
            "title": "",
            "source": "",
            "characteristics": {},
        }
        for line in lines[1:]:
            if " = " not in line or not line.startswith("!"):
                continue
            key, value = line[1:].split(" = ", 1)
            if key == "Sample_title":
                record["title"] = value.strip()
            elif key == "Sample_source_name_ch1":
                record["source"] = value.strip()
            elif key == "Sample_characteristics_ch1":
                if ":" in value:
                    name, parsed_value = value.split(":", 1)
                    record["characteristics"][name.strip()] = (
                        parsed_value.strip()
                    )
                else:
                    record["characteristics"][
                        f"characteristic_{len(record['characteristics']) + 1}"
                    ] = value.strip()
        records.append(record)
    return records


def numeric_fraction(values: list[str]) -> float:
    usable = [
        value
        for value in values
        if value and value.lower() not in {"na", "n/a", "unknown"}
    ]
    if not usable:
        return 0.0
    numeric = sum(
        bool(
            re.fullmatch(
                r"[<>]?\s*-?\d+(?:\.\d+)?",
                value,
            )
        )
        for value in usable
    )
    return numeric / len(usable)


def summarize_dataset(
    accession: str,
    title: str,
    records: list[dict],
) -> tuple[dict, list[dict], list[dict]]:
    characteristic_values: dict[str, list[str]] = {}
    sample_rows = []
    for record in records:
        for name, value in record["characteristics"].items():
            characteristic_values.setdefault(name, []).append(value)
            sample_rows.append(
                {
                    "accession": accession,
                    "gsm": record["gsm"],
                    "sample_title": record["title"],
                    "characteristic": name,
                    "value": value,
                }
            )

    characteristic_rows = []
    for name, values in sorted(characteristic_values.items()):
        nonempty = [
            value
            for value in values
            if value and value.lower() not in {"na", "n/a", "unknown"}
        ]
        characteristic_rows.append(
            {
                "accession": accession,
                "characteristic": name,
                "nonempty": len(nonempty),
                "unique_values": len(set(nonempty)),
                "numeric_fraction": numeric_fraction(nonempty),
                "example_values": "; ".join(nonempty[:5]),
            }
        )

    names_lower = [name.lower() for name in characteristic_values]
    curated_names = [
        name
        for name in characteristic_values
        if is_meaningful_characteristic(name)
    ]
    combined = "\n".join(
        [
            title,
            *[
                f"{name}: {'; '.join(values[:20])}"
                for name, values in characteristic_values.items()
            ],
        ]
    ).lower()

    def has_key(*terms: str) -> bool:
        return any(
            term in name for name in names_lower for term in terms
        )

    def has_numeric_key(*terms: str) -> bool:
        for name, values in characteristic_values.items():
            if any(term in name.lower() for term in terms):
                if numeric_fraction(values) >= 0.7:
                    return True
        return False

    summary = {
        "accession": accession,
        "title": title,
        "n_samples": len(records),
        "n_characteristics": len(characteristic_values),
        "gfr_numeric": has_numeric_key("egfr", "gfr"),
        "gfr_keyword": "gfr" in combined,
        "creatinine_keyword": "creatinine" in combined,
        "fibrosis_keyword": any(
            term in combined
            for term in ["fibrosis", "ifta", "interstitial fibrosis"]
        ),
        "progression_keyword": any(
            term in combined
            for term in [
                "progression",
                "graft loss",
                "allograft loss",
                "renal survival",
                "kidney survival",
                "outcome",
            ]
        ),
        "proteinuria_keyword": any(
            term in combined
            for term in ["proteinuria", "albuminuria"]
        ),
        "age_key": has_key("age"),
        "sex_key": has_key("sex", "gender"),
        "disease_key": has_key(
            "diagnosis",
            "disease",
            "group",
            "condition",
        ),
        "curated_characteristic_count": len(curated_names),
        "curated_characteristic_names": "; ".join(sorted(curated_names)),
    }
    return summary, characteristic_rows, sample_rows


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    title_map = load_title_map()
    summary_rows = []
    characteristic_rows = []
    sample_rows = []
    fetch_status = {}
    for accession in SHORTLIST:
        print(f"SAMPLES {accession}", flush=True)
        try:
            path, status = download_sample_metadata(accession)
            fetch_status[accession] = status
            records = parse_sample_records(
                path.read_text(encoding="utf-8")
            )
            summary, characteristics, samples = summarize_dataset(
                accession,
                title_map.get(accession, ""),
                records,
            )
            summary_rows.append(summary)
            characteristic_rows.extend(characteristics)
            sample_rows.extend(samples)
        except Exception as error:
            fetch_status[accession] = f"ERROR: {error}"
            summary_rows.append(
                {
                    "accession": accession,
                    "title": title_map.get(accession, ""),
                    "n_samples": 0,
                    "n_characteristics": 0,
                    "error": str(error),
                }
            )
        time.sleep(0.4)

    summary_frame = pd.DataFrame(summary_rows)
    summary_frame.to_csv(
        OUT / "geo_sample_metadata_audit_summary.csv",
        index=False,
    )
    pd.DataFrame(characteristic_rows).to_csv(
        OUT / "geo_sample_characteristics_summary.csv",
        index=False,
    )
    pd.DataFrame(sample_rows).to_csv(
        OUT / "geo_sample_characteristics_long.csv",
        index=False,
    )
    (OUT / "geo_sample_metadata_fetch_status.json").write_text(
        json.dumps(fetch_status, indent=2),
        encoding="utf-8",
    )
    printable = summary_frame.copy()
    printable["title"] = printable["title"].str.slice(0, 70)
    print(printable.to_string(index=False))


if __name__ == "__main__":
    main()
