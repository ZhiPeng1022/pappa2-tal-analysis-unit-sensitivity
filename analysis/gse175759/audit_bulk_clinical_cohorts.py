"""Search GEO for bulk kidney cohorts with clinical outcome variables."""

from __future__ import annotations

import json
import re
import time
from pathlib import Path

import pandas as pd
import requests


ROOT = Path(r"C:\Users\Elsa\Documents\Codex\2026-10-15\ji")
OUT = ROOT / "work" / "bulk_clinical_cohort_audit"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
MAX_RESULTS_PER_QUERY = 250

COMMON = (
    '"Homo sapiens"[Organism] AND '
    '(kidney[All Fields] OR renal[All Fields]) AND '
    'gse[Entry Type] AND '
    '("expression profiling by array"[DataSet Type] OR '
    '"expression profiling by high throughput sequencing"[DataSet Type])'
)

QUERIES = {
    "gfr": (
        '"glomerular filtration rate"[All Fields] OR '
        'GFR[All Fields] OR eGFR[All Fields]'
    ),
    "creatinine": "creatinine[All Fields]",
    "fibrosis": (
        'fibrosis[All Fields] OR "interstitial fibrosis"[All Fields]'
    ),
    "progression": (
        '"disease progression"[All Fields] OR progression[All Fields] OR '
        '"graft loss"[All Fields] OR "renal function decline"[All Fields]'
    ),
    "ckd": (
        '"chronic kidney disease"[All Fields] OR CKD[All Fields] OR '
        '"chronic renal disease"[All Fields]'
    ),
    "dkd": (
        '"diabetic kidney disease"[All Fields] OR '
        '"diabetic nephropathy"[All Fields]'
    ),
    "iga": '"IgA nephropathy"[All Fields]',
    "lupus": '"lupus nephritis"[All Fields]',
    "transplant": (
        '"kidney transplant"[All Fields] OR "renal allograft"[All Fields] OR '
        '"kidney allograft"[All Fields]'
    ),
    "biopsy": (
        'biopsy[All Fields] AND '
        '("kidney disease"[All Fields] OR renal[All Fields])'
    ),
    "proteinuria": (
        'proteinuria[All Fields] OR albuminuria[All Fields]'
    ),
}

EXCLUSION_TERMS = [
    "single cell",
    "single-cell",
    "single nucleus",
    "single-nucleus",
    "snrna",
    "scrna",
    "spatial transcript",
    "visium",
    "xenium",
    "cosmx",
    "methylation",
    "microrna",
    "mirna",
    "exosome",
    "extracellular vesicle",
    "urine",
    "urinary",
    "peripheral blood",
    "pbmc",
    "cell line",
    "organoid",
]

CLINICAL_PATTERNS = {
    "gfr": re.compile(
        r"\b(?:egfr|gfr|glomerular filtration rate)\b",
        re.IGNORECASE,
    ),
    "fibrosis": re.compile(
        r"\b(?:fibrosis|interstitial fibrosis|ifta)\b",
        re.IGNORECASE,
    ),
    "progression": re.compile(
        r"\b(?:progression|progressive|decline|graft loss|"
        r"allograft loss)\b",
        re.IGNORECASE,
    ),
    "creatinine": re.compile(r"\bcreatinine\b", re.IGNORECASE),
    "proteinuria": re.compile(
        r"\b(?:proteinuria|albuminuria)\b",
        re.IGNORECASE,
    ),
    "transplant": re.compile(
        r"\b(?:transplant|allograft)\b",
        re.IGNORECASE,
    ),
    "ckd": re.compile(
        r"\b(?:chronic kidney disease|chronic renal disease|ckd)\b",
        re.IGNORECASE,
    ),
    "dkd": re.compile(
        r"\b(?:diabetic kidney disease|diabetic nephropathy)\b",
        re.IGNORECASE,
    ),
    "iga": re.compile(r"\b(?:iga nephropathy|igan)\b", re.IGNORECASE),
    "lupus": re.compile(r"\blupus nephritis\b", re.IGNORECASE),
    "biopsy": re.compile(r"\bbiopsy\b", re.IGNORECASE),
}


def search_ids(query: str) -> tuple[list[str], int]:
    response = requests.get(
        f"{EUTILS}/esearch.fcgi",
        params={
            "db": "gds",
            "term": query,
            "retmax": MAX_RESULTS_PER_QUERY,
            "retmode": "json",
            "sort": "relevance",
        },
        timeout=60,
    )
    response.raise_for_status()
    payload = response.json()["esearchresult"]
    return payload.get("idlist", []), int(payload.get("count", 0))


def fetch_summaries(ids: list[str]) -> dict[str, dict]:
    summaries = {}
    for start in range(0, len(ids), 100):
        batch = ids[start : start + 100]
        response = requests.post(
            f"{EUTILS}/esummary.fcgi",
            data={
                "db": "gds",
                "id": ",".join(batch),
                "retmode": "json",
            },
            timeout=120,
        )
        response.raise_for_status()
        result = response.json()["result"]
        summaries.update(
            {
                uid: result[uid]
                for uid in result.get("uids", [])
            }
        )
        time.sleep(0.4)
    return summaries


def clinical_terms(text: str) -> list[str]:
    return [
        name
        for name, pattern in CLINICAL_PATTERNS.items()
        if pattern.search(text)
    ]


def exclusion_terms(text: str) -> list[str]:
    lowered = text.lower()
    return [term for term in EXCLUSION_TERMS if term in lowered]


def score_record(record: dict) -> int:
    terms = set(record["clinical_terms"])
    weights = {
        "gfr": 5,
        "fibrosis": 4,
        "progression": 4,
        "creatinine": 2,
        "proteinuria": 2,
        "transplant": 2,
        "ckd": 2,
        "dkd": 2,
        "iga": 2,
        "lupus": 2,
        "biopsy": 1,
    }
    score = sum(weights.get(term, 0) for term in terms)
    score += min(record["n_samples"], 300) / 100
    score -= 4 * len(record["exclusion_terms"])
    return round(score, 2)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    query_hits = {}
    id_queries: dict[str, list[str]] = {}
    for name, clinical_query in QUERIES.items():
        query = f"{COMMON} AND ({clinical_query})"
        ids, count = search_ids(query)
        query_hits[name] = {
            "query": query,
            "reported_count": count,
            "retrieved": len(ids),
        }
        for uid in ids:
            id_queries.setdefault(uid, []).append(name)
        print(
            f"QUERY {name}: reported={count:,} retrieved={len(ids)}",
            flush=True,
        )
        time.sleep(0.4)

    all_ids = sorted(
        id_queries,
        key=lambda uid: int(uid),
        reverse=True,
    )
    summaries = fetch_summaries(all_ids)
    records = []
    for uid in all_ids:
        item = summaries.get(uid)
        if not item or item.get("entrytype") != "GSE":
            continue
        title = item.get("title", "")
        summary = item.get("summary", "")
        combined = f"{title}\n{summary}"
        record = {
            "accession": item.get("accession", ""),
            "uid": uid,
            "title": title,
            "summary": summary,
            "taxon": item.get("taxon", ""),
            "gdstype": item.get("gdstype", ""),
            "n_samples": int(item.get("n_samples") or 0),
            "platform": item.get("gpl", ""),
            "supplementary_file": item.get("suppfile", ""),
            "pubmed_ids": ";".join(item.get("pubmedids", [])),
            "query_hits": ";".join(id_queries[uid]),
            "clinical_terms": ";".join(clinical_terms(combined)),
            "exclusion_terms": ";".join(exclusion_terms(combined)),
        }
        record["candidate_score"] = score_record(record)
        records.append(record)

    raw = pd.DataFrame(records)
    raw = raw.sort_values(
        ["candidate_score", "n_samples"],
        ascending=[False, False],
    )
    raw.to_csv(OUT / "geo_bulk_kidney_search_raw.csv", index=False)

    candidates = raw[
        raw["taxon"].eq("Homo sapiens")
        & raw["n_samples"].ge(20)
        & raw["exclusion_terms"].eq("")
        & raw["clinical_terms"].ne("")
    ].copy()
    candidates.to_csv(
        OUT / "geo_bulk_kidney_clinical_candidates.csv",
        index=False,
    )
    summary_payload = {
        "queries": query_hits,
        "unique_series": len(raw),
        "filtered_candidates": len(candidates),
        "filter_rule": (
            "human GSE; n>=20; no single-cell/spatial/urine/blood/"
            "methylation exclusion term; at least one clinical term"
        ),
    }
    (OUT / "geo_bulk_kidney_search_summary.json").write_text(
        json.dumps(summary_payload, indent=2),
        encoding="utf-8",
    )
    print(
        f"DONE unique_series={len(raw)} "
        f"filtered_candidates={len(candidates)}",
        flush=True,
    )
    print(candidates.head(30).to_string(index=False))


if __name__ == "__main__":
    main()
