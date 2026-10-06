"""Verify reference metadata against Europe PMC by DOI."""

import json
import sys
import time
import urllib.parse
import urllib.request


DOIS = [
    "10.1016/S0140-6736(20)30045-3",
    "10.1038/s41581-018-0052-0",
    "10.1038/nrneph.2015.3",
    "10.1172/JCI45161",
    "10.1038/s41586-023-05769-3",
    "10.1681/ASN.V103628",
    "10.1172/JCI111189",
    "10.1007/s12079-015-0259-9",
    "10.15252/emmm.201506106",
    "10.1210/en.2011-0036",
    "10.1038/s41587-020-0605-1",
    "10.1172/jci.insight.123151",
    "10.1038/s41587-021-01033-z",
    "10.34067/KID.0000000602",
    "10.1126/sciadv.adv8918",
]


def fetch(doi):
    params = {
        "query": 'DOI:"%s"' % doi,
        "format": "json",
        "resultType": "core",
        "pageSize": "1",
    }
    url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + urllib.parse.urlencode(params)
    with urllib.request.urlopen(url, timeout=30) as response:
        return json.load(response)


def main():
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    for doi in DOIS:
        try:
            data = fetch(doi)
            results = data["resultList"]["result"]
            if not results:
                print("%s => NOT FOUND" % doi)
                continue
            item = results[0]
            info = item.get("journalInfo", {})
            journal = info.get("journal", {}).get("title")
            print(
                "%s | %s | %s %s;%s(%s):%s | PMID %s"
                % (
                    doi,
                    item.get("title"),
                    journal,
                    info.get("yearOfPublication"),
                    info.get("volume"),
                    info.get("issue"),
                    item.get("pageInfo"),
                    item.get("pmid"),
                )
            )
        except Exception as exc:  # noqa: BLE001
            print("%s => ERROR %s" % (doi, exc))
        time.sleep(0.3)


if __name__ == "__main__":
    main()
