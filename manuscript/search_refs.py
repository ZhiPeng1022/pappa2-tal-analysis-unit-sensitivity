"""Search OpenAlex (and fall back to CrossRef) for citation candidates.

Usage:
    python search_refs.py "query text" [--limit 8] [--year 2020]
    python search_refs.py --doi 10.1038/s41587-021-01033-z
"""

import argparse
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request


MAILTO = "research@example.org"
USER_AGENT = f"PAPPA2-citation-check/1.0 (mailto:{MAILTO})"

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def fetch_json(url, timeout=40, attempts=3):
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    last_error = None
    for attempt in range(1, attempts + 1):
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                return json.load(response)
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            last_error = exc
            if attempt < attempts:
                time.sleep(2 * attempt)
    raise last_error


def short_authors(authorships, max_authors=3):
    names = []
    for item in authorships[:max_authors]:
        author = item.get("author") or {}
        name = author.get("display_name")
        if name:
            names.append(name)
    if not names:
        return "Unknown authors"
    suffix = " et al." if len(authorships) > max_authors else ""
    return ", ".join(names) + suffix


def format_openalex(work):
    ids = work.get("ids") or {}
    location = work.get("primary_location") or {}
    source = location.get("source") or {}
    journal = source.get("display_name") or "No journal"
    doi = (work.get("doi") or "").replace("https://doi.org/", "")
    pmid = (ids.get("pmid") or "").replace("https://pubmed.ncbi.nlm.nih.gov/", "")
    return {
        "title": work.get("display_name") or work.get("title") or "No title",
        "authors": short_authors(work.get("authorships") or []),
        "journal": journal,
        "year": work.get("publication_year"),
        "date": work.get("publication_date"),
        "doi": doi,
        "pmid": pmid,
        "citations": work.get("cited_by_count"),
        "type": work.get("type"),
    }


def search_openalex(query, limit, year=None):
    params = {
        "search": query,
        "per-page": str(limit),
        "mailto": MAILTO,
        "select": (
            "id,doi,title,display_name,publication_year,publication_date,"
            "primary_location,authorships,ids,cited_by_count,type"
        ),
    }
    if year:
        params["filter"] = f"publication_year:{year}"
    url = "https://api.openalex.org/works?" + urllib.parse.urlencode(params)
    data = fetch_json(url)
    return [format_openalex(work) for work in data.get("results", [])]


def search_crossref(query, limit, year=None):
    params = {
        "query.bibliographic": query,
        "rows": str(limit),
        "mailto": MAILTO,
    }
    if year:
        params["filter"] = f"from-pub-date:{year}-01-01,until-pub-date:{year}-12-31"
    url = "https://api.crossref.org/works?" + urllib.parse.urlencode(params)
    data = fetch_json(url)
    results = []
    for item in (data.get("message") or {}).get("items", []):
        authors = []
        for author in (item.get("author") or [])[:3]:
            given = author.get("given") or ""
            family = author.get("family") or ""
            name = (given + " " + family).strip() or author.get("name")
            if name:
                authors.append(name)
        results.append(
            {
                "title": (item.get("title") or ["No title"])[0],
                "authors": ", ".join(authors) or "Unknown authors",
                "journal": (item.get("container-title") or ["No journal"])[0],
                "year": ((item.get("issued") or {}).get("date-parts") or [[None]])[0][0],
                "date": None,
                "doi": item.get("DOI") or "",
                "pmid": "",
                "citations": item.get("is-referenced-by-count"),
                "type": item.get("type"),
            }
        )
    return results


def resolve_doi(doi):
    url = "https://api.crossref.org/works/" + urllib.parse.quote(doi)
    data = fetch_json(url)
    item = data["message"]
    authors = []
    for author in (item.get("author") or [])[:6]:
        given = author.get("given") or ""
        family = author.get("family") or ""
        name = (given + " " + family).strip() or author.get("name")
        if name:
            authors.append(name)
    return {
        "title": (item.get("title") or ["No title"])[0],
        "authors": ", ".join(authors) or "Unknown authors",
        "journal": (item.get("container-title") or ["No journal"])[0],
        "year": ((item.get("issued") or {}).get("date-parts") or [[None]])[0][0],
        "doi": item.get("DOI") or doi,
        "volume": item.get("volume"),
        "issue": item.get("issue"),
        "page": item.get("page"),
        "type": item.get("type"),
    }


def print_record(index, record):
    print(f"[{index}] {record['title']}")
    print(f"    Authors : {record['authors']}")
    print(
        "    Source  : {journal} ({year})  DOI: {doi}  PMID: {pmid}  cited: {citations}".format(
            journal=record.get("journal"),
            year=record.get("year"),
            doi=record.get("doi") or "-",
            pmid=record.get("pmid") or "-",
            citations=record.get("citations"),
        )
    )
    if record.get("volume") or record.get("page"):
        print(
            "    Biblio  : vol {volume}, issue {issue}, pages {page}".format(
                volume=record.get("volume") or "-",
                issue=record.get("issue") or "-",
                page=record.get("page") or "-",
            )
        )
    print()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("query", nargs="?", default="")
    parser.add_argument("--limit", type=int, default=6)
    parser.add_argument("--year", type=int, default=None)
    parser.add_argument("--doi", default=None)
    args = parser.parse_args()

    if args.doi:
        print_record(1, resolve_doi(args.doi))
        return

    if not args.query:
        parser.error("provide a query or --doi")

    print(f"=== OpenAlex: {args.query} ===")
    try:
        records = search_openalex(args.query, args.limit, args.year)
    except (urllib.error.URLError, TimeoutError, OSError) as exc:
        print(f"OpenAlex failed: {exc}")
        records = []

    if not records:
        print(f"=== CrossRef fallback: {args.query} ===")
        try:
            records = search_crossref(args.query, args.limit, args.year)
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            print(f"CrossRef failed: {exc}")
            records = []

    for index, record in enumerate(records, start=1):
        print_record(index, record)


if __name__ == "__main__":
    main()
