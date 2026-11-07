"""Download GSE175759 per-sample expression files from GEO."""

from __future__ import annotations

import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd


ROOT = Path(r"C:\Users\Elsa\Documents\Codex\2026-10-15\ji")
AUDIT = ROOT / "work" / "bulk_clinical_cohort_audit"
FILELIST = AUDIT / "GSE175759_filelist.txt"
OUT = AUDIT / "GSE175759_samples"
WORKERS = 6


def load_files() -> pd.DataFrame:
    rows = []
    for line in FILELIST.read_text(encoding="utf-8").splitlines():
        if not line.startswith("File\t"):
            continue
        fields = line.split("\t")
        file_name = fields[1]
        gsm = file_name.split("_", 1)[0]
        expected_size = int(fields[3])
        rows.append(
            {
                "gsm": gsm,
                "file_name": file_name,
                "expected_size": expected_size,
            }
        )
    return pd.DataFrame(rows)


def download_one(row: dict, index: int) -> dict:
    gsm = row["gsm"]
    file_name = row["file_name"]
    expected_size = int(row["expected_size"])
    output = OUT / file_name
    if output.exists() and output.stat().st_size == expected_size:
        return {
            "gsm": gsm,
            "file_name": file_name,
            "size": expected_size,
            "status": "cached",
        }
    ip = "130.14.250.7" if index % 2 == 0 else "130.14.250.31"
    url = (
        "https://ftp.ncbi.nlm.nih.gov/geo/samples/"
        f"{gsm[:-3]}nnn/{gsm}/suppl/{file_name}"
    )
    command = [
        "curl.exe",
        "-4",
        "--ssl-no-revoke",
        "-fsSL",
        "--http1.1",
        "--connect-timeout",
        "30",
        "--speed-limit",
        "512",
        "--speed-time",
        "120",
        "--retry",
        "3",
        "--retry-delay",
        "3",
        "--resolve",
        f"ftp.ncbi.nlm.nih.gov:443:{ip}",
        "--max-time",
        "300",
        url,
        "-o",
        str(output),
    ]
    subprocess.run(
        command,
        check=True,
        capture_output=True,
        timeout=320,
    )
    size = output.stat().st_size
    if size != expected_size:
        raise IOError(
            f"{file_name}: expected {expected_size}, received {size}"
        )
    return {
        "gsm": gsm,
        "file_name": file_name,
        "size": size,
        "status": f"downloaded via {ip}",
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    files = load_files()
    if len(files) != 90:
        raise ValueError(f"Expected 90 sample files, found {len(files)}")
    rows = []
    print(f"START files={len(files)} workers={WORKERS}", flush=True)
    with ThreadPoolExecutor(max_workers=WORKERS) as executor:
        futures = {
            executor.submit(download_one, row, index): row["gsm"]
            for index, row in files.reset_index(drop=True).iterrows()
        }
        completed = 0
        for future in as_completed(futures):
            completed += 1
            gsm = futures[future]
            try:
                rows.append(future.result())
            except Exception as error:
                rows.append(
                    {
                        "gsm": gsm,
                        "file_name": "",
                        "size": 0,
                        "status": f"ERROR: {error}",
                    }
                )
            if completed % 5 == 0 or completed == len(futures):
                print(
                    f"PROGRESS {completed}/{len(futures)}",
                    flush=True,
                )
    report = pd.DataFrame(rows).sort_values("gsm")
    report.to_csv(
        AUDIT / "GSE175759_sample_download_status.csv",
        index=False,
    )
    failures = report[
        report["status"].str.startswith("ERROR")
    ]
    print(
        f"DONE downloaded_or_cached={len(report) - len(failures)} "
        f"failures={len(failures)}",
        flush=True,
    )
    if not failures.empty:
        print(failures.to_string(index=False))


if __name__ == "__main__":
    main()
