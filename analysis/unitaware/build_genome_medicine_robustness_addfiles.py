"""Build Additional files 28 and 29 for the Genome Medicine draft."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(r"C:\Users\Elsa\Documents\Codex\2026-10-15\ji")
SOURCE = ROOT / "outputs" / "unitaware_robustness"
OUT = ROOT / "outputs" / "genome_medicine_additional_v4"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    thresholds = pd.read_csv(SOURCE / "threshold_sensitivity.csv")
    thresholds.to_csv(OUT / "Additional file 28.csv", index=False)

    summary = pd.read_csv(SOURCE / "pairing_permutation_summary.csv")
    null_runs = pd.read_csv(SOURCE / "pairing_permutation_null.csv")
    with pd.ExcelWriter(
        OUT / "Additional file 29.xlsx",
        engine="openpyxl",
    ) as writer:
        summary.to_excel(writer, sheet_name="summary", index=False)
        null_runs.to_excel(writer, sheet_name="null_runs", index=False)
    print(OUT / "Additional file 28.csv")
    print(OUT / "Additional file 29.xlsx")


if __name__ == "__main__":
    main()
