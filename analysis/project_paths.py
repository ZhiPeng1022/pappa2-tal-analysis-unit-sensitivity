import os
from pathlib import Path


PROJECT_ROOT = Path(
    os.environ.get(
        "PAPPA2_PROJECT_ROOT",
        Path(__file__).resolve().parents[2],
    )
)
DATA_ROOT = Path(
    os.environ.get(
        "PAPPA2_DATA_ROOT",
        PROJECT_ROOT / "data" / "RAMP3",
    )
)
OUTPUT_ROOT = Path(
    os.environ.get(
        "PAPPA2_OUTPUT_ROOT",
        PROJECT_ROOT / "outputs",
    )
)
WORK_ROOT = Path(
    os.environ.get(
        "PAPPA2_WORK_ROOT",
        PROJECT_ROOT / "work",
    )
)
