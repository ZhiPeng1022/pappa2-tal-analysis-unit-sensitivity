"""Basic programmatic QA for generated figure PNGs."""

import os
from pathlib import Path

import matplotlib.image as mpimg
import numpy as np


PROJECT_ROOT = Path(
    os.environ.get(
        "PAPPA2_PROJECT_ROOT",
        Path(__file__).resolve().parents[2],
    )
)
OUT = Path(
    os.environ.get(
        "PAPPA2_OUTPUT_ROOT",
        PROJECT_ROOT / "outputs",
    )
)
FILES = [
    "Fig1_study_design.png",
    "Fig2_TAL_neighborhood_remodeling.png",
    "Fig3_PAPPA2_pseudobulk.png",
    "Fig4_two_level_association.png",
    "Fig5_KPMP_eGFR_specificity.png",
    "Fig6_KPMP_spatial_crossmodality.png",
    "Fig7_external_GFR_associations.png",
    "Fig8_bounded_model.png",
]


def main():
    for name in FILES:
        path = OUT / name
        image = mpimg.imread(path)
        rgb = image[:, :, :3]
        non_white = (rgb < 0.98).any(axis=2).mean()
        print(f"{name}")
        print(f"  size          : {rgb.shape[1]} x {rgb.shape[0]} px")
        print(f"  mean / std    : {rgb.mean():.3f} / {rgb.std():.3f}")
        print(f"  non-white     : {non_white * 100:.1f}%")
        print(f"  blank         : {'YES' if non_white < 0.005 else 'no'}")


if __name__ == "__main__":
    main()
