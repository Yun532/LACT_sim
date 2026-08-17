#!/usr/bin/env python3
"""Plot light-collector angular response from scan_light_collector_angular_response."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("csv", help="CSV from scan_light_collector_angular_response")
    parser.add_argument(
        "--output",
        default=None,
        help="output image path; defaults to <csv stem>.png",
    )
    parser.add_argument("--dpi", type=int, default=350)
    args = parser.parse_args()

    df = pd.read_csv(args.csv)
    angle = df["angle_deg"].to_numpy(float)
    input_photons = df["input_photons"].to_numpy(float)
    hit_photons = df["hit_sipm_photons"].to_numpy(float)
    weighted_photons = df["weighted_photons"].to_numpy(float)
    geom = df["geometric_acceptance"].to_numpy(float)
    weighted = df["weighted_acceptance"].to_numpy(float)
    mean_weight = df["mean_collector_weight"].to_numpy(float)
    mean_refl = df["mean_reflections"].to_numpy(float)

    out = Path(args.output) if args.output else Path(args.csv).with_suffix(".png")
    out.parent.mkdir(parents=True, exist_ok=True)

    plt.rcParams.update({
        "font.size": 11,
        "axes.linewidth": 1.1,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
    })

    fig, axes = plt.subplots(2, 1, figsize=(7.2, 7.4), sharex=True)

    ax = axes[0]
    ax.plot(angle, hit_photons, color="#1f77b4", lw=2.0, label="hit SiPM photons")
    ax.plot(angle, weighted_photons, color="#d62728", lw=2.0,
            label="weighted photons")
    ax.fill_between(angle, 0, input_photons, color="0.92", label="input photons")
    ax.set_ylabel("photons per angle")
    ax.legend(frameon=True, framealpha=1.0, edgecolor="0.2")
    ax.grid(alpha=0.25, lw=0.6)

    ax = axes[1]
    ax.plot(angle, geom, color="#1f77b4", lw=2.0, label="geometric acceptance")
    ax.plot(angle, weighted, color="#d62728", lw=2.0, label="weighted acceptance")
    ax.plot(angle, mean_weight, color="#2ca02c", lw=1.8, ls="--",
            label="mean weight of accepted photons")
    ax.set_xlabel("incidence angle relative to collector axis (deg)")
    ax.set_ylabel("fraction / weight")
    ax.set_ylim(bottom=0.0)
    ax.grid(alpha=0.25, lw=0.6)

    ax2 = ax.twinx()
    ax2.plot(angle, mean_refl, color="#9467bd", lw=1.8, ls=":",
             label="mean reflections")
    ax2.set_ylabel("mean reflections")

    lines, labels = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines + lines2, labels + labels2, frameon=True,
              framealpha=1.0, edgecolor="0.2", loc="best")

    fig.suptitle("Square Bezier Light Collector Angular Response", y=0.985)
    fig.tight_layout()
    fig.savefig(out, dpi=args.dpi)
    print(f"saved {out}")


if __name__ == "__main__":
    main()
