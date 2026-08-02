#!/usr/bin/env python3

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def read_csv(root: Path, case: str, name: str) -> pd.DataFrame:
    return pd.read_csv(root / case / name)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    voltage_waveform = read_csv(
        args.root, "electronics_voltage", "waveform.csv"
    )
    pe_images = read_csv(args.root, "electronics_pe", "images.csv")
    nsb_waveform = read_csv(args.root, "electronics_nsb", "waveform.csv")

    voltage_pixel = int(
        voltage_waveform.groupby("pixel_id")["sample_value"].max().idxmax()
    )
    voltage = voltage_waveform[voltage_waveform.pixel_id == voltage_pixel]
    nsb_pixel = int(
        nsb_waveform.groupby("pixel_id")["sample_value"].max().idxmax()
    )
    nsb = nsb_waveform[nsb_waveform.pixel_id == nsb_pixel]
    active_pe = pe_images[
        (pe_images.primary_total_pe > 0) | (pe_images.fired_total_pe > 0)
    ]

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))

    axes[0].plot(
        voltage.time_center_ns,
        voltage.sample_value,
        marker="o",
        linewidth=1.8,
    )
    axes[0].set_xlabel("time [ns]")
    axes[0].set_ylabel("sample [mV]")
    axes[0].set_title(f"voltage mode, pixel {voltage_pixel}")
    axes[0].grid(alpha=0.25)

    x = range(len(active_pe))
    axes[1].bar(
        [value - 0.18 for value in x],
        active_pe.primary_total_pe,
        width=0.36,
        label="primary p.e.",
    )
    axes[1].bar(
        [value + 0.18 for value in x],
        active_pe.fired_total_pe,
        width=0.36,
        label="fired p.e.",
    )
    axes[1].set_xticks(list(x), active_pe.pixel_id.astype(int))
    axes[1].set_xlabel("pixel id")
    axes[1].set_ylabel("integrated p.e.")
    axes[1].set_title("p.e. mode after microcell saturation")
    axes[1].legend(frameon=False)
    axes[1].grid(axis="y", alpha=0.25)

    axes[2].step(
        nsb.time_center_ns,
        nsb.sample_value,
        where="mid",
        linewidth=1.7,
    )
    axes[2].set_xlabel("time [ns]")
    axes[2].set_ylabel("sample [mV]")
    axes[2].set_title(f"pure NSB, max pixel {nsb_pixel}")
    axes[2].grid(alpha=0.25)

    fig.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=180, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
