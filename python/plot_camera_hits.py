#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle


def read_camera_pixels(path):
    pixels = []
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            pixels.append({
                "id": int(row["id"]),
                "x_m": float(row["x_m"]),
                "y_m": float(row["y_m"]),
                "size_m": float(row["size_m"]),
            })
    return pixels


def main():
    parser = argparse.ArgumentParser(description="Plot camera photon counts per pixel.")
    parser.add_argument("--hits", required=True, help="surface hit CSV from run_optical_sim")
    parser.add_argument("--camera-csv", required=True, help="camera pixel CSV")
    parser.add_argument("--output", default="camera_hits.png")
    parser.add_argument("--summary-csv", default=None, help="optional per-pixel summary CSV")
    parser.add_argument("--dpi", type=int, default=300)
    args = parser.parse_args()

    hits = pd.read_csv(args.hits)
    if "pixel_id" not in hits.columns:
        raise SystemExit("hit CSV has no pixel_id column; rerun with camera.enabled=true")
    valid = hits[(hits["hit_camera"] == 1) & (hits["pixel_id"] >= 0)]
    weights = valid["weight"].to_numpy(float) * valid["relative_efficiency"].to_numpy(float)
    counts = {}
    for pid, w in zip(valid["pixel_id"].to_numpy(int), weights):
        counts[pid] = counts.get(pid, 0.0) + w

    pixels = read_camera_pixels(Path(args.camera_csv))
    pixel_by_id = {p["id"]: p for p in pixels}
    patches = []
    values = []
    for p in pixels:
        size = p["size_m"]
        patches.append(Rectangle((p["x_m"] - 0.5 * size, p["y_m"] - 0.5 * size), size, size))
        values.append(counts.get(p["id"], 0.0))

    total_weight = float(np.sum(weights)) if len(weights) else 0.0
    filled_pixels = sum(1 for v in values if v > 0.0)
    brightest_id = max(counts, key=counts.get) if counts else -1
    brightest_count = counts.get(brightest_id, 0.0)
    brightest_pixel = pixel_by_id.get(brightest_id)

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.labelsize": 10,
        "axes.titlesize": 11,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "figure.dpi": args.dpi,
        "savefig.dpi": args.dpi,
    })
    fig, ax = plt.subplots(figsize=(6.8, 6.2))
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_under((1.0, 1.0, 1.0, 0.0))
    positive = [v for v in values if v > 0.0]

    collection = PatchCollection(
        patches,
        cmap=cmap,
        edgecolor=(0.15, 0.15, 0.15, 0.25),
        linewidth=0.25,
    )
    collection.set_array(np.asarray(values, dtype=float))
    collection.set_clim(vmin=0.5 if positive else 0.5, vmax=max(positive) if positive else 1.0)
    ax.add_collection(collection)
    ax.set_aspect("equal", adjustable="box")
    ax.autoscale_view()
    ax.set_xlabel("camera x [m]")
    ax.set_ylabel("camera y [m]")
    ax.set_title("Camera photon image")
    ax.grid(True, alpha=0.18, linewidth=0.5)
    cbar = fig.colorbar(collection, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("weighted photons")
    ax.text(
        0.02,
        0.02,
        f"pixels = {len(pixels)}\nfilled = {len(counts)}\nphotons = {weights.sum():.0f}",
        transform=ax.transAxes,
        fontsize=8,
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="0.82", alpha=0.9),
    )
    fig.tight_layout()
    fig.savefig(args.output, bbox_inches="tight")
    if args.summary_csv:
        with open(args.summary_csv, "w", newline="") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=["pixel_id", "x_m", "y_m", "weighted_count"],
            )
            writer.writeheader()
            for p, value in zip(pixels, values):
                writer.writerow({
                    "pixel_id": p["id"],
                    "x_m": p["x_m"],
                    "y_m": p["y_m"],
                    "weighted_count": value,
                })

    print(f"Camera pixels = {len(pixels)}")
    print(f"Filled pixels = {filled_pixels}")
    print(f"Weighted photons in camera = {total_weight:.6f}")
    if brightest_pixel is not None:
        print(
            "Brightest pixel = "
            f"{brightest_id} at "
            f"({brightest_pixel['x_m']:.6f}, {brightest_pixel['y_m']:.6f}) m, "
            f"count = {brightest_count:.6f}"
        )
    print(f"Saved camera image = {args.output}")
    if args.summary_csv:
        print(f"Saved camera summary = {args.summary_csv}")
    if "agg" not in plt.get_backend().lower():
        plt.show()


if __name__ == "__main__":
    main()
