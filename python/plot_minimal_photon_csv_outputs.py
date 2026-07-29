#!/usr/bin/env python3
"""Plot the two pure-optics outputs from a minimal PhotonCsv."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def plot_whiteboard(hits_path: Path, output: Path) -> None:
    rows = read_rows(hits_path)
    u = np.asarray([float(row["u_m"]) for row in rows], dtype=float)
    v = np.asarray([float(row["v_m"]) for row in rows], dtype=float)
    figure, axis = plt.subplots(figsize=(7.2, 6.2))
    image = axis.hist2d(
        u * 100.0,
        v * 100.0,
        bins=160,
        range=((-65.0, 65.0), (-65.0, 65.0)),
        cmap="viridis",
    )
    figure.colorbar(image[3], ax=axis, label="Photons / bin")
    axis.set(
        xlabel="Focal-plane u (cm)",
        ylabel="Focal-plane v (cm)",
        title="Pure optical whiteboard",
        aspect="equal",
    )
    figure.tight_layout()
    figure.savefig(output, dpi=220)
    plt.close(figure)


def plot_camera(
    camera_path: Path,
    pixels_path: Path,
    field: str,
    label: str,
    title: str,
    output: Path,
) -> None:
    camera = read_rows(camera_path)
    values = {
        int(row["pixel_id"]): float(row[field]) for row in read_rows(pixels_path)
    }
    patches = []
    colors = []
    for row in camera:
        pixel_id = int(row["id"])
        x = float(row["x_m"]) * 100.0
        y = float(row["y_m"]) * 100.0
        size = float(row["size_m"]) * 100.0
        patches.append(Rectangle((x - size / 2.0, y - size / 2.0), size, size))
        colors.append(values.get(pixel_id, 0.0))

    figure, axis = plt.subplots(figsize=(7.2, 6.2))
    collection = PatchCollection(
        patches,
        array=np.asarray(colors, dtype=float),
        cmap="viridis",
        edgecolor="none",
    )
    axis.add_collection(collection)
    axis.autoscale_view()
    axis.set(
        xlabel="LACT focal-plane u (cm)",
        ylabel="LACT focal-plane v (cm)",
        title=title,
        aspect="equal",
    )
    figure.colorbar(collection, ax=axis, label=label)
    figure.tight_layout()
    figure.savefig(output, dpi=220)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mode",
        choices=("optics",),
        default="optics",
        help="write the whiteboard and raw photon-count camera diagnostics",
    )
    parser.add_argument("--hits", type=Path, required=True)
    parser.add_argument("--photon-pixels", type=Path, required=True)
    parser.add_argument("--camera", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    plot_whiteboard(
        args.hits, args.output_dir / "minimal_photons_whiteboard.png"
    )
    plot_camera(
        args.camera,
        args.photon_pixels,
        "photon_count",
        "Photons",
        "Ideal camera photon count",
        args.output_dir / "minimal_photons_camera_counts.png",
    )


if __name__ == "__main__":
    main()
