#!/usr/bin/env python3
"""Plot the minimal-CSV outputs.

The whiteboard and ideal photon-count figures are optical diagnostics.  The
full expected-p.e. image is rendered by pylast's camera plotting function, so
its axes, camera orientation, inactive-pixel treatment, and angular units are
the same as the normal CORSIKA/ROOT analysis path.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

from pylast.visualize.visualize import plot_camera_image


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
        xlabel="Camera x (cm)",
        ylabel="Camera y (cm)",
        title=title,
        aspect="equal",
    )
    figure.colorbar(collection, ax=axis, label=label)
    figure.tight_layout()
    figure.savefig(output, dpi=220)
    plt.close(figure)


def plot_pylast_camera(
    camera_path: Path,
    pixels_path: Path,
    output: Path,
    focal_length_m: float,
) -> None:
    """Render the full camera response using pylast's native display logic."""
    camera = read_rows(camera_path)
    by_id = {
        int(row["pixel_id"]): float(row["pe"]) for row in read_rows(pixels_path)
    }
    pix_x_m = np.asarray([float(row["x_m"]) for row in camera], dtype=float)
    pix_y_m = np.asarray([float(row["y_m"]) for row in camera], dtype=float)
    image_pe = np.asarray(
        [by_id.get(int(row["id"]), 0.0) for row in camera], dtype=float
    )
    pixel_size_m = float(camera[0]["size_m"])

    # Match LactEventSource::load_telescopes: LACT_sim stores focal-plane hit
    # coordinates, while pyLAST displays source-offset camera coordinates.
    pix_x_deg = np.degrees(np.arctan2(-pix_x_m, focal_length_m))
    pix_y_deg = np.degrees(np.arctan2(-pix_y_m, focal_length_m))
    pixel_size_deg = float(
        np.degrees(np.arctan2(pixel_size_m, focal_length_m))
    )
    mask = image_pe > 0.0

    axis = plot_camera_image(
        pix_x_deg,
        pix_y_deg,
        pixel_size_deg,
        image_pe,
        mask=mask,
        vmin=0.0,
        vmax=float(np.max(image_pe)) if np.any(mask) else 1.0,
        title="Full camera response (pylast convention, no waveform)",
        pixel_shape="square",
    )
    axis.text(
        0.02,
        0.02,
        f"Total expected p.e. = {image_pe.sum():.2f}",
        transform=axis.transAxes,
        ha="left",
        va="bottom",
        fontsize=9,
        bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none"},
    )
    figure = axis.figure
    figure.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--hits", type=Path, required=True)
    parser.add_argument("--photon-pixels", type=Path, required=True)
    parser.add_argument("--pe-pixels", type=Path, required=True)
    parser.add_argument("--camera", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--focal-length-m",
        type=float,
        default=8.0,
        help="effective focal length used by pylast's angular camera view",
    )
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
    plot_pylast_camera(
        args.camera,
        args.pe_pixels,
        args.output_dir / "minimal_photons_camera_pe.png",
        args.focal_length_m,
    )


if __name__ == "__main__":
    main()
