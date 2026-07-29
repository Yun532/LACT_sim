#!/usr/bin/env python3
"""Plot a PhotonCsv ROOT camera image through pyLAST."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np

def as_array(values) -> np.ndarray:
    if hasattr(values, "to_value"):
        return np.asarray(values.to_value(), dtype=float)
    return np.asarray(values, dtype=float)


def select_event(source, event_id: int | None):
    if event_id is None:
        return source[0]
    for index in range(len(source)):
        event = source[index]
        if int(event.event_id) == event_id:
            return event
    raise RuntimeError(f"event {event_id} is not present in the ROOT file")


def coordinates_for_view(
    pix_x_deg: np.ndarray,
    pix_y_deg: np.ndarray,
    coordinate_view: str,
) -> tuple[np.ndarray, np.ndarray, str, str, str]:
    """Return coordinates in the argument order expected by plot_camera_image.

    pyLAST's plot_camera_image places its second coordinate on the horizontal
    axis and its first coordinate on the vertical axis. LactEventSource has
    already converted canonical LACT_sim focal-plane coordinates with
    pix_x=-v and pix_y=-u.
    """
    if coordinate_view == "lact-uv":
        lact_u_deg = -pix_y_deg
        lact_v_deg = -pix_x_deg
        return (
            lact_v_deg,
            lact_u_deg,
            "LACT_sim u [deg]",
            "LACT_sim v [deg]",
            "LACT u=-pyLAST pix_y, v=-pyLAST pix_x",
        )
    return (
        pix_x_deg,
        pix_y_deg,
        "pyLAST pix_y [deg]",
        "pyLAST pix_x [deg]",
        "pyLAST display: horizontal=pix_y, vertical=pix_x",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root_file", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--event-id", type=int)
    parser.add_argument("--telescope-id", type=int, default=19)
    parser.add_argument(
        "--coordinate-view",
        choices=("pylast", "lact-uv"),
        default="pylast",
        help=(
            "pylast uses the current source-offset camera display; lact-uv "
            "plots the same pyLAST event in LACT_sim output-plane u/v"
        ),
    )
    args = parser.parse_args()

    from pylast.io import LactEventSource
    from pylast.visualize.visualize import plot_camera_image

    source = LactEventSource(str(args.root_file), max_events=-1)
    event = select_event(source, args.event_id)
    telescope_id = args.telescope_id
    if telescope_id not in event.simulation.tels:
        raise RuntimeError(
            f"telescope {telescope_id} has no simulated image; "
            f"available={sorted(event.simulation.tels)}"
        )

    telescope = source.subarray.tels[telescope_id]
    geometry = telescope.camera.geometry
    image_pe = np.asarray(
        event.simulation.tels[telescope_id].true_image, dtype=float
    )
    focal_length_m = float(telescope.optics.effective_focal_length)
    if not np.isfinite(focal_length_m) or focal_length_m <= 0.0:
        focal_length_m = float(telescope.optics.equivalent_focal_length)
    pix_x_deg = np.degrees(
        np.arctan2(as_array(geometry.pix_x), focal_length_m)
    )
    pix_y_deg = np.degrees(
        np.arctan2(as_array(geometry.pix_y), focal_length_m)
    )
    pixel_size_m = float(np.sqrt(np.median(as_array(geometry.pix_area))))
    pixel_size_deg = float(np.degrees(np.arctan2(pixel_size_m, focal_length_m)))
    mask = image_pe > 0.0
    plot_x, plot_y, xlabel, ylabel, coordinate_note = coordinates_for_view(
        pix_x_deg,
        pix_y_deg,
        args.coordinate_view,
    )

    axis = plot_camera_image(
        plot_x,
        plot_y,
        pixel_size_deg,
        image_pe,
        mask=mask,
        vmin=0.0,
        vmax=float(image_pe.max()) if np.any(mask) else 1.0,
        title=(
            f"LACT ROOT full camera: event {int(event.event_id)}, "
            f"telescope {telescope_id}"
        ),
        pixel_shape="square",
    )
    axis.set_xlabel(xlabel)
    axis.set_ylabel(ylabel)
    axis.text(
        0.02,
        0.02,
        f"Total = {image_pe.sum():.0f} p.e.\n{coordinate_note}",
        transform=axis.transAxes,
        ha="left",
        va="bottom",
        bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none"},
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    axis.figure.savefig(args.output, dpi=220, bbox_inches="tight")
    plt.close(axis.figure)
    print(args.output)


if __name__ == "__main__":
    main()
