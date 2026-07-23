#!/usr/bin/env python3
"""Compare a minimal PhotonCsv camera image with direct CORSIKA in pyLAST.

The comparison aligns pixels by their physical camera coordinates, reports raw
and shape-only metrics, and renders both normalized images with pyLAST's native
camera plotting convention.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np

from pylast.io import LactEventSource
from pylast.visualize.visualize import plot_camera_image


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def quantity_values(values, unit: str) -> np.ndarray:
    if hasattr(values, "to_value"):
        return np.asarray(values.to_value(unit), dtype=float)
    return np.asarray(values, dtype=float)


def read_root_observation(
    path: Path,
    event_id: int,
    telescope_id: int,
    camera_rows: list[dict[str, str]],
    value_branch: str = "image_cherenkov_pe",
) -> np.ndarray:
    """Read the ROOT pixel-id/value pairs without assuming an array order."""
    import ROOT

    root_file = ROOT.TFile.Open(str(path))
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"cannot open ROOT file: {path}")
    tree = root_file.Get("observations")
    if tree is None:
        root_file.Close()
        raise RuntimeError("ROOT file has no observations tree")
    values_by_id: dict[int, float] | None = None
    for entry in range(tree.GetEntries()):
        tree.GetEntry(entry)
        if int(tree.event_id) != event_id or int(tree.telescope_id) != telescope_id:
            continue
        values_by_id = {
            int(pixel_id): float(value)
            for pixel_id, value in zip(
                tree.pixel_id, getattr(tree, value_branch)
            )
        }
        break
    root_file.Close()
    if values_by_id is None:
        raise RuntimeError(
            f"ROOT observation event={event_id}, telescope={telescope_id} not found"
        )
    return np.asarray(
        [values_by_id.get(int(row["id"]), 0.0) for row in camera_rows],
        dtype=float,
    )


def weighted_shape(
    image: np.ndarray, display_x: np.ndarray, display_y: np.ndarray
) -> dict[str, float]:
    total = float(image.sum())
    if total <= 0.0:
        raise RuntimeError("cannot measure an empty image")
    x_mean = float(np.dot(image, display_x) / total)
    y_mean = float(np.dot(image, display_y) / total)
    dx = display_x - x_mean
    dy = display_y - y_mean
    covariance = np.asarray(
        [
            [np.dot(image, dx * dx), np.dot(image, dx * dy)],
            [np.dot(image, dx * dy), np.dot(image, dy * dy)],
        ],
        dtype=float,
    ) / total
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    major = eigenvectors[:, int(np.argmax(eigenvalues))]
    angle_deg = float(np.degrees(np.arctan2(major[1], major[0])))
    return {
        "centroid_display_x_deg": x_mean,
        "centroid_display_y_deg": y_mean,
        "length_deg": float(np.sqrt(max(eigenvalues))),
        "width_deg": float(np.sqrt(min(eigenvalues))),
        "major_axis_angle_deg_mod180": angle_deg % 180.0,
    }


def angular_difference_180(first: float, second: float) -> float:
    difference = abs(first - second) % 180.0
    return min(difference, 180.0 - difference)


def render_pylast(
    pix_x_deg: np.ndarray,
    pix_y_deg: np.ndarray,
    image: np.ndarray,
    pixel_size_deg: float,
    title: str,
    output: Path,
    vmax_fraction_percent: float,
) -> None:
    fraction_percent = image / image.sum() * 100.0
    mask = fraction_percent > 0.0
    axis = plot_camera_image(
        pix_x_deg,
        pix_y_deg,
        pixel_size_deg,
        fraction_percent,
        mask=mask,
        vmin=0.0,
        vmax=vmax_fraction_percent,
        title=title,
        pixel_shape="square",
    )
    colorbar_axis = axis.figure.axes[-1]
    colorbar_axis.set_ylabel("Fraction of image [%]")
    axis.text(
        0.02,
        0.02,
        f"Raw total = {image.sum():.2f} p.e.",
        transform=axis.transAxes,
        ha="left",
        va="bottom",
        fontsize=9,
        bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none"},
    )
    axis.figure.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(axis.figure)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--corsika-root", required=True, type=Path)
    parser.add_argument("--minimal-pixels", required=True, type=Path)
    parser.add_argument("--camera", required=True, type=Path)
    parser.add_argument("--telescope-id", type=int, default=19)
    parser.add_argument("--focal-length-m", type=float, default=8.0)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    source = LactEventSource(str(args.corsika_root), max_events=1)
    event = source[0]
    if args.telescope_id not in event.simulation.tels:
        raise RuntimeError(
            f"telescope {args.telescope_id} is absent; "
            f"available={sorted(event.simulation.tels)}"
        )
    pylast_actual = np.asarray(
        event.simulation.tels[args.telescope_id].true_image, dtype=float
    )
    geometry = source.subarray.tels[args.telescope_id].camera.geometry
    actual_x_m = quantity_values(geometry.pix_x, "m")
    actual_y_m = quantity_values(geometry.pix_y, "m")

    camera_rows = read_rows(args.camera)
    root_actual_config_order = read_root_observation(
        args.corsika_root,
        int(event.event_id),
        args.telescope_id,
        camera_rows,
    )
    minimal_by_id = {
        int(row["pixel_id"]): float(row["pe"])
        for row in read_rows(args.minimal_pixels)
    }
    config_x_m = np.asarray([float(row["x_m"]) for row in camera_rows])
    config_y_m = np.asarray([float(row["y_m"]) for row in camera_rows])
    config_pe = np.asarray(
        [minimal_by_id.get(int(row["id"]), 0.0) for row in camera_rows]
    )

    if actual_x_m.size != config_x_m.size:
        raise RuntimeError("pyLAST and configuration camera sizes differ")
    # LactEventSource deliberately converts LACT_sim focal-plane coordinates
    # to pyLAST source-offset coordinates by negating both transverse axes.
    maximum_coordinate_error_m = float(
        np.max(
            np.hypot(
                actual_x_m + config_x_m,
                actual_y_m + config_y_m,
            )
        )
    )
    if maximum_coordinate_error_m > 1.0e-6:
        raise RuntimeError(
            "pyLAST camera order/coordinate convention does not match the "
            f"configuration: maximum error={maximum_coordinate_error_m} m"
        )
    actual = root_actual_config_order
    minimal = config_pe

    pix_x_deg = np.degrees(np.arctan2(actual_x_m, args.focal_length_m))
    pix_y_deg = np.degrees(np.arctan2(actual_y_m, args.focal_length_m))
    display_x_deg = pix_y_deg
    display_y_deg = pix_x_deg
    pixel_size_m = float(camera_rows[0]["size_m"])
    pixel_size_deg = float(
        np.degrees(np.arctan2(pixel_size_m, args.focal_length_m))
    )

    actual_shape = weighted_shape(actual, display_x_deg, display_y_deg)
    minimal_shape = weighted_shape(minimal, display_x_deg, display_y_deg)
    actual_fraction = actual / actual.sum()
    minimal_fraction = minimal / minimal.sum()
    cosine = float(
        np.dot(actual_fraction, minimal_fraction)
        / np.sqrt(
            np.dot(actual_fraction, actual_fraction)
            * np.dot(minimal_fraction, minimal_fraction)
        )
    )
    correlation = float(np.corrcoef(actual_fraction, minimal_fraction)[0, 1])
    best_scale = float(np.dot(actual, minimal) / np.dot(minimal, minimal))
    scaled_residual = actual - best_scale * minimal

    summary = {
        "event_id": int(event.event_id),
        "telescope_id": args.telescope_id,
        "comparison_scope": (
            "shape-only metrics compare normalized images; raw intensity is "
            "not expected to match because the minimal CSV uses one photon per "
            "CORSIKA bunch row and one configured wavelength"
        ),
        "corsika_root_signal_branch": "observations.image_cherenkov_pe",
        "pixel_alignment": {
            "pixels": int(actual.size),
            "maximum_coordinate_error_m": maximum_coordinate_error_m,
            "pylast_true_image_equals_rounded_root_cherenkov_mapping": bool(
                np.array_equal(pylast_actual, np.rint(actual))
            ),
            "pylast_true_image_vs_root_cherenkov_mapping_correlation": float(
                np.corrcoef(pylast_actual, actual)[0, 1]
            ),
            "camera_coordinate_rule": (
                "pyLAST pix_x=-LACT_sim x_m, pix_y=-LACT_sim y_m; pixel-id "
                "array order is unchanged"
            ),
        },
        "corsika": {
            "total_pe": float(actual.sum()),
            "nonzero_pixels": int(np.count_nonzero(actual)),
            **actual_shape,
        },
        "minimal_csv": {
            "total_expected_pe": float(minimal.sum()),
            "nonzero_pixels": int(np.count_nonzero(minimal)),
            **minimal_shape,
        },
        "shape_comparison": {
            "cosine_similarity": cosine,
            "pearson_pixel_correlation": correlation,
            "centroid_separation_deg": float(
                np.hypot(
                    actual_shape["centroid_display_x_deg"]
                    - minimal_shape["centroid_display_x_deg"],
                    actual_shape["centroid_display_y_deg"]
                    - minimal_shape["centroid_display_y_deg"],
                )
            ),
            "major_axis_angle_difference_deg": angular_difference_180(
                actual_shape["major_axis_angle_deg_mod180"],
                minimal_shape["major_axis_angle_deg_mod180"],
            ),
            "best_raw_intensity_scale_corsika_per_minimal": best_scale,
            "relative_l2_residual_after_best_scale": float(
                np.linalg.norm(scaled_residual) / np.linalg.norm(actual)
            ),
        },
    }

    shared_vmax_fraction_percent = float(
        max((actual / actual.sum()).max(), (minimal / minimal.sum()).max())
        * 100.0
    )
    render_pylast(
        pix_x_deg,
        pix_y_deg,
        actual,
        pixel_size_deg,
        f"Direct CORSIKA — event {event.event_id}, telescope {args.telescope_id}",
        args.output_dir / "event1909_tel19_direct_corsika_pylast.png",
        shared_vmax_fraction_percent,
    )
    render_pylast(
        pix_x_deg,
        pix_y_deg,
        minimal,
        pixel_size_deg,
        f"Minimal six-column CSV — event {event.event_id}, telescope {args.telescope_id}",
        args.output_dir / "event1909_tel19_minimal_csv_pylast.png",
        shared_vmax_fraction_percent,
    )
    (args.output_dir / "event1909_tel19_comparison.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
