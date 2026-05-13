#!/usr/bin/env python3
"""Convert a simple ROOT `photons` tree into LACT_sim PhotonCsv.

The expected ROOT branches are:

    photon_x, photon_y, photon_z, photon_u, photon_v, photon_weight, photon_lambda

For the standalone laser/ROOT test, the input coordinates are treated as a
CORSIKA-like NWU frame in centimeters.  Rays are projected onto the selected
telescope local entrance plane (local z=0) and can optionally be recentered on
the photon-bundle centroid before being passed to the normal PhotonCsv raytrace.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np
import uproot


OUTPUT_COLUMNS = [
    "x_m",
    "y_m",
    "z_m",
    "dir_x",
    "dir_y",
    "dir_z",
    "time_ns",
    "wavelength_nm",
    "weight",
    "multiplicity",
    "event_id",
    "telescope_id",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert ROOT photons tree to LACT_sim local-frame PhotonCsv."
    )
    parser.add_argument("input_root")
    parser.add_argument("output_csv")
    parser.add_argument("--tree", default="photons")
    parser.add_argument("--entry-start", type=int, default=0)
    parser.add_argument("--entry-stop", type=int, default=None)
    parser.add_argument("--x-column", default="photon_x")
    parser.add_argument("--y-column", default="photon_y")
    parser.add_argument("--z-column", default="photon_z")
    parser.add_argument("--u-column", default="photon_u")
    parser.add_argument("--v-column", default="photon_v")
    parser.add_argument("--weight-column", default="photon_weight")
    parser.add_argument("--wavelength-column", default="photon_lambda")
    parser.add_argument("--position-scale", type=float, default=0.01, help="cm -> m by default")
    parser.add_argument("--telescope-x-cm", type=float, default=-14464.5)
    parser.add_argument("--telescope-y-cm", type=float, default=-2873.45)
    parser.add_argument("--telescope-z-cm", type=float, default=119.7)
    parser.add_argument("--telescope-az-deg", type=float, default=100.156)
    parser.add_argument("--telescope-zenith-deg", type=float, default=25.802)
    parser.add_argument(
        "--azimuth-convention",
        choices=["corsika_nwu"],
        default="corsika_nwu",
        help="CORSIKA x=north, y=west, z=up; az=90 points east (-y).",
    )
    parser.add_argument(
        "--coordinate-mode",
        choices=["pointing_transform", "subtract_telescope_only"],
        default="pointing_transform",
        help=(
            "pointing_transform rotates global ROOT coordinates into the telescope frame; "
            "subtract_telescope_only subtracts telX/telY and assumes photon_u/v are already "
            "in the telescope-interface frame."
        ),
    )
    parser.add_argument(
        "--center-mode",
        choices=["none", "mean", "median"],
        default="median",
        help="recenter projected entrance-plane coordinates for standalone optical tests",
    )
    parser.add_argument(
        "--max-radius-m",
        type=float,
        default=8.0,
        help="keep only projected photons within this radius after centering; <=0 keeps all",
    )
    parser.add_argument(
        "--min-source-z-cm",
        type=float,
        default=100.0,
        help="skip rows below this z; filters metadata rows in the provided test.root",
    )
    parser.add_argument("--event-id", type=int, default=0)
    parser.add_argument("--telescope-id", type=int, default=0)
    parser.add_argument("--metadata-output", default=None)
    args = parser.parse_args()
    return args


def telescope_axes(az_deg: float, zenith_deg: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    az = math.radians(az_deg)
    el = math.radians(90.0 - zenith_deg)

    # CORSIKA/IACT NWU convention used in docs/coordinate_systems.md:
    # +x north, +y west, +z up; az=90 deg points east, i.e. -y.
    z_axis = np.array([
        math.cos(el) * math.cos(az),
        -math.cos(el) * math.sin(az),
        math.sin(el),
    ])
    z_axis = z_axis / np.linalg.norm(z_axis)
    x_axis = np.array([-math.sin(az), -math.cos(az), 0.0])
    x_axis = x_axis / np.linalg.norm(x_axis)
    y_axis = np.cross(z_axis, x_axis)
    y_axis = y_axis / np.linalg.norm(y_axis)
    return x_axis, y_axis, z_axis


def read_arrays(args: argparse.Namespace) -> dict[str, np.ndarray]:
    tree = uproot.open(args.input_root)[args.tree]
    columns = [
        args.x_column,
        args.y_column,
        args.z_column,
        args.u_column,
        args.v_column,
        args.weight_column,
        args.wavelength_column,
    ]
    arrays = tree.arrays(
        columns,
        library="np",
        entry_start=args.entry_start,
        entry_stop=args.entry_stop,
    )
    return {key: arrays[key].astype(float) for key in columns}


def main() -> None:
    args = parse_args()
    arrays = read_arrays(args)

    x = arrays[args.x_column]
    y = arrays[args.y_column]
    z = arrays[args.z_column]
    u = arrays[args.u_column]
    v = arrays[args.v_column]
    weight = arrays[args.weight_column]
    wavelength = arrays[args.wavelength_column]

    finite = (
        np.isfinite(x)
        & np.isfinite(y)
        & np.isfinite(z)
        & np.isfinite(u)
        & np.isfinite(v)
        & np.isfinite(weight)
        & np.isfinite(wavelength)
    )
    direction_ok = (u * u + v * v) <= 1.0
    physical = (z >= args.min_source_z_cm) & (weight > 0.0) & (wavelength > 0.0)
    valid = finite & direction_ok & physical

    x = x[valid]
    y = y[valid]
    z = z[valid]
    u = u[valid]
    v = v[valid]
    weight = weight[valid]
    wavelength = wavelength[valid]

    w = -np.sqrt(np.maximum(0.0, 1.0 - u * u - v * v))
    if args.coordinate_mode == "subtract_telescope_only":
        local_pos_m = np.stack([
            (x - args.telescope_x_cm) * args.position_scale,
            (y - args.telescope_y_cm) * args.position_scale,
            (z - args.telescope_z_cm) * args.position_scale,
        ], axis=1)
        local_dir = np.stack([u, v, w], axis=1)
        local_dir /= np.linalg.norm(local_dir, axis=1)[:, None]
    else:
        global_pos_cm = np.stack([x, y, z], axis=1)
        global_dir = np.stack([u, v, w], axis=1)

        x_axis, y_axis, z_axis = telescope_axes(args.telescope_az_deg, args.telescope_zenith_deg)
        local_rotation = np.stack([x_axis, y_axis, z_axis], axis=0)
        telescope_cm = np.array([
            args.telescope_x_cm,
            args.telescope_y_cm,
            args.telescope_z_cm,
        ])

        local_pos_m = (global_pos_cm - telescope_cm) @ local_rotation.T * args.position_scale
        local_dir = global_dir @ local_rotation.T
        local_dir /= np.linalg.norm(local_dir, axis=1)[:, None]

    dz = local_dir[:, 2]
    projectable = np.abs(dz) > 1e-12
    local_pos_m = local_pos_m[projectable]
    local_dir = local_dir[projectable]
    weight = weight[projectable]
    wavelength = wavelength[projectable]

    t = -local_pos_m[:, 2] / local_dir[:, 2]
    entrance_x = local_pos_m[:, 0] + local_dir[:, 0] * t
    entrance_y = local_pos_m[:, 1] + local_dir[:, 1] * t

    if args.center_mode == "mean":
        center_x = float(np.mean(entrance_x))
        center_y = float(np.mean(entrance_y))
    elif args.center_mode == "median":
        center_x = float(np.median(entrance_x))
        center_y = float(np.median(entrance_y))
    else:
        center_x = 0.0
        center_y = 0.0

    source_x = entrance_x - center_x
    source_y = entrance_y - center_y
    keep = np.ones(source_x.shape, dtype=bool)
    if args.max_radius_m > 0.0:
        keep = np.hypot(source_x, source_y) <= args.max_radius_m

    out_path = Path(args.output_csv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    kept = int(np.count_nonzero(keep))
    with out_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS)
        writer.writeheader()
        for sx, sy, direction, wl, mult in zip(
            source_x[keep],
            source_y[keep],
            local_dir[keep],
            wavelength[keep],
            weight[keep],
        ):
            writer.writerow({
                "x_m": f"{sx:.9g}",
                "y_m": f"{sy:.9g}",
                "z_m": "0",
                "dir_x": f"{direction[0]:.9g}",
                "dir_y": f"{direction[1]:.9g}",
                "dir_z": f"{direction[2]:.9g}",
                "time_ns": "0",
                "wavelength_nm": f"{wl:.9g}",
                "weight": "1",
                "multiplicity": f"{mult:.9g}",
                "event_id": str(args.event_id),
                "telescope_id": str(args.telescope_id),
            })

    if kept > 0:
        mean_dir = (
            float(np.mean(local_dir[keep, 0])),
            float(np.mean(local_dir[keep, 1])),
            float(np.mean(local_dir[keep, 2])),
        )
    else:
        mean_dir = (float("nan"), float("nan"), float("nan"))

    lines = [
        "ROOT photons -> LACT_sim PhotonCsv",
        "===================================",
        f"input_root: {args.input_root}",
        f"tree: {args.tree}",
        f"input_rows_read: {len(arrays[args.x_column])}",
        f"valid_physical_rows: {len(x)}",
        f"projectable_rows: {len(source_x)}",
        f"kept_rows: {kept}",
        f"kept_weight_sum: {float(np.sum(weight[keep])):.9g}",
        f"coordinate_mode: {args.coordinate_mode}",
        f"center_mode: {args.center_mode}",
        f"projected_center_m: ({center_x:.9g}, {center_y:.9g})",
        f"max_radius_m: {args.max_radius_m:g}",
        f"mean_local_direction_kept: ({mean_dir[0]:.9g}, "
        f"{mean_dir[1]:.9g}, {mean_dir[2]:.9g})",
        f"output_csv: {out_path}",
        "",
        "Note: center_mode=median recenters the photon bundle on the telescope entrance plane.",
        "This is intended for the standalone optical/camera response test when the ROOT file",
        "does not carry the full array/core metadata needed to place the telescope absolutely.",
    ]
    if args.metadata_output:
        meta_path = Path(args.metadata_output)
        meta_path.parent.mkdir(parents=True, exist_ok=True)
        meta_path.write_text("\n".join(lines) + "\n")

    for line in lines:
        print(line)


if __name__ == "__main__":
    main()
