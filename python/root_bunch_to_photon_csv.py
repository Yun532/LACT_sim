#!/usr/bin/env python3
"""Convert a ROOT photon-bunch tree into LACT_sim PhotonCsv."""

import argparse
import csv
import math
from pathlib import Path

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


def downward_dir_z(cx, cy):
    return -math.sqrt(max(0.0, 1.0 - cx * cx - cy * cy))


def main():
    parser = argparse.ArgumentParser(
        description="Convert ROOT photon bunch tree to LACT_sim PhotonCsv."
    )
    parser.add_argument("input_root")
    parser.add_argument("output_csv")
    parser.add_argument("--tree", default="bunch")
    parser.add_argument("--entry-start", type=int, default=None)
    parser.add_argument("--entry-stop", type=int, default=None)
    parser.add_argument("--x-column", default="bunch_x")
    parser.add_argument("--y-column", default="bunch_y")
    parser.add_argument("--cx-column", default="cx")
    parser.add_argument("--cy-column", default="cy")
    parser.add_argument("--time-column", default="time")
    parser.add_argument("--wavelength-column", default="lambda")
    parser.add_argument("--multiplicity-column", default="nbunch")
    parser.add_argument("--event-column", default="runid")
    parser.add_argument("--telescope-column", default="itel")
    parser.add_argument(
        "--position-scale",
        type=float,
        default=1.0,
        help="multiply ROOT x/y values by this factor; default assumes meters",
    )
    parser.add_argument(
        "--replace-zero-wavelength-nm",
        type=float,
        default=None,
        help="replace wavelength<=0 with this value; omit to preserve ROOT values",
    )
    args = parser.parse_args()

    tree = uproot.open(args.input_root)[args.tree]
    columns = [
        args.x_column,
        args.y_column,
        args.cx_column,
        args.cy_column,
        args.time_column,
        args.wavelength_column,
        args.multiplicity_column,
        args.event_column,
        args.telescope_column,
    ]
    arrays = tree.arrays(columns, library="np", entry_start=args.entry_start, entry_stop=args.entry_stop)

    out_path = Path(args.output_csv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    n_rows = len(arrays[args.x_column])
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=OUTPUT_COLUMNS)
        writer.writeheader()
        for i in range(n_rows):
            cx = float(arrays[args.cx_column][i])
            cy = float(arrays[args.cy_column][i])
            wavelength = float(arrays[args.wavelength_column][i])
            if args.replace_zero_wavelength_nm is not None and wavelength <= 0.0:
                wavelength = args.replace_zero_wavelength_nm
            writer.writerow({
                "x_m": float(arrays[args.x_column][i]) * args.position_scale,
                "y_m": float(arrays[args.y_column][i]) * args.position_scale,
                "z_m": 0.0,
                "dir_x": cx,
                "dir_y": cy,
                "dir_z": downward_dir_z(cx, cy),
                "time_ns": float(arrays[args.time_column][i]),
                "wavelength_nm": wavelength,
                "weight": 1.0,
                "multiplicity": float(arrays[args.multiplicity_column][i]),
                "event_id": int(arrays[args.event_column][i]),
                "telescope_id": int(arrays[args.telescope_column][i]),
            })

    print(f"wrote {out_path} with {n_rows} photon bunches")
    print(f"tree entries = {tree.num_entries}")


if __name__ == "__main__":
    main()
