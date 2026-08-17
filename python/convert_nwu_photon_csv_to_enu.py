#!/usr/bin/env python3
"""Convert PhotonCsv vectors from CORSIKA NWU to east-start ENU."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


REQUIRED = ("x_m", "y_m", "z_m", "dir_x", "dir_y", "dir_z")


def converted(value: float) -> str:
    return format(value, ".17g")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Convert +x North, +y West, +z Up PhotonCsv rows to "
            "+x East, +y North, +z Up."
        )
    )
    parser.add_argument("input_csv", type=Path)
    parser.add_argument("output_csv", type=Path)
    parser.add_argument(
        "--nwu-pointing-az-deg",
        type=float,
        help="print the matching ENU east-start telescope pointing azimuth",
    )
    args = parser.parse_args()

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    rows = 0
    with args.input_csv.open(newline="", encoding="utf-8") as src:
        reader = csv.DictReader(src)
        if reader.fieldnames is None:
            raise SystemExit("input CSV has no header")
        missing = [name for name in REQUIRED if name not in reader.fieldnames]
        if missing:
            raise SystemExit(f"missing required columns: {', '.join(missing)}")

        with args.output_csv.open("w", newline="", encoding="utf-8") as dst:
            writer = csv.DictWriter(dst, fieldnames=reader.fieldnames)
            writer.writeheader()
            for row in reader:
                x_north = float(row["x_m"])
                y_west = float(row["y_m"])
                dir_north = float(row["dir_x"])
                dir_west = float(row["dir_y"])
                row["x_m"] = converted(-y_west)
                row["y_m"] = converted(x_north)
                row["dir_x"] = converted(-dir_west)
                row["dir_y"] = converted(dir_north)
                writer.writerow(row)
                rows += 1

    print(f"rows={rows}")
    print(f"output={args.output_csv}")
    if args.nwu_pointing_az_deg is not None:
        enu_az = (90.0 - args.nwu_pointing_az_deg) % 360.0
        print(f"enu_east_start_pointing_az_deg={enu_az:.12g}")


if __name__ == "__main__":
    main()
