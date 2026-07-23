#!/usr/bin/env python3
"""Make a six-column, one-row-per-ray PhotonCsv from CORSIKA NWU rows.

The input may be the detailed output of eventio_to_photon_csv.  Each selected
row is treated as one representative photon; bunch multiplicity, wavelength,
time, and event metadata are intentionally not copied.  Coordinates are kept
in CORSIKA NWU so the simulator's normal ``corsika_nwu_relative`` adapter is
the single source of coordinate conversion.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


OUTPUT_COLUMNS = ("x_m", "y_m", "z_m", "dir_x", "dir_y", "dir_z")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_csv", type=Path)
    parser.add_argument("output_csv", type=Path)
    parser.add_argument("--telescope-id", type=int)
    parser.add_argument("--max-rows", type=int)
    args = parser.parse_args()

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    with args.input_csv.open(newline="", encoding="utf-8") as src:
        reader = csv.DictReader(src)
        if reader.fieldnames is None:
            raise SystemExit("input CSV has no header")
        required = set(OUTPUT_COLUMNS)
        missing = sorted(required.difference(reader.fieldnames))
        if missing:
            raise SystemExit(f"missing required columns: {', '.join(missing)}")

        with args.output_csv.open("w", newline="", encoding="utf-8") as dst:
            writer = csv.DictWriter(dst, fieldnames=OUTPUT_COLUMNS)
            writer.writeheader()
            for row in reader:
                if args.telescope_id is not None:
                    if "telescope_id" not in row:
                        raise SystemExit("--telescope-id requires a telescope_id column")
                    if int(row["telescope_id"]) != args.telescope_id:
                        continue

                writer.writerow(
                    {
                        column: row[column] for column in OUTPUT_COLUMNS
                    }
                )
                written += 1
                if args.max_rows is not None and written >= args.max_rows:
                    break

    if written == 0:
        raise SystemExit("no rows selected")
    print(f"rows={written}")
    print(f"output={args.output_csv}")


if __name__ == "__main__":
    main()
