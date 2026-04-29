#!/usr/bin/env python3
"""Normalize an external photon table into LACT_sim PhotonCsv format.

This is a lightweight adapter layer for early CORSIKA-style workflows: keep the
main C++ ray tracer fixed, and translate external columns into the common photon
CSV stream.
"""
import argparse
import csv
import math
from pathlib import Path


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


def parse_mapping(text):
    mapping = {}
    if not text:
        return mapping
    for item in text.split(","):
        item = item.strip()
        if not item:
            continue
        if "=" not in item:
            raise SystemExit(f"bad mapping item '{item}', expected output=input")
        out_col, in_col = item.split("=", 1)
        out_col = out_col.strip()
        in_col = in_col.strip()
        if out_col not in OUTPUT_COLUMNS:
            raise SystemExit(f"unsupported output column in mapping: {out_col}")
        mapping[out_col] = in_col
    return mapping


def parse_defaults(text):
    defaults = {
        "time_ns": "0",
        "wavelength_nm": "400",
        "weight": "1",
        "multiplicity": "1",
        "event_id": "0",
        "telescope_id": "0",
    }
    if not text:
        return defaults
    for item in text.split(","):
        item = item.strip()
        if not item:
            continue
        if "=" not in item:
            raise SystemExit(f"bad default item '{item}', expected column=value")
        key, value = item.split("=", 1)
        key = key.strip()
        if key not in OUTPUT_COLUMNS:
            raise SystemExit(f"unsupported default column: {key}")
        defaults[key] = value.strip()
    return defaults


def main():
    parser = argparse.ArgumentParser(
        description="Convert an external photon CSV into LACT_sim PhotonCsv format."
    )
    parser.add_argument("input_csv")
    parser.add_argument("output_csv")
    parser.add_argument(
        "--map",
        default="",
        help=(
            "comma-separated output=input mapping, e.g. "
            "x_m=x,y_m=y,z_m=z,dir_x=ux,dir_y=uy,dir_z=uz"
        ),
    )
    parser.add_argument(
        "--defaults",
        default="",
        help="comma-separated defaults for optional columns, e.g. wavelength_nm=420,weight=1",
    )
    parser.add_argument(
        "--scale-position",
        type=float,
        default=1.0,
        help="multiply mapped x_m/y_m/z_m values by this factor",
    )
    parser.add_argument(
        "--scale-time",
        type=float,
        default=1.0,
        help="multiply mapped time_ns values by this factor",
    )
    parser.add_argument(
        "--normalize-direction",
        action="store_true",
        help="normalize dir_x/dir_y/dir_z before writing",
    )
    parser.add_argument(
        "--fail-on-nonfinite",
        action="store_true",
        help="reject rows with non-finite numeric values",
    )
    args = parser.parse_args()

    mapping = parse_mapping(args.map)
    defaults = parse_defaults(args.defaults)
    required = ["x_m", "y_m", "z_m", "dir_x", "dir_y", "dir_z"]
    missing_required = [col for col in required if col not in mapping]
    if missing_required:
        raise SystemExit(f"missing required mappings: {missing_required}")

    out_path = Path(args.output_csv)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    n_rows = 0
    n_bad = 0
    with open(args.input_csv, newline="") as fin, open(out_path, "w", newline="") as fout:
        reader = csv.DictReader(fin)
        if reader.fieldnames is None:
            raise SystemExit("input CSV has no header")
        missing_inputs = sorted({mapping[col] for col in mapping}.difference(reader.fieldnames))
        if missing_inputs:
            raise SystemExit(f"input CSV missing mapped columns: {missing_inputs}")
        writer = csv.DictWriter(fout, fieldnames=OUTPUT_COLUMNS)
        writer.writeheader()
        for line_no, row in enumerate(reader, start=2):
            out = {}
            for col in OUTPUT_COLUMNS:
                if col in mapping:
                    value = row[mapping[col]]
                else:
                    value = defaults.get(col, "")
                if col in {"x_m", "y_m", "z_m"}:
                    value = str(float(value) * args.scale_position)
                elif col == "time_ns":
                    value = str(float(value) * args.scale_time)
                out[col] = value

            numeric_cols = [
                "x_m", "y_m", "z_m", "dir_x", "dir_y", "dir_z", "time_ns",
                "wavelength_nm", "weight", "multiplicity",
            ]
            try:
                numbers = {col: float(out[col]) for col in numeric_cols}
                event_id = int(float(out["event_id"]))
                telescope_id = int(float(out["telescope_id"]))
            except ValueError as exc:
                raise SystemExit(f"line {line_no}: invalid numeric value: {exc}") from exc

            if args.fail_on_nonfinite:
                if not all(math.isfinite(v) for v in numbers.values()):
                    raise SystemExit(f"line {line_no}: non-finite numeric value")

            if args.normalize_direction:
                norm = math.sqrt(
                    numbers["dir_x"] ** 2 + numbers["dir_y"] ** 2 + numbers["dir_z"] ** 2
                )
                if norm <= 0.0 or not math.isfinite(norm):
                    raise SystemExit(f"line {line_no}: invalid direction norm {norm}")
                out["dir_x"] = f"{numbers['dir_x'] / norm:.12g}"
                out["dir_y"] = f"{numbers['dir_y'] / norm:.12g}"
                out["dir_z"] = f"{numbers['dir_z'] / norm:.12g}"

            out["event_id"] = str(event_id)
            out["telescope_id"] = str(telescope_id)
            writer.writerow(out)
            n_rows += 1

    print(f"wrote {out_path} with {n_rows} photons")
    if n_bad:
        print(f"skipped {n_bad} bad photons")


if __name__ == "__main__":
    main()
