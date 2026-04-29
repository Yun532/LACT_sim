#!/usr/bin/env python3
"""Convert LACT camera mapping with phy_x/phy_y in cm to simulator camera CSV."""
import argparse
import csv
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(
        description="Convert camera mapping CSV to id,x_m,y_m,shape,size_m format."
    )
    parser.add_argument("input_csv", help="input mapping CSV with s_number, phy_x, phy_y")
    parser.add_argument("output_csv", help="output camera pixel CSV")
    parser.add_argument("--pixel-size-cm", type=float, default=2.44)
    parser.add_argument("--shape", default="Square")
    parser.add_argument(
        "--id-column",
        default="s_number",
        help="column used as pixel id; s_number keeps the physical sensor numbering",
    )
    args = parser.parse_args()

    out_path = Path(args.output_csv)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    with open(args.input_csv, newline="") as fin:
        reader = csv.DictReader(fin)
        required = {args.id_column, "phy_x", "phy_y"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"missing required columns: {sorted(missing)}")
        for row in reader:
            rows.append({
                "id": int(row[args.id_column]),
                "x_m": float(row["phy_x"]) * 0.01,
                "y_m": float(row["phy_y"]) * 0.01,
                "shape": args.shape,
                "size_m": args.pixel_size_cm * 0.01,
            })

    rows.sort(key=lambda r: r["id"])
    with open(out_path, "w", newline="") as fout:
        writer = csv.DictWriter(fout, fieldnames=["id", "x_m", "y_m", "shape", "size_m"])
        writer.writeheader()
        writer.writerows(rows)

    print(f"wrote {out_path} with {len(rows)} pixels")


if __name__ == "__main__":
    main()
