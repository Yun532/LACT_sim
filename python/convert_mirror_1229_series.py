#!/usr/bin/env python3
import csv
import math
import sys
from pathlib import Path


def normalize(v):
    n = math.sqrt(sum(x * x for x in v))
    if n <= 0.0:
        raise RuntimeError("zero-length normal")
    return [x / n for x in v]


def convert_one(path, elevation_deg, swap_xy=True):
    rows = []
    with open(path, newline="") as f:
        reader = csv.reader(f)
        for facet_id, row in enumerate(reader):
            if not row:
                continue
            values = [float(x) for x in row]
            if len(values) != 8:
                raise RuntimeError(f"{path}: expected 8 columns, got {len(values)}")
            center_mm = values[2:5]
            curvature_center_mm = values[5:8]
            if swap_xy:
                center_mm = [center_mm[1], center_mm[0], center_mm[2]]
                curvature_center_mm = [curvature_center_mm[1], curvature_center_mm[0], curvature_center_mm[2]]
            center_m = [x * 0.001 for x in center_mm]
            curvature_center_m = [x * 0.001 for x in curvature_center_mm]
            normal = normalize([
                curvature_center_m[i] - center_m[i]
                for i in range(3)
            ])
            rows.append({
                "elevation_deg": elevation_deg,
                "id": facet_id,
                "cell_index": int(round(values[0])),
                "ring_index": int(round(values[1])),
                "center_x": center_m[0],
                "center_y": center_m[1],
                "center_z": center_m[2],
                "normal_x": normal[0],
                "normal_y": normal[1],
                "normal_z": normal[2],
                "surface_type": "Spherical",
                "radius_of_curvature": 16.0,
                "aperture_shape": "Hexagon",
                "size1": 0.8,
                "size2": 0.0,
                "aperture_rotation_rad": 0.0,
            })
    return rows


def main():
    if len(sys.argv) != 3:
        print("usage: convert_mirror_1229_series.py <input_dir> <output_csv>")
        return 2

    input_dir = Path(sys.argv[1])
    output_csv = Path(sys.argv[2])
    angles = list(range(0, 91, 10))

    all_rows = []
    for angle in angles:
        path = input_dir / f"mirror_1229_{angle}_base10.csv"
        all_rows.extend(convert_one(path, angle, swap_xy=True))

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(output_csv, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "elevation_deg",
                "id",
                "cell_index",
                "ring_index",
                "center_x",
                "center_y",
                "center_z",
                "normal_x",
                "normal_y",
                "normal_z",
                "surface_type",
                "radius_of_curvature",
                "aperture_shape",
                "size1",
                "size2",
                "aperture_rotation_rad",
            ],
        )
        writer.writeheader()
        writer.writerows(all_rows)

    print(f"wrote {output_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
