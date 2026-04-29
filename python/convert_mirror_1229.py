import csv
import math
import sys


def norm(v):
    return math.sqrt(sum(x * x for x in v))


def main():
    if len(sys.argv) not in (3, 4):
        print("usage: convert_mirror_1229.py <input.csv> <output.csv> [--swap-xy]")
        return 2

    input_path = sys.argv[1]
    output_path = sys.argv[2]
    swap_xy = len(sys.argv) == 4 and sys.argv[3] == "--swap-xy"
    if len(sys.argv) == 4 and not swap_xy:
        print(f"unknown option: {sys.argv[3]}")
        return 2

    with open(input_path, newline="") as fin, open(output_path, "w", newline="") as fout:
        reader = csv.reader(fin)
        writer = csv.writer(fout)
        writer.writerow([
            "id",
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
            "reflectivity_scale",
            "roughness_sigma_rad",
            "misalign_sigma_rad",
        ])

        for out_id, row in enumerate(reader):
            if not row:
                continue
            values = [float(x) for x in row]
            if len(values) != 8:
                raise RuntimeError(f"expected 8 columns, got {len(values)}")

            center_mm = values[2:5]
            curvature_center_mm = values[5:8]
            if swap_xy:
                center_mm = [center_mm[1], center_mm[0], center_mm[2]]
                curvature_center_mm = [
                    curvature_center_mm[1],
                    curvature_center_mm[0],
                    curvature_center_mm[2],
                ]
            normal = [
                curvature_center_mm[i] - center_mm[i]
                for i in range(3)
            ]
            radius_mm = norm(normal)
            if radius_mm <= 0.0:
                raise RuntimeError(f"row {out_id}: zero curvature radius")
            normal = [x / radius_mm for x in normal]

            writer.writerow([
                out_id,
                center_mm[0] * 0.001,
                center_mm[1] * 0.001,
                center_mm[2] * 0.001,
                normal[0],
                normal[1],
                normal[2],
                "Spherical",
                16.0,
                "Hexagon",
                0.8,
                0.0,
                0.0,
                1.0,
                0.0,
                0.0,
            ])

    print(f"wrote {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
