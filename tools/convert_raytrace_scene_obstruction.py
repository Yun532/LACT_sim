#!/usr/bin/env python3
"""Convert an external raytrace_scene package into LACT obstruction primitives."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("package_dir", type=Path, help="directory containing raytrace_scene_final.json")
    parser.add_argument("--output-csv", type=Path, required=True, help="LACT obstruction primitive CSV")
    parser.add_argument("--output-cfg", type=Path, required=True, help="LACT obstruction config")
    parser.add_argument(
        "--include-camera-body",
        action="store_true",
        help="include camera_body as a solid obstruction; usually false for focal-plane tests",
    )
    parser.add_argument(
        "--unit-scale",
        type=float,
        default=0.001,
        help="source model unit to meter scale; default converts mm to m",
    )
    parser.add_argument(
        "--local-lens-z-m",
        type=float,
        default=-16.0,
        help="LACT local z coordinate of the lens/mirror reference plane",
    )
    return parser.parse_args()


def read_json(path: Path) -> list[dict]:
    with path.open(encoding="utf-8") as handle:
        data = json.load(handle)
    return list(data["objects"])


def read_lens_reference(package_dir: Path) -> tuple[float, float, float]:
    path = package_dir / "final_primitives.csv"
    if not path.exists():
        return 0.0, 0.0, -7577.2917669755
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row.get("role") == "lens_reference":
                return (
                    float(row["center_x"]),
                    float(row["center_y"]),
                    float(row["center_z"]),
                )
    return 0.0, 0.0, -7577.2917669755


def num(value, fallback=0.0) -> float:
    if value in ("", None):
        return fallback
    return float(value)


def make_transform(lens_center, unit_scale: float, local_lens_z_m: float):
    lx, ly, lz = lens_center

    def point(x, y, z):
        return (
            (float(x) - lx) * unit_scale,
            (float(y) - ly) * unit_scale,
            local_lens_z_m + (float(z) - lz) * unit_scale,
        )

    def length(v):
        return float(v) * unit_scale

    return point, length


def add_box(rows, name, center, half):
    rows.append(
        {
            "type": "box",
            "name": name,
            "x0_m": center[0],
            "y0_m": center[1],
            "z0_m": center[2],
            "x1_m": center[0],
            "y1_m": center[1],
            "z1_m": center[2],
            "radius_m": 0.0,
            "half_x_m": half[0],
            "half_y_m": half[1],
            "half_z_m": half[2],
        }
    )


def convert_scene(objects: list[dict], point, length, include_camera_body: bool) -> tuple[list[dict], dict]:
    rows: list[dict] = []
    stats = {
        "support_cylinders": 0,
        "gap_box_ring_boxes": 0,
        "camera_body_included": 0,
        "camera_body_skipped": 0,
        "void_prisms_skipped": 0,
    }

    hole = next((o for o in objects if o.get("role") == "camera_adapter_hole"), None)
    hole_half_x = length(abs(num(hole.get("bbox_max_x"))) if hole else 0.0)
    hole_half_y = length(abs(num(hole.get("bbox_max_y"))) if hole else 0.0)

    for obj in objects:
        if not obj.get("enabled", True):
            continue
        role = obj.get("role", "")
        material = obj.get("material_id", "default")
        typ = obj.get("type", "")

        if material == "void" or role == "camera_adapter_hole":
            stats["void_prisms_skipped"] += 1
            continue

        if typ == "support_cylinder":
            p0 = point(obj["p0_x"], obj["p0_y"], obj["p0_z"])
            p1 = point(obj["p1_x"], obj["p1_y"], obj["p1_z"])
            rows.append(
                {
                    "type": "cylinder",
                    "name": f"{role}_{obj['object_id']}",
                    "x0_m": p0[0],
                    "y0_m": p0[1],
                    "z0_m": p0[2],
                    "x1_m": p1[0],
                    "y1_m": p1[1],
                    "z1_m": p1[2],
                    "radius_m": length(obj["radius"]),
                    "half_x_m": 0.0,
                    "half_y_m": 0.0,
                    "half_z_m": 0.0,
                }
            )
            stats["support_cylinders"] += 1
            continue

        if typ == "box" and role == "camera_support_gap_box":
            xmin = length(obj["bbox_min_x"])
            xmax = length(obj["bbox_max_x"])
            ymin = length(obj["bbox_min_y"])
            ymax = length(obj["bbox_max_y"])
            zmin = point(0.0, 0.0, obj["bbox_min_z"])[2]
            zmax = point(0.0, 0.0, obj["bbox_max_z"])[2]
            outer_x = max(abs(xmin), abs(xmax))
            outer_y = max(abs(ymin), abs(ymax))
            inner_x = min(hole_half_x, outer_x)
            inner_y = min(hole_half_y, outer_y)
            zc = 0.5 * (zmin + zmax)
            hz = 0.5 * (zmax - zmin)

            # Four bars around the central adapter aperture. This is a
            # conservative portable approximation of the box-minus-octagon CSG.
            hx_side = 0.5 * (outer_x - inner_x)
            hy_side = outer_y
            x_side = 0.5 * (outer_x + inner_x)
            add_box(rows, f"{role}_right", (x_side, 0.0, zc), (hx_side, hy_side, hz))
            add_box(rows, f"{role}_left", (-x_side, 0.0, zc), (hx_side, hy_side, hz))
            hx_top = inner_x
            hy_top = 0.5 * (outer_y - inner_y)
            y_top = 0.5 * (outer_y + inner_y)
            add_box(rows, f"{role}_top", (0.0, y_top, zc), (hx_top, hy_top, hz))
            add_box(rows, f"{role}_bottom", (0.0, -y_top, zc), (hx_top, hy_top, hz))
            stats["gap_box_ring_boxes"] += 4
            continue

        if typ == "polygon_prism" and role == "camera_body":
            if not include_camera_body:
                stats["camera_body_skipped"] += 1
                continue
            center = point(obj["center_x"], obj["center_y"], obj["center_z"])
            half = (
                length(abs(num(obj["bbox_max_x"]))),
                length(abs(num(obj["bbox_max_y"]))),
                0.5 * length(obj["height"]),
            )
            add_box(rows, "camera_body_aabb", center, half)
            stats["camera_body_included"] += 1

    return rows, stats


def write_csv(path: Path, rows: list[dict], source: Path, stats: dict, lens_center) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        handle.write("# format=LACT_obstruction_primitives_v1\n")
        handle.write(f"# source_package={source}\n")
        handle.write("# source_units=mm\n")
        handle.write("# target_units=m\n")
        handle.write("# coordinate_transform=x/y centered on lens_reference, z(lens_reference)=-16 m\n")
        handle.write(f"# lens_reference_center={lens_center[0]},{lens_center[1]},{lens_center[2]}\n")
        for key, value in stats.items():
            handle.write(f"# {key}={value}\n")
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "type",
                "name",
                "x0_m",
                "y0_m",
                "z0_m",
                "x1_m",
                "y1_m",
                "z1_m",
                "radius_m",
                "half_x_m",
                "half_y_m",
                "half_z_m",
            ],
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_cfg(path: Path, csv_path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rel_csv = csv_path.name if csv_path.parent == path.parent else str(csv_path)
    path.write_text(
        "\n".join(
            [
                "# External 3D obstruction model converted from raytrace_scene_final.json.",
                "# The model is evaluated in LACT telescope-local coordinates.",
                "enabled=true",
                "mode=primitives",
                f"primitives_csv={rel_csv}",
                "check_incoming=true",
                "check_reflected=true",
                "",
            ]
        ),
        encoding="utf-8",
    )


def main() -> None:
    args = parse_args()
    scene_path = args.package_dir / "raytrace_scene_final.json"
    objects = read_json(scene_path)
    lens_center = read_lens_reference(args.package_dir)
    point, length = make_transform(lens_center, args.unit_scale, args.local_lens_z_m)
    rows, stats = convert_scene(objects, point, length, args.include_camera_body)
    write_csv(args.output_csv, rows, args.package_dir, stats, lens_center)
    write_cfg(args.output_cfg, args.output_csv)
    print(f"wrote {args.output_csv}")
    print(f"wrote {args.output_cfg}")
    print(f"primitives={len(rows)}")
    for key, value in stats.items():
        print(f"{key}={value}")


if __name__ == "__main__":
    main()
