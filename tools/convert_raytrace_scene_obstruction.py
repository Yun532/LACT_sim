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
        default=True,
        help="include camera_body as a solid obstruction; enabled by default",
    )
    parser.add_argument(
        "--exclude-camera-body",
        action="store_false",
        dest="include_camera_body",
        help="exclude camera_body from the generated obstruction CSV",
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


def empty_row() -> dict:
    return {
        "type": "",
        "name": "",
        "role": "",
        "material_id": "default",
        "x0_m": "",
        "y0_m": "",
        "z0_m": "",
        "x1_m": "",
        "y1_m": "",
        "z1_m": "",
        "center_x_m": "",
        "center_y_m": "",
        "center_z_m": "",
        "radius_m": "",
        "height_m": "",
        "rotation_rad": "",
        "sides": "",
        "half_x_m": "",
        "half_y_m": "",
        "half_z_m": "",
        "bbox_min_x_m": "",
        "bbox_min_y_m": "",
        "bbox_min_z_m": "",
        "bbox_max_x_m": "",
        "bbox_max_y_m": "",
        "bbox_max_z_m": "",
        "hole_radius_m": "",
        "hole_rotation_rad": "",
        "hole_sides": "",
    }


def bbox_to_local(obj, point):
    bmin = point(obj["bbox_min_x"], obj["bbox_min_y"], obj["bbox_min_z"])
    bmax = point(obj["bbox_max_x"], obj["bbox_max_y"], obj["bbox_max_z"])
    return (
        (min(bmin[0], bmax[0]), min(bmin[1], bmax[1]), min(bmin[2], bmax[2])),
        (max(bmin[0], bmax[0]), max(bmin[1], bmax[1]), max(bmin[2], bmax[2])),
    )


def convert_scene(objects: list[dict], point, length, include_camera_body: bool) -> tuple[list[dict], dict]:
    rows: list[dict] = []
    stats = {
        "support_cylinders": 0,
        "gap_boxes": 0,
        "adapter_holes": 0,
        "camera_body_included": 0,
        "camera_body_skipped": 0,
    }

    hole = next((o for o in objects if o.get("role") == "camera_adapter_hole"), None)
    hole_radius = length(hole["radius"]) if hole else 0.0
    hole_rotation = float(hole.get("rotation", 0.0)) if hole else 0.0
    hole_sides = int(float(hole.get("sides", 0))) if hole else 0

    for obj in objects:
        if not obj.get("enabled", True):
            continue
        role = obj.get("role", "")
        material = obj.get("material_id", "default")
        typ = obj.get("type", "")

        if material == "void" or role == "camera_adapter_hole":
            stats["adapter_holes"] += 1
            continue

        if typ == "support_cylinder":
            p0 = point(obj["p0_x"], obj["p0_y"], obj["p0_z"])
            p1 = point(obj["p1_x"], obj["p1_y"], obj["p1_z"])
            row = empty_row()
            row.update({
                "type": "cylinder",
                "name": f"{role}_{obj['object_id']}",
                "role": role,
                "material_id": material,
                "x0_m": p0[0],
                "y0_m": p0[1],
                "z0_m": p0[2],
                "x1_m": p1[0],
                "y1_m": p1[1],
                "z1_m": p1[2],
                "radius_m": length(obj["radius"]),
            })
            rows.append(row)
            stats["support_cylinders"] += 1
            continue

        if typ == "box" and role == "camera_support_gap_box":
            bmin, bmax = bbox_to_local(obj, point)
            center = tuple(0.5 * (lo + hi) for lo, hi in zip(bmin, bmax))
            half = tuple(0.5 * (hi - lo) for lo, hi in zip(bmin, bmax))
            row = empty_row()
            row.update({
                "type": "box",
                "name": role,
                "role": role,
                "material_id": material,
                "center_x_m": center[0],
                "center_y_m": center[1],
                "center_z_m": center[2],
                "x0_m": center[0],
                "y0_m": center[1],
                "z0_m": center[2],
                "x1_m": center[0],
                "y1_m": center[1],
                "z1_m": center[2],
                "half_x_m": half[0],
                "half_y_m": half[1],
                "half_z_m": half[2],
                "bbox_min_x_m": bmin[0],
                "bbox_min_y_m": bmin[1],
                "bbox_min_z_m": bmin[2],
                "bbox_max_x_m": bmax[0],
                "bbox_max_y_m": bmax[1],
                "bbox_max_z_m": bmax[2],
                "hole_radius_m": hole_radius,
                "hole_rotation_rad": hole_rotation,
                "hole_sides": hole_sides,
            })
            rows.append(row)
            stats["gap_boxes"] += 1
            continue

        if typ == "polygon_prism" and role == "camera_body":
            if not include_camera_body:
                stats["camera_body_skipped"] += 1
                continue
            center = point(obj["center_x"], obj["center_y"], obj["center_z"])
            row = empty_row()
            bmin, bmax = bbox_to_local(obj, point)
            row.update({
                "type": "polygon_prism",
                "name": role,
                "role": role,
                "material_id": material,
                "center_x_m": center[0],
                "center_y_m": center[1],
                "center_z_m": center[2],
                "x0_m": center[0],
                "y0_m": center[1],
                "z0_m": center[2],
                "x1_m": center[0],
                "y1_m": center[1],
                "z1_m": center[2],
                "radius_m": length(obj["radius"]),
                "height_m": length(obj["height"]),
                "rotation_rad": obj["rotation"],
                "sides": int(float(obj["sides"])),
                "bbox_min_x_m": bmin[0],
                "bbox_min_y_m": bmin[1],
                "bbox_min_z_m": bmin[2],
                "bbox_max_x_m": bmax[0],
                "bbox_max_y_m": bmax[1],
                "bbox_max_z_m": bmax[2],
            })
            rows.append(row)
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
                "role",
                "material_id",
                "x0_m",
                "y0_m",
                "z0_m",
                "x1_m",
                "y1_m",
                "z1_m",
                "center_x_m",
                "center_y_m",
                "center_z_m",
                "radius_m",
                "height_m",
                "rotation_rad",
                "sides",
                "half_x_m",
                "half_y_m",
                "half_z_m",
                "bbox_min_x_m",
                "bbox_min_y_m",
                "bbox_min_z_m",
                "bbox_max_x_m",
                "bbox_max_y_m",
                "bbox_max_z_m",
                "hole_radius_m",
                "hole_rotation_rad",
                "hole_sides",
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
