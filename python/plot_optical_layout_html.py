#!/usr/bin/env python3
"""Write an interactive self-contained 3D HTML view of the LACT optics."""

from __future__ import annotations

import argparse
import csv
import json
import math
from html import escape
from pathlib import Path

import numpy as np

from config_io import (
    apply_telescope_frame_to_facets,
    expand_component_config,
    load_facets_from_config,
    parse_vec3,
    point_to_global,
    rotate_local_vector,
    source_direction_from_config,
    telescope_frame_from_config,
)


def local_frame(normal):
    n = normal / np.linalg.norm(normal)
    ref = np.array([0.0, 0.0, 1.0]) if abs(n[2]) < 0.9 else np.array([0.0, 1.0, 0.0])
    u = np.cross(ref, n)
    u = u / np.linalg.norm(u)
    v = np.cross(n, u)
    return u, v, n


def source_coordinate_frame_from_config(cfg):
    """Return the physical input-adapter basis used by OpticalSimCommon.cpp."""
    raw_name = cfg.get("source.coordinate_frame", "telescope_local").strip().lower()
    aliases = {
        "local": "telescope_local",
        "optical_local": "telescope_local",
        "corsika_iact": "corsika_nwu_relative",
        "corsika": "corsika_nwu_relative",
        "simtelarray": "corsika_nwu_relative",
        "corsika_global": "corsika_nwu_global",
        "enu_relative": "enu_east_relative",
        "east_north_up_relative": "enu_east_relative",
        "east_start_relative": "enu_east_relative",
        "enu_global": "enu_east_global",
        "east_north_up_global": "enu_east_global",
        "east_start_global": "enu_east_global",
        "generic_global": "lact_generic_global",
        "array_global": "lact_generic_global",
        "global": "lact_generic_global",
    }
    name = aliases.get(raw_name, raw_name)
    origin = parse_vec3(cfg.get("telescope.position_m"), [0.0, 0.0, 0.0])
    az = math.radians(float(cfg.get("telescope.pointing_az_deg", 0.0)))
    el = math.radians(float(cfg.get("telescope.pointing_el_deg", 90.0)))
    sin_az, cos_az = math.sin(az), math.cos(az)
    sin_el, cos_el = math.sin(el), math.cos(el)

    if name in {"corsika_nwu_relative", "corsika_nwu_global"}:
        frame = {
            "origin": origin,
            "x_axis": np.array([-sin_el * cos_az, sin_el * sin_az, cos_el]),
            "y_axis": np.array([-sin_az, -cos_az, 0.0]),
            "z_axis": np.array([cos_el * cos_az, -cos_el * sin_az, sin_el]),
        }
        input_axes = [
            (np.array([1.0, 0.0, 0.0]), "input +x / magnetic North", "#ff5c5c"),
            (np.array([0.0, 1.0, 0.0]), "input +y / West", "#54d67a"),
            (np.array([0.0, 0.0, 1.0]), "input +z / Up", "#4aa3ff"),
        ]
        description = "CORSIKA NWU: +x magnetic North, +y West, +z Up"
    elif name in {"enu_east_relative", "enu_east_global"}:
        frame = {
            "origin": origin,
            "x_axis": np.array([-sin_el * cos_az, -sin_el * sin_az, cos_el]),
            "y_axis": np.array([sin_az, -cos_az, 0.0]),
            "z_axis": np.array([cos_el * cos_az, cos_el * sin_az, sin_el]),
        }
        input_axes = [
            (np.array([1.0, 0.0, 0.0]), "input +x / East", "#ff5c5c"),
            (np.array([0.0, 1.0, 0.0]), "input +y / North", "#54d67a"),
            (np.array([0.0, 0.0, 1.0]), "input +z / Up", "#4aa3ff"),
        ]
        description = "ENU: +x East, +y North, +z Up"
    elif name == "telescope_local":
        frame = {
            "origin": origin,
            "x_axis": np.array([1.0, 0.0, 0.0]),
            "y_axis": np.array([0.0, 1.0, 0.0]),
            "z_axis": np.array([0.0, 0.0, 1.0]),
        }
        input_axes = [
            (frame["x_axis"], "input +x / telescope local x", "#ff5c5c"),
            (frame["y_axis"], "input +y / telescope local y", "#54d67a"),
            (frame["z_axis"], "input +z / boresight", "#4aa3ff"),
        ]
        description = "Input is already in telescope-local optical coordinates"
    elif name == "lact_generic_global":
        frame = telescope_frame_from_config(cfg)
        input_axes = [
            (np.array([1.0, 0.0, 0.0]), "input +x / generic global X", "#ff5c5c"),
            (np.array([0.0, 1.0, 0.0]), "input +y / generic global Y", "#54d67a"),
            (np.array([0.0, 0.0, 1.0]), "input +z / Up", "#4aa3ff"),
        ]
        description = "Legacy LACT generic global XY"
    else:
        raise ValueError(f"unsupported source.coordinate_frame: {raw_name}")

    for key in ("x_axis", "y_axis", "z_axis"):
        frame[key] = frame[key] / np.linalg.norm(frame[key])
    return frame, name, input_axes, description


def aperture_polygon(facet):
    center = facet["center"]
    u = np.array(facet.get("u_axis", local_frame(facet["normal"])[0]), dtype=float)
    v = np.array(facet.get("v_axis", local_frame(facet["normal"])[1]), dtype=float)
    u = u / np.linalg.norm(u)
    v = v / np.linalg.norm(v)
    size = facet["size1"]
    if facet["shape"] == "hexagon":
        a = 0.5 * size
        local = np.array([
            [a, a / math.sqrt(3.0)],
            [0.0, 2.0 * a / math.sqrt(3.0)],
            [-a, a / math.sqrt(3.0)],
            [-a, -a / math.sqrt(3.0)],
            [0.0, -2.0 * a / math.sqrt(3.0)],
            [a, -a / math.sqrt(3.0)],
        ])
    else:
        h = 0.5 * size
        local = np.array([[h, h], [-h, h], [-h, -h], [h, -h]])
    rot = facet.get("rotation", 0.0)
    if rot:
        c, s = math.cos(rot), math.sin(rot)
        local = local @ np.array([[c, -s], [s, c]]).T
    return [list((center + x * u + y * v).astype(float)) for x, y in local]


def regular_prism_edges(center, radius, height, rotation, sides):
    if radius <= 0 or height <= 0 or sides < 3:
        return []
    angles = rotation + np.arange(sides) * (2.0 * math.pi / sides)
    bottom = [
        [center[0] + radius * math.cos(a), center[1] + radius * math.sin(a), center[2] - 0.5 * height]
        for a in angles
    ]
    top = [[p[0], p[1], center[2] + 0.5 * height] for p in bottom]
    lines = [bottom + [bottom[0]], top + [top[0]]]
    for a, b in zip(bottom, top):
        lines.append([a, b])
    return lines


def box_edges(center, half):
    x, y, z = center
    hx, hy, hz = half
    c = [
        [x - hx, y - hy, z - hz], [x + hx, y - hy, z - hz],
        [x + hx, y + hy, z - hz], [x - hx, y + hy, z - hz],
        [x - hx, y - hy, z + hz], [x + hx, y - hy, z + hz],
        [x + hx, y + hy, z + hz], [x - hx, y + hy, z + hz],
    ]
    return [
        [c[0], c[1], c[2], c[3], c[0]],
        [c[4], c[5], c[6], c[7], c[4]],
        [c[0], c[4]], [c[1], c[5]], [c[2], c[6]], [c[3], c[7]],
    ]


def transform_line(points, frame):
    return [list(point_to_global(np.array(p, dtype=float), frame).astype(float)) for p in points]


def load_input_local_photon_lines(path, frame, stride, length_m):
    """Convert saved local input photons to world-space direction segments."""
    required = {
        "input_local_x_m", "input_local_y_m", "input_local_z_m",
        "input_local_dir_x", "input_local_dir_y", "input_local_dir_z",
    }
    lines = []
    with Path(path).open(newline="") as handle:
        rows = csv.DictReader(handle)
        available = set(rows.fieldnames or [])
        missing = sorted(required - available)
        if missing:
            raise ValueError(
                "input photon CSV is missing: " + ", ".join(missing) +
                "; enable output.whiteboard_input_photon=true when tracing"
            )
        for index, row in enumerate(rows):
            if index % stride:
                continue
            local_pos = np.array([
                float(row["input_local_x_m"]),
                float(row["input_local_y_m"]),
                float(row["input_local_z_m"]),
            ])
            local_dir = np.array([
                float(row["input_local_dir_x"]),
                float(row["input_local_dir_y"]),
                float(row["input_local_dir_z"]),
            ])
            norm = np.linalg.norm(local_dir)
            if not np.isfinite(norm) or norm == 0.0:
                continue
            start = point_to_global(local_pos, frame)
            direction = rotate_local_vector(local_dir / norm, frame)
            lines.append({
                "points": [list(start.astype(float)),
                           list((start + length_m * direction).astype(float))],
                "color": "#d98cff",
                "role": "input_photon",
                "name": "input photon (telescope local)",
                "width": 2,
            })
    return lines


def load_primitives(path: Path):
    primitives = []
    with path.open(newline="") as handle:
        rows = csv.DictReader(line for line in handle if not line.lstrip().startswith("#"))
        for row in rows:
            def f(key, default=0.0):
                value = row.get(key, "")
                return default if value in ("", None) else float(value)

            primitives.append({
                "type": row.get("type", "").lower(),
                "name": row.get("name", ""),
                "role": row.get("role", ""),
                "p0": [f("x0_m"), f("y0_m"), f("z0_m")],
                "p1": [f("x1_m"), f("y1_m"), f("z1_m")],
                "center": [f("center_x_m", f("x0_m")), f("center_y_m", f("y0_m")), f("center_z_m", f("z0_m"))],
                "radius": f("radius_m"),
                "height": f("height_m"),
                "rotation": f("rotation_rad"),
                "sides": int(f("sides")),
                "half": [f("half_x_m"), f("half_y_m"), f("half_z_m")],
                "hole_radius": f("hole_radius_m"),
                "hole_rotation": f("hole_rotation_rad"),
                "hole_sides": int(f("hole_sides")),
            })
    return primitives


def resolve_primitives_path(cfg):
    path = cfg.get("obstruction.primitives_csv")
    if not path:
        raise SystemExit("config does not define obstruction.primitives_csv")
    return Path(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--ray-stride", type=int, default=12)
    parser.add_argument(
        "--input-photon-csv",
        default="",
        help="whiteboard CSV with input_local_* columns to overlay",
    )
    parser.add_argument(
        "--input-photon-stride",
        type=int,
        default=1,
        help="draw one input-photon direction segment every N CSV rows",
    )
    parser.add_argument(
        "--input-photon-length-m",
        type=float,
        default=3.0,
        help="world-space length of each overlaid input-photon direction segment",
    )
    parser.add_argument(
        "--show-coordinate-axes",
        action="store_true",
        help="draw input, telescope-local, camera, sky, and photon directions",
    )
    parser.add_argument(
        "--coordinate-frame-mode",
        choices=("layout", "source"),
        default="layout",
        help=(
            "layout uses the generic optical-layout placement; source uses the exact "
            "source.coordinate_frame adapter basis (recommended for coordinate explanations)"
        ),
    )
    parser.add_argument(
        "--highlight-primitives",
        default="",
        help="comma-separated obstruction primitive names to highlight and label",
    )
    args = parser.parse_args()
    if args.input_photon_stride <= 0:
        raise SystemExit("--input-photon-stride must be positive")
    if args.input_photon_length_m <= 0.0:
        raise SystemExit("--input-photon-length-m must be positive")
    highlighted = {name.strip() for name in args.highlight_primitives.split(",") if name.strip()}

    cfg_path = Path(args.config).resolve()
    cfg, _ = expand_component_config(cfg_path)
    if args.coordinate_frame_mode == "source":
        frame, source_frame_name, input_axes, source_frame_description = (
            source_coordinate_frame_from_config(cfg)
        )
    else:
        frame = telescope_frame_from_config(cfg)
        source_frame_name = cfg.get("source.coordinate_frame", "telescope_local")
        input_axes = [
            (np.array([1.0, 0.0, 0.0]), "global +X", "#ff5c5c"),
            (np.array([0.0, 1.0, 0.0]), "global +Y", "#54d67a"),
            (np.array([0.0, 0.0, 1.0]), "global +Z / Up", "#4aa3ff"),
        ]
        source_frame_description = "Generic optical-layout placement"
    facets = apply_telescope_frame_to_facets(load_facets_from_config(cfg_path, cfg), frame)

    local_plane_point = parse_vec3(cfg.get("output.plane_point"), [0, 0, -8])
    local_plane_normal = parse_vec3(cfg.get("output.plane_normal"), [0, 0, -1])
    local_plane_normal = local_plane_normal / np.linalg.norm(local_plane_normal)
    plane_point = point_to_global(local_plane_point, frame)
    plane_normal = rotate_local_vector(local_plane_normal, frame)
    local_u = parse_vec3(cfg["output.plane_u_axis"]) if "output.plane_u_axis" in cfg else None
    local_v = parse_vec3(cfg["output.plane_v_axis"]) if "output.plane_v_axis" in cfg else None
    if local_u is None or local_v is None:
        u, v, _ = local_frame(plane_normal)
    else:
        u = rotate_local_vector(local_u, frame)
        v = rotate_local_vector(local_v, frame)
        u = u / np.linalg.norm(u)
        v = v / np.linalg.norm(v)
    plane_radius = 0.8
    plane = [
        list((plane_point - plane_radius * u - plane_radius * v).astype(float)),
        list((plane_point + plane_radius * u - plane_radius * v).astype(float)),
        list((plane_point + plane_radius * u + plane_radius * v).astype(float)),
        list((plane_point - plane_radius * u + plane_radius * v).astype(float)),
    ]

    polygons = [{"points": aperture_polygon(f), "color": "#75aadb", "role": "mirror"} for f in facets]
    polygons.append({"points": plane, "color": "#e05a47", "role": "output_plane"})

    lines = []
    labels = []
    if args.show_coordinate_axes:
        focal_length_m = float(cfg.get("telescope.focal_length_m", 8.0))
        local_mirror_vertex = local_plane_point + focal_length_m * local_plane_normal
        mirror_origin = point_to_global(local_mirror_vertex, frame)
        labels.append({
            "point": list(mirror_origin.astype(float)),
            "text": f"mirror vertex / input record plane z={local_mirror_vertex[2]:g} m",
            "color": "#75aadb",
        })

        def add_axis(start, direction, length, label_text, color, role="coordinate_axis"):
            direction = direction / np.linalg.norm(direction)
            end = start + length * direction
            lines.append({
                "points": [list(start.astype(float)), list(end.astype(float))],
                "color": color,
                "role": role,
                "name": label_text,
                "width": 5,
                "arrow": True,
            })
            labels.append({
                "point": list(end.astype(float)),
                "text": label_text,
                "color": color,
            })

        for direction, label_text, color in input_axes:
            add_axis(mirror_origin, direction, 3.2, label_text, color, "input_axis")

        if args.coordinate_frame_mode == "source" and source_frame_name != "telescope_local":
            local_x_label = "local +x / increasing elevation"
            local_y_label = "local +y / azimuth N -> E"
        else:
            local_x_label = "local +x"
            local_y_label = "local +y"
        local_axes = [
            (frame["x_axis"], local_x_label, "#ff9f43"),
            (frame["y_axis"], local_y_label, "#c678dd"),
            (frame["z_axis"], "local +z / boresight -> sky", "#33d6ff"),
        ]
        for direction, label_text, color in local_axes:
            add_axis(mirror_origin, direction, 4.4, label_text, color, "telescope_local_axis")

        add_axis(plane_point, u, 2.2, "camera +u = output x_m", "#ffd84d", "camera_axis")
        add_axis(plane_point, v, 2.2, "camera +v = output y_m", "#68e88b", "camera_axis")
        add_axis(
            plane_point,
            plane_normal,
            2.0,
            "camera normal = local -z",
            "#ff6b8a",
            "camera_normal",
        )

        boresight = frame["z_axis"] / np.linalg.norm(frame["z_axis"])
        sky_point = mirror_origin + 9.5 * boresight
        labels.append({
            "point": list(sky_point.astype(float)),
            "text": "SKY / telescope pointing",
            "color": "#33d6ff",
        })
        photon_start = mirror_origin + 9.5 * boresight
        photon_end = mirror_origin + 0.8 * boresight
        lines.append({
            "points": [list(photon_start.astype(float)), list(photon_end.astype(float))],
            "color": "#ffffff",
            "role": "incoming_photon_direction",
            "name": "on-axis incoming photon: -local z",
            "width": 4,
            "arrow": True,
        })
        labels.append({
            "point": list((photon_start + 0.25 * (photon_end - photon_start)).astype(float)),
            "text": "incoming photon travels along -local z",
            "color": "#ffffff",
        })
        lines.append({
            "points": [list(mirror_origin.astype(float)), list(plane_point.astype(float))],
            "color": "#ff74d4",
            "role": "reflected_ray_direction",
            "name": "reflected ray toward camera",
            "width": 3,
            "arrow": True,
        })
        labels.append({
            "point": list((mirror_origin + 0.58 * (plane_point - mirror_origin)).astype(float)),
            "text": "reflected light -> camera",
            "color": "#ff74d4",
        })
    for p in load_primitives(resolve_primitives_path(cfg)):
        role = p["role"]
        is_highlighted = p["name"] in highlighted
        color = "#ff3b30" if is_highlighted else ("#f08a24" if role == "support_strut" else "#28a95f")
        width = 6 if is_highlighted else 3
        if p["type"] == "cylinder":
            a = point_to_global(np.array(p["p0"]), frame)
            b = point_to_global(np.array(p["p1"]), frame)
            lines.append({
                "points": [list(a.astype(float)), list(b.astype(float))],
                "color": color,
                "role": role,
                "name": p["name"],
                "width": width,
            })
            if is_highlighted:
                label_pos = 0.5 * (a + b)
                labels.append({"point": list(label_pos.astype(float)), "text": p["name"], "color": color})
        elif p["type"] in {"box", "aabb"}:
            center = np.array(p["center"])
            for edge in box_edges(center, p["half"]):
                lines.append({"points": transform_line(edge, frame), "color": color, "role": role, "name": p["name"], "width": 3 if is_highlighted else 2})
            if p["hole_radius"] > 0:
                for edge in regular_prism_edges(center, p["hole_radius"], p["half"][2] * 2, p["hole_rotation"], p["hole_sides"]):
                    lines.append({"points": transform_line(edge, frame), "color": "#e02f45", "role": "camera_adapter_hole", "width": 2})
            if is_highlighted:
                labels.append({"point": list(point_to_global(center, frame).astype(float)), "text": p["name"], "color": color})

        elif p["type"] == "polygon_prism":
            center = np.array(p["center"])
            for edge in regular_prism_edges(center, p["radius"], p["height"], p["rotation"], p["sides"]):
                lines.append({"points": transform_line(edge, frame), "color": color if is_highlighted else "#00b9c7", "role": role, "name": p["name"], "width": 3 if is_highlighted else 2})
            if is_highlighted:
                labels.append({"point": list(point_to_global(center, frame).astype(float)), "text": p["name"], "color": color})

    input_photon_count = 0
    if args.input_photon_csv:
        input_lines = load_input_local_photon_lines(
            args.input_photon_csv,
            frame,
            args.input_photon_stride,
            args.input_photon_length_m,
        )
        input_photon_count = len(input_lines)
        lines.extend(input_lines)

    all_points = [pt for poly in polygons for pt in poly["points"]] + [pt for line in lines for pt in line["points"]]
    arr = np.array(all_points, dtype=float)
    bounds = [arr.min(axis=0).tolist(), arr.max(axis=0).tolist()]
    data = {"polygons": polygons, "lines": lines, "labels": labels, "bounds": bounds}

    input_photon_legend = (
        f" · Purple: {input_photon_count} input-photon directions"
        if input_photon_count else ""
    )
    pointing_az = float(cfg.get("telescope.pointing_az_deg", 0.0))
    pointing_el = float(cfg.get("telescope.pointing_el_deg", 90.0))
    coordinate_summary = (
        f"<b>Input:</b> {escape(source_frame_description)} "
        f"(<code>{escape(source_frame_name)}</code>)<br>"
        f"<b>Pointing:</b> az={pointing_az:g} deg, elevation={pointing_el:g} deg; "
        "local +z points to the sky, incoming photons travel along local -z.<br>"
        "<b>Camera/output:</b> u=local +x, v=local +y; stored "
        "<code>x_m=u</code>, <code>y_m=v</code>.<br>"
        "<b>pyLAST:</b> <code>pix_x=-u</code>, <code>pix_y=-v</code>; current camera plot "
        "uses horizontal=-v and vertical=-u."
    )
    html = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><title>LACT 3D Coordinate System</title>
<style>
html,body{{margin:0;height:100%;background:#101216;color:#e8edf2;font-family:Segoe UI,Arial,sans-serif}}
#bar{{position:fixed;z-index:2;left:12px;top:10px;max-width:min(760px,calc(100vw - 48px));padding:10px 12px;background:rgba(16,18,22,.88);border:1px solid #39414d;border-radius:8px;font-size:13px;line-height:1.45;box-shadow:0 6px 24px #0008}}
#bar .title{{font-size:16px;font-weight:700;color:#fff;margin-bottom:4px}}
#bar .controls{{color:#aeb8c5;margin-top:5px}}
#bar code{{color:#ffd84d}}
canvas{{width:100vw;height:100vh;display:block;cursor:grab}}canvas:active{{cursor:grabbing}}
</style>
</head><body><canvas id="view"></canvas><div id="bar"><div class="title">LACT physical 3D coordinate model</div>{coordinate_summary}<div class="controls">Drag to rotate · Wheel to zoom · Blue facets: mirrors · Red plane: camera/output · Orange/green/cyan: telescope structure{input_photon_legend}</div></div>
<script>
const data = {json.dumps(data)};
const canvas=document.getElementById('view'),ctx=canvas.getContext('2d');let rx=-0.7,ry=0.8,zoom=1,drag=false,lx=0,ly=0;
const b=data.bounds,center=[(b[0][0]+b[1][0])/2,(b[0][1]+b[1][1])/2,(b[0][2]+b[1][2])/2],span=Math.max(b[1][0]-b[0][0],b[1][1]-b[0][1],b[1][2]-b[0][2])||1;
function resize(){{canvas.width=innerWidth*devicePixelRatio;canvas.height=innerHeight*devicePixelRatio;draw();}}
function rot(p){{let x=(p[0]-center[0])/span,y=(p[1]-center[1])/span,z=(p[2]-center[2])/span;const cy=Math.cos(ry),sy=Math.sin(ry),cx=Math.cos(rx),sx=Math.sin(rx);let x1=x*cy+z*sy,z1=-x*sy+z*cy;let y1=y*cx-z1*sx,z2=y*sx+z1*cx;return[x1,y1,z2];}}
function proj(p){{const r=rot(p),s=Math.min(canvas.width,canvas.height)*0.78*zoom;return[canvas.width/2+r[0]*s,canvas.height/2-r[1]*s,r[2]];}}
function line(points,color,width=1,closed=false,arrow=false){{if(!points.length)return;ctx.beginPath();let p=proj(points[0]);ctx.moveTo(p[0],p[1]);for(let i=1;i<points.length;i++){{p=proj(points[i]);ctx.lineTo(p[0],p[1]);}}if(closed)ctx.closePath();ctx.strokeStyle=color;ctx.lineWidth=width*devicePixelRatio;ctx.stroke();if(arrow&&points.length>1){{const a=proj(points[points.length-2]),b=proj(points[points.length-1]),angle=Math.atan2(b[1]-a[1],b[0]-a[0]),size=(7+width)*devicePixelRatio;ctx.beginPath();ctx.moveTo(b[0],b[1]);ctx.lineTo(b[0]-size*Math.cos(angle-.48),b[1]-size*Math.sin(angle-.48));ctx.lineTo(b[0]-size*Math.cos(angle+.48),b[1]-size*Math.sin(angle+.48));ctx.closePath();ctx.fillStyle=color;ctx.fill();}}}}
function poly(points,color){{ctx.beginPath();let p=proj(points[0]);ctx.moveTo(p[0],p[1]);for(let i=1;i<points.length;i++){{p=proj(points[i]);ctx.lineTo(p[0],p[1]);}}ctx.closePath();ctx.fillStyle=color+'44';ctx.strokeStyle=color;ctx.lineWidth=devicePixelRatio;ctx.fill();ctx.stroke();}}
function label(item){{const p=proj(item.point);ctx.font=`${{13*devicePixelRatio}}px Segoe UI,Arial,sans-serif`;ctx.textBaseline='middle';const pad=4*devicePixelRatio,w=ctx.measureText(item.text).width+2*pad,h=20*devicePixelRatio;ctx.fillStyle='rgba(16,18,22,.78)';ctx.strokeStyle=item.color;ctx.lineWidth=1.5*devicePixelRatio;ctx.fillRect(p[0]+6*devicePixelRatio,p[1]-h/2,w,h);ctx.strokeRect(p[0]+6*devicePixelRatio,p[1]-h/2,w,h);ctx.fillStyle='#fff';ctx.fillText(item.text,p[0]+6*devicePixelRatio+pad,p[1]);}}
function draw(){{ctx.fillStyle='#101216';ctx.fillRect(0,0,canvas.width,canvas.height);for(const p of data.polygons)poly(p.points,p.color);for(const l of data.lines)line(l.points,l.color,l.width||1,false,!!l.arrow);for(const item of data.labels)label(item);}}
canvas.addEventListener('mousedown',e=>{{drag=true;lx=e.clientX;ly=e.clientY;}});addEventListener('mouseup',()=>drag=false);addEventListener('mousemove',e=>{{if(!drag)return;ry+=(e.clientX-lx)*.008;rx+=(e.clientY-ly)*.008;lx=e.clientX;ly=e.clientY;draw();}});canvas.addEventListener('wheel',e=>{{e.preventDefault();zoom*=Math.exp(-e.deltaY*.001);zoom=Math.max(.2,Math.min(5,zoom));draw();}},{{passive:false}});addEventListener('resize',resize);resize();
</script></body></html>"""
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(html, encoding="utf-8")
    print(f"Saved HTML layout = {output}")


if __name__ == "__main__":
    main()
