#!/usr/bin/env python3
import argparse
import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection

from config_io import (
    apply_telescope_frame_to_facets,
    expand_component_config,
    load_facets_csv,
    load_facets_from_config,
    point_to_global,
    parse_vec3,
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
    v = v / np.linalg.norm(v)
    return u, v, n


def aperture_polygon(facet):
    center = facet["center"]
    if "u_axis" in facet and "v_axis" in facet:
        u = np.array(facet["u_axis"], dtype=float)
        v = np.array(facet["v_axis"], dtype=float)
        u = u / np.linalg.norm(u)
        v = v / np.linalg.norm(v)
    else:
        u, v, _ = local_frame(facet["normal"])
    shape = facet["shape"]
    size = facet["size1"]

    if shape == "hexagon":
        a = 0.5 * size
        local = np.array([
            [a, a / math.sqrt(3.0)],
            [0.0, 2.0 * a / math.sqrt(3.0)],
            [-a, a / math.sqrt(3.0)],
            [-a, -a / math.sqrt(3.0)],
            [0.0, -2.0 * a / math.sqrt(3.0)],
            [a, -a / math.sqrt(3.0)],
        ])
    elif shape == "square":
        h = 0.5 * size
        local = np.array([[h, h], [-h, h], [-h, -h], [h, -h]])
    else:
        r = size
        angles = np.linspace(0, 2 * math.pi, 48, endpoint=False)
        local = np.column_stack([r * np.cos(angles), r * np.sin(angles)])

    rot = facet.get("rotation", 0.0)
    if rot != 0.0:
        c = math.cos(rot)
        s = math.sin(rot)
        R = np.array([[c, -s], [s, c]])
        local = local @ R.T

    return np.array([center + x * u + y * v for x, y in local])


def reflect_direction(in_dir, normal):
    in_dir = in_dir / np.linalg.norm(in_dir)
    normal = normal / np.linalg.norm(normal)
    out = in_dir - 2.0 * np.dot(in_dir, normal) * normal
    return out / np.linalg.norm(out)


def intersect_plane(point, direction, plane_point, plane_normal):
    denom = np.dot(direction, plane_normal)
    if abs(denom) < 1e-14:
        return None
    t = np.dot(plane_point - point, plane_normal) / denom
    if t <= 0:
        return None
    return point + t * direction


def make_plane_patch(plane_point, plane_normal, radius):
    u, v, _ = local_frame(plane_normal)
    return make_oriented_plane_patch(plane_point, u, v, radius)


def make_oriented_plane_patch(plane_point, u_axis, v_axis, radius):
    u = u_axis / np.linalg.norm(u_axis)
    v = v_axis / np.linalg.norm(v_axis)
    s = radius
    corners = np.array([
        plane_point + (-s) * u + (-s) * v,
        plane_point + s * u + (-s) * v,
        plane_point + s * u + s * v,
        plane_point + (-s) * u + s * v,
    ])
    return corners


def load_obstruction_cells(mask_csv):
    mask_csv = Path(mask_csv)
    meta = {}
    cells = []
    with mask_csv.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                item = line[1:].strip()
                if "=" in item:
                    key, value = item.split("=", 1)
                    meta[key.strip()] = value.strip()
                continue
            if line.lower().startswith("ix,"):
                continue
            parts = [p.strip() for p in line.split(",")]
            if len(parts) < 2:
                continue
            cells.append((int(parts[0]), int(parts[1])))

    required = ("x_min_m", "y_min_m", "cell_size_m", "plane_z_m")
    missing = [key for key in required if key not in meta]
    if missing:
        raise ValueError(f"obstruction mask is missing metadata: {', '.join(missing)}")

    x_min = float(meta["x_min_m"])
    y_min = float(meta["y_min_m"])
    cell = float(meta["cell_size_m"])
    z = float(meta["plane_z_m"])
    points = np.array(
        [[x_min + (ix + 0.5) * cell, y_min + (iy + 0.5) * cell, z] for ix, iy in cells],
        dtype=float,
    )
    return points, meta


def load_obstruction_primitives(primitives_csv):
    primitives = []
    with Path(primitives_csv).open(newline="") as handle:
        rows = csv.DictReader(line for line in handle if not line.lstrip().startswith("#"))
        for row in rows:
            def value(key, default=0.0):
                text = row.get(key, "")
                if text in ("", None):
                    return default
                return float(text)
            primitives.append(
                {
                    "type": row.get("type", "").lower(),
                    "name": row.get("name", ""),
                    "role": row.get("role", ""),
                    "p0": np.array([value("x0_m"), value("y0_m"), value("z0_m")]),
                    "p1": np.array([value("x1_m"), value("y1_m"), value("z1_m")]),
                    "center": np.array([
                        value("center_x_m", value("x0_m")),
                        value("center_y_m", value("y0_m")),
                        value("center_z_m", value("z0_m")),
                    ]),
                    "radius": value("radius_m"),
                    "height": value("height_m"),
                    "rotation": value("rotation_rad"),
                    "sides": int(value("sides")),
                    "half": np.array([value("half_x_m"), value("half_y_m"), value("half_z_m")]),
                    "hole_radius": value("hole_radius_m"),
                    "hole_rotation": value("hole_rotation_rad"),
                    "hole_sides": int(value("hole_sides")),
                }
            )
    return primitives


def box_corners(center, half):
    x, y, z = center
    hx, hy, hz = half
    return np.array([
        [x - hx, y - hy, z - hz],
        [x + hx, y - hy, z - hz],
        [x + hx, y + hy, z - hz],
        [x - hx, y + hy, z - hz],
        [x - hx, y - hy, z + hz],
        [x + hx, y - hy, z + hz],
        [x + hx, y + hy, z + hz],
        [x - hx, y + hy, z + hz],
    ])


def box_faces(center, half):
    c = box_corners(center, half)
    return [
        c[[0, 1, 2, 3]],
        c[[4, 5, 6, 7]],
        c[[0, 1, 5, 4]],
        c[[2, 3, 7, 6]],
        c[[1, 2, 6, 5]],
        c[[0, 3, 7, 4]],
    ]


def box_edges(center, half):
    c = box_corners(center, half)
    return [
        c[[0, 1, 2, 3, 0]],
        c[[4, 5, 6, 7, 4]],
        c[[0, 4]],
        c[[1, 5]],
        c[[2, 6]],
        c[[3, 7]],
    ]


def regular_prism_faces(center, radius, height, rotation, sides):
    if radius <= 0 or height <= 0 or sides < 3:
        return []
    angles = rotation + np.arange(sides) * (2 * math.pi / sides)
    xy = np.column_stack([
        center[0] + radius * np.cos(angles),
        center[1] + radius * np.sin(angles),
    ])
    bottom = np.column_stack([xy, np.full(sides, center[2] - 0.5 * height)])
    top = np.column_stack([xy, np.full(sides, center[2] + 0.5 * height)])
    faces = [bottom, top]
    for i in range(sides):
        j = (i + 1) % sides
        faces.append(np.array([bottom[i], bottom[j], top[j], top[i]]))
    return faces


def regular_prism_edges(center, radius, height, rotation, sides):
    faces = regular_prism_faces(center, radius, height, rotation, sides)
    if len(faces) < 2:
        return []
    bottom = faces[0]
    top = faces[1]
    edges = [np.vstack([bottom, bottom[0]]), np.vstack([top, top[0]])]
    for a, b in zip(bottom, top):
        edges.append(np.array([a, b]))
    return edges


def transform_points(points, frame):
    if frame is None:
        return np.array(points, dtype=float)
    return np.array([point_to_global(np.array(p, dtype=float), frame) for p in points])


def set_equal_3d(ax, points, pad=0.08):
    pts = np.asarray(points)
    mins = pts.min(axis=0)
    maxs = pts.max(axis=0)
    center = 0.5 * (mins + maxs)
    span = max(maxs - mins)
    if span <= 0:
        span = 1.0
    span *= (1.0 + pad)
    half = 0.5 * span
    ax.set_xlim(center[0] - half, center[0] + half)
    ax.set_ylim(center[1] - half, center[1] + half)
    ax.set_zlim(center[2] - half, center[2] + half)
    try:
        ax.set_box_aspect((1, 1, 1))
    except Exception:
        pass


def main():
    parser = argparse.ArgumentParser(
        description="Plot a publication-style 3D view of mirror facets and the output plane."
    )
    parser.add_argument("--config", help="run_optical_sim key=value config")
    parser.add_argument("--mirror-csv", help="MirrorFacet CSV; overrides config mirror.csv_path")
    parser.add_argument("--elevation-deg", type=float, default=None, help="override telescope/mirror elevation state")
    parser.add_argument("--output", default="optical_layout_3d.png", help="output image path")
    parser.add_argument("--format", default=None, help="matplotlib format, e.g. png or pdf")
    parser.add_argument("--dpi", type=int, default=300, help="output DPI")
    parser.add_argument("--ray-stride", type=int, default=1, help="draw every Nth center ray")
    parser.add_argument("--normal-scale", type=float, default=0.35, help="normal arrow length in m")
    parser.add_argument("--plane-radius", type=float, default=None, help="output plane half-size in m")
    parser.add_argument(
        "--show-obstruction",
        action="store_true",
        help="overlay obstruction mask cells from obstruction.config",
    )
    parser.add_argument(
        "--obstruction-mask",
        help="explicit obstruction mask CSV; overrides config obstruction.mask_csv",
    )
    parser.add_argument(
        "--obstruction-primitives",
        help="explicit obstruction primitives CSV; overrides config obstruction.primitives_csv",
    )
    parser.add_argument(
        "--obstruction-stride",
        type=int,
        default=1,
        help="draw every Nth obstruction mask cell",
    )
    parser.add_argument(
        "--show-camera-axes",
        action="store_true",
        help="draw output-plane +u/+v image axes on the static layout",
    )
    parser.add_argument("--view", default="32,-58", help="elev,azim camera view")
    args = parser.parse_args()

    cfg_path = Path(args.config).resolve() if args.config else None
    cfg = {}
    if cfg_path:
        cfg, _ = expand_component_config(cfg_path)
        if args.elevation_deg is not None:
            cfg["telescope.pointing_el_deg"] = str(args.elevation_deg)
            cfg["mirror.series_elevation_deg"] = str(args.elevation_deg)
    frame = telescope_frame_from_config(cfg) if cfg else None

    plane_point = parse_vec3(cfg.get("output.plane_point"), [0.0, 0.0, 0.0])
    plane_normal = parse_vec3(cfg.get("output.plane_normal"), [0.0, 0.0, 1.0])
    plane_u_axis = parse_vec3(cfg["output.plane_u_axis"]) if "output.plane_u_axis" in cfg else None
    plane_v_axis = parse_vec3(cfg["output.plane_v_axis"]) if "output.plane_v_axis" in cfg else None
    beam_dir = source_direction_from_config(cfg)
    if frame is not None:
        plane_point = point_to_global(plane_point, frame)
        plane_normal = rotate_local_vector(plane_normal, frame)
        if plane_u_axis is not None:
            plane_u_axis = rotate_local_vector(plane_u_axis, frame)
        if plane_v_axis is not None:
            plane_v_axis = rotate_local_vector(plane_v_axis, frame)
        beam_dir = rotate_local_vector(beam_dir, frame)
    plane_normal = plane_normal / np.linalg.norm(plane_normal)
    if plane_u_axis is None or plane_v_axis is None:
        plane_u_axis, plane_v_axis, _ = local_frame(plane_normal)
    else:
        plane_u_axis = plane_u_axis / np.linalg.norm(plane_u_axis)
        plane_v_axis = plane_v_axis / np.linalg.norm(plane_v_axis)
    beam_dir = beam_dir / np.linalg.norm(beam_dir)

    if args.mirror_csv:
        facets = load_facets_csv(Path(args.mirror_csv).resolve())
    elif cfg_path:
        facets = load_facets_from_config(cfg_path, cfg)
    else:
        raise SystemExit("mirror geometry is required via --mirror-csv or --config")
    if frame is not None:
        facets = apply_telescope_frame_to_facets(facets, frame)
    centers = np.array([f["center"] for f in facets])
    radial = np.linalg.norm(centers[:, :2], axis=1)

    obstruction_points = None
    obstruction_primitives = []
    if args.show_obstruction or args.obstruction_mask or args.obstruction_primitives:
        mode = cfg.get("obstruction.mode", "mask").lower()
        obstruction_enabled = cfg.get("obstruction.enabled", "false").lower() in {
            "1",
            "true",
            "yes",
            "on",
        }
        if not obstruction_enabled and not args.obstruction_mask and not args.obstruction_primitives:
            raise SystemExit("config obstruction.enabled is false; use --obstruction-mask to force drawing")
        if args.obstruction_primitives or mode == "primitives":
            primitives_path = args.obstruction_primitives or cfg.get("obstruction.primitives_csv")
            if not primitives_path:
                raise SystemExit("primitive obstruction drawing requires obstruction.primitives_csv")
            obstruction_primitives = load_obstruction_primitives(primitives_path)
        else:
            mask_path = args.obstruction_mask or cfg.get("obstruction.mask_csv")
            if not mask_path:
                raise SystemExit("--show-obstruction requires obstruction.mask_csv in config")
            obstruction_points, _ = load_obstruction_cells(mask_path)
            stride = max(1, args.obstruction_stride)
            obstruction_points = obstruction_points[::stride]
            if frame is not None:
                obstruction_points = np.array([point_to_global(p, frame) for p in obstruction_points])

    polygons = [aperture_polygon(f) for f in facets]
    rays = []
    for i, facet in enumerate(facets):
        if args.ray_stride <= 0 or i % args.ray_stride != 0:
            continue
        out_dir = reflect_direction(beam_dir, facet["normal"])
        hit = intersect_plane(facet["center"], out_dir, plane_point, plane_normal)
        if hit is not None:
            rays.append(np.array([facet["center"], hit]))

    if args.plane_radius is None:
        xy_span = max(np.ptp(centers[:, 0]), np.ptp(centers[:, 1]), 1.0)
        plane_radius = max(0.7, 0.18 * xy_span)
    else:
        plane_radius = args.plane_radius
    plane_patch = make_oriented_plane_patch(plane_point, plane_u_axis, plane_v_axis, plane_radius)

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.labelsize": 10,
        "axes.titlesize": 11,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "figure.dpi": args.dpi,
        "savefig.dpi": args.dpi,
    })

    fig = plt.figure(figsize=(7.2, 6.4))
    ax = fig.add_subplot(111, projection="3d")
    ax.set_facecolor("white")

    cmap = plt.get_cmap("viridis")
    norm = plt.Normalize(radial.min(), radial.max())
    facecolors = [cmap(norm(r)) for r in radial]
    poly = Poly3DCollection(
        polygons,
        facecolors=facecolors,
        edgecolors=(0.15, 0.15, 0.15, 0.68),
        linewidths=0.35,
        alpha=0.93,
    )
    ax.add_collection3d(poly)

    plane = Poly3DCollection(
        [plane_patch],
        facecolors=(0.78, 0.84, 0.92, 0.25),
        edgecolors=(0.12, 0.22, 0.35, 0.85),
        linewidths=1.2,
    )
    ax.add_collection3d(plane)

    if rays:
        ray_collection = Line3DCollection(rays, colors=(0.82, 0.22, 0.16, 0.32), linewidths=0.55)
        ax.add_collection3d(ray_collection)

    normal_starts = centers
    normal_dirs = np.array([f["normal"] for f in facets])
    ax.quiver(
        normal_starts[:, 0],
        normal_starts[:, 1],
        normal_starts[:, 2],
        normal_dirs[:, 0],
        normal_dirs[:, 1],
        normal_dirs[:, 2],
        length=args.normal_scale,
        normalize=True,
        color=(0.05, 0.05, 0.05, 0.45),
        linewidth=0.45,
        arrow_length_ratio=0.25,
    )

    ax.scatter(
        [plane_point[0]],
        [plane_point[1]],
        [plane_point[2]],
        s=32,
        marker="o",
        color="#B2182B",
        depthshade=False,
        label="Output plane center",
    )

    if args.show_camera_axes:
        axis_length = max(0.75, plane_radius * 1.45)
        for direction, color, label in (
            (plane_u_axis, "#F2B701", "+u / image +x"),
            (plane_v_axis, "#2CA25F", "+v / image +y"),
            (plane_normal, "#B2182B", "camera normal / toward mirror"),
        ):
            end = plane_point + direction * axis_length
            ax.quiver(
                [plane_point[0]],
                [plane_point[1]],
                [plane_point[2]],
                [direction[0]],
                [direction[1]],
                [direction[2]],
                length=axis_length,
                normalize=True,
                color=color,
                linewidth=2.0,
                arrow_length_ratio=0.18,
                label=label,
            )
            ax.text(
                end[0],
                end[1],
                end[2],
                label,
                color=color,
                fontsize=8,
                weight="bold",
            )
        global_z = np.array([0.0, 0.0, 1.0])
        global_z_start = plane_point - 0.18 * plane_u_axis - 0.18 * plane_v_axis
        global_z_end = global_z_start + global_z * axis_length
        ax.quiver(
            [global_z_start[0]],
            [global_z_start[1]],
            [global_z_start[2]],
            [global_z[0]],
            [global_z[1]],
            [global_z[2]],
            length=axis_length,
            normalize=True,
            color="#111111",
            linewidth=1.8,
            arrow_length_ratio=0.18,
            label="+global z / up",
        )
        ax.text(
            global_z_end[0],
            global_z_end[1],
            global_z_end[2],
            "+global z / up",
            color="#111111",
            fontsize=8,
            weight="bold",
        )

    if obstruction_points is not None and len(obstruction_points) > 0:
        ax.scatter(
            obstruction_points[:, 0],
            obstruction_points[:, 1],
            obstruction_points[:, 2],
            s=1.2,
            marker="s",
            color=(0.02, 0.02, 0.02, 0.34),
            depthshade=False,
            label="Support obstruction mask",
        )
    if obstruction_primitives:
        cylinder_segments = []
        cylinder_widths = []
        solid_polys = []
        wire_segments = []
        wire_colors = []
        wire_widths = []
        for primitive in obstruction_primitives:
            if primitive["type"] == "cylinder":
                cylinder_segments.append(transform_points([primitive["p0"], primitive["p1"]], frame))
                cylinder_widths.append(max(0.8, primitive["radius"] * 22.0))
            elif primitive["type"] in {"box", "aabb"}:
                faces = [transform_points(face, frame) for face in box_faces(primitive["center"], primitive["half"])]
                solid_polys.extend(faces)
                for edge in box_edges(primitive["center"], primitive["half"]):
                    wire_segments.append(transform_points(edge, frame))
                    wire_colors.append((0.10, 0.44, 0.22, 0.95))
                    wire_widths.append(1.1)
                if primitive["hole_radius"] > 0 and primitive["hole_sides"] >= 3:
                    for edge in regular_prism_edges(
                            primitive["center"],
                            primitive["hole_radius"],
                            primitive["half"][2] * 2.0,
                            primitive["hole_rotation"],
                            primitive["hole_sides"],
                    ):
                        wire_segments.append(transform_points(edge, frame))
                        wire_colors.append((0.78, 0.05, 0.14, 1.0))
                        wire_widths.append(1.4)
            elif primitive["type"] == "polygon_prism":
                faces = [
                    transform_points(face, frame)
                    for face in regular_prism_faces(
                        primitive["center"],
                        primitive["radius"],
                        primitive["height"],
                        primitive["rotation"],
                        primitive["sides"],
                    )
                ]
                solid_polys.extend(faces)
                for edge in regular_prism_edges(
                    primitive["center"],
                    primitive["radius"],
                    primitive["height"],
                    primitive["rotation"],
                    primitive["sides"],
                ):
                    wire_segments.append(transform_points(edge, frame))
                    wire_colors.append((0.0, 0.58, 0.68, 1.0))
                    wire_widths.append(1.2)
        if cylinder_segments:
            ax.add_collection3d(
                Line3DCollection(
                    cylinder_segments,
                    colors=(0.03, 0.03, 0.03, 0.72),
                    linewidths=cylinder_widths,
                    label="3D support obstruction",
                ),
                autolim=False,
            )
        if solid_polys:
            ax.add_collection3d(
                Poly3DCollection(
                    solid_polys,
                    facecolors=(0.08, 0.08, 0.08, 0.13),
                    edgecolors=(0.03, 0.03, 0.03, 0.30),
                    linewidths=0.35,
                ),
                autolim=False,
            )
        if wire_segments:
            ax.add_collection3d(
                Line3DCollection(
                    wire_segments,
                    colors=wire_colors,
                    linewidths=wire_widths,
                ),
                autolim=False,
            )

    all_points = []
    all_points.extend(centers)
    all_points.extend(np.concatenate(polygons))
    all_points.extend(plane_patch)
    if rays:
        all_points.extend(np.concatenate(rays))
    if obstruction_points is not None and len(obstruction_points) > 0:
        all_points.extend(obstruction_points)
    if obstruction_primitives:
        for primitive in obstruction_primitives:
            all_points.extend(transform_points([primitive["p0"], primitive["p1"]], frame))
            if primitive["type"] in {"box", "aabb"}:
                all_points.extend(transform_points(box_corners(primitive["center"], primitive["half"]), frame))
            if primitive["type"] == "polygon_prism":
                for face in regular_prism_faces(
                    primitive["center"],
                    primitive["radius"],
                    primitive["height"],
                    primitive["rotation"],
                    primitive["sides"],
                ):
                    all_points.extend(transform_points(face, frame))
    set_equal_3d(ax, np.array(all_points))

    elev, azim = [float(x) for x in args.view.split(",", 1)]
    ax.view_init(elev=elev, azim=azim)
    ax.set_xlabel("x [m]", labelpad=8)
    ax.set_ylabel("y [m]", labelpad=8)
    ax.set_zlabel("z [m]", labelpad=8)
    ax.set_title("Mirror facets and optical output plane")

    mappable = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    mappable.set_array(radial)
    cbar = fig.colorbar(mappable, ax=ax, shrink=0.62, pad=0.08)
    cbar.set_label("facet radial position [m]")
    if (obstruction_points is not None and len(obstruction_points) > 0) or obstruction_primitives:
        ax.legend(loc="upper right", frameon=True, framealpha=0.9, fontsize=8)

    ax.text2D(
        0.02,
        0.02,
        (
            f"N facets = {len(facets)}\n"
            f"El = {float(cfg.get('mirror.series_elevation_deg', cfg.get('telescope.pointing_el_deg', 90.0))):.1f} deg\n"
            f"Output plane = ({plane_point[0]:.2f}, {plane_point[1]:.2f}, {plane_point[2]:.2f}) m"
        ),
        transform=ax.transAxes,
        fontsize=8,
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="0.82", alpha=0.88),
    )

    ax.grid(True, alpha=0.18)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.set_facecolor((1.0, 1.0, 1.0, 0.0))
        axis.pane.set_edgecolor((0.82, 0.82, 0.82, 1.0))

    fig.tight_layout()
    output = Path(args.output)
    fig.savefig(output, format=args.format, bbox_inches="tight")
    print(f"Saved 3D layout = {output}")

    if "agg" not in plt.get_backend().lower():
        plt.show()


if __name__ == "__main__":
    main()
