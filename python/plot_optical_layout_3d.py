#!/usr/bin/env python3
import argparse
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
    with Path(primitives_csv).open() as handle:
        reader = None
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if reader is None:
                reader = [item.strip() for item in line.split(",")]
                continue
            cells = [item.strip() for item in line.split(",")]
            if len(cells) < 12:
                continue
            primitives.append(
                {
                    "type": cells[0].lower(),
                    "name": cells[1],
                    "p0": np.array([float(cells[2]), float(cells[3]), float(cells[4])]),
                    "p1": np.array([float(cells[5]), float(cells[6]), float(cells[7])]),
                    "radius": float(cells[8]),
                    "half": np.array([float(cells[9]), float(cells[10]), float(cells[11])]),
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
    beam_dir = source_direction_from_config(cfg)
    if frame is not None:
        plane_point = point_to_global(plane_point, frame)
        plane_normal = rotate_local_vector(plane_normal, frame)
        beam_dir = rotate_local_vector(beam_dir, frame)
    plane_normal = plane_normal / np.linalg.norm(plane_normal)
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
            if frame is not None:
                for primitive in obstruction_primitives:
                    primitive["p0"] = point_to_global(primitive["p0"], frame)
                    primitive["p1"] = point_to_global(primitive["p1"], frame)
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
    plane_patch = make_plane_patch(plane_point, plane_normal, plane_radius)

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
        segments = []
        widths = []
        box_polys = []
        for primitive in obstruction_primitives:
            if primitive["type"] == "cylinder":
                segments.append(np.array([primitive["p0"], primitive["p1"]]))
                widths.append(max(0.8, primitive["radius"] * 22.0))
            elif primitive["type"] in {"box", "aabb"}:
                box_polys.extend(box_faces(primitive["p0"], primitive["half"]))
        if segments:
            ax.add_collection3d(
                Line3DCollection(
                    segments,
                    colors=(0.03, 0.03, 0.03, 0.72),
                    linewidths=widths,
                    label="3D support obstruction",
                )
            )
        if box_polys:
            ax.add_collection3d(
                Poly3DCollection(
                    box_polys,
                    facecolors=(0.04, 0.04, 0.04, 0.25),
                    edgecolors=(0.03, 0.03, 0.03, 0.68),
                    linewidths=0.6,
                )
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
            all_points.append(primitive["p0"])
            all_points.append(primitive["p1"])
            if primitive["type"] in {"box", "aabb"}:
                all_points.extend(box_corners(primitive["p0"], primitive["half"]))
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
