#!/usr/bin/env python3
"""Plot one real parallel-light C++ run without changing its coordinates."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection, PatchCollection
from matplotlib.patches import Rectangle

from build_coordinate_diagnostics_html import read_geometry


INK = "#17212b"
MUTED = "#687684"
GRID = "#d9e0e6"
STRUCTURE = "#7f95a8"
REFLECTION = "#d99a16"
INCOMING = "#2378a7"
BLOCKED_IN = "#c63f52"
BLOCKED_OUT = "#7755a6"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(line for line in handle if not line.lstrip().startswith("#")))


def choose_font() -> None:
    candidates = ["Microsoft YaHei", "Microsoft YaHei UI", "SimHei", "Noto Sans CJK SC"]
    installed = {font.name for font in mpl.font_manager.fontManager.ttflist}
    for candidate in candidates:
        if candidate in installed:
            mpl.rcParams["font.sans-serif"] = [candidate, "DejaVu Sans"]
            break
    mpl.rcParams["axes.unicode_minus"] = False
    mpl.rcParams["figure.dpi"] = 150
    mpl.rcParams["savefig.dpi"] = 180
    mpl.rcParams["font.size"] = 10
    mpl.rcParams["axes.titleweight"] = "normal"


def expand_local(point: list[float], basis: dict[str, list[float]]) -> np.ndarray:
    point_array = np.asarray(point, dtype=float)
    matrix = np.asarray([basis["x"], basis["y"], basis["z"]], dtype=float).T
    return matrix @ point_array


def global_to_local(points: np.ndarray, basis: dict[str, list[float]]) -> np.ndarray:
    matrix = np.asarray([basis["x"], basis["y"], basis["z"]], dtype=float)
    return points @ matrix.T


def set_equal_3d(ax, points: np.ndarray, padding: float = 0.08) -> None:
    low = points.min(axis=0)
    high = points.max(axis=0)
    center = (low + high) / 2.0
    half = max(high - low) * (0.5 + padding)
    ax.set_xlim(center[0] - half, center[0] + half)
    ax.set_ylim(center[1] - half, center[1] + half)
    ax.set_zlim(center[2] - half, center[2] + half)
    ax.set_box_aspect((1, 1, 1))


def configure_axis(ax, xlabel: str, ylabel: str) -> None:
    ax.set_xlabel(xlabel, color=INK)
    ax.set_ylabel(ylabel, color=INK)
    ax.tick_params(colors=MUTED, labelsize=8)
    ax.grid(True, color=GRID, linewidth=0.55)
    ax.set_facecolor("#ffffff")


def plot_optics(case: dict, ideal: list[dict], primitives: list[dict], output: Path) -> None:
    basis = case["basis"]
    reflection_global = np.asarray(case["full_reflection_points_global_m"], dtype=float)
    blocked_in_global = np.asarray(
        case["full_incoming_blocked_mirror_endpoints_global_m"], dtype=float
    )
    blocked_out_global = np.asarray(
        case["full_reflected_blocked_surface_endpoints_global_m"], dtype=float
    )
    reflection_local = global_to_local(reflection_global, basis)
    blocked_in_local = global_to_local(blocked_in_global, basis)
    blocked_out_local = global_to_local(blocked_out_global, basis)
    output_uv = np.asarray(case["full_output_uv_m"], dtype=float)

    fig = plt.figure(figsize=(15.2, 8.6), facecolor="white")
    grid = fig.add_gridspec(2, 2, width_ratios=(1.42, 1.0), hspace=0.28, wspace=0.2)
    ax3d = fig.add_subplot(grid[:, 0], projection="3d")
    ax_mirror = fig.add_subplot(grid[0, 1])
    ax_output = fig.add_subplot(grid[1, 1])

    # Full ideal mirror and support/obstruction geometry, transformed by the logged basis.
    structure_points: list[np.ndarray] = []
    for facet in ideal:
        polygon = np.asarray([expand_local(point, basis) for point in facet["polygon"]])
        polygon = np.vstack([polygon, polygon[0]])
        ax3d.plot(polygon[:, 0], polygon[:, 1], polygon[:, 2], color="#9aa8b4", lw=0.55, alpha=0.7)
        structure_points.extend(polygon)
        local_polygon = np.asarray(facet["polygon"])
        local_polygon = np.vstack([local_polygon, local_polygon[0]])
        ax_mirror.plot(local_polygon[:, 0], local_polygon[:, 1], color="#b8c3cc", lw=0.5, zorder=1)
    for primitive in primitives:
        segment = np.asarray([expand_local(point, basis) for point in primitive["points"]])
        ax3d.plot(segment[:, 0], segment[:, 1], segment[:, 2], color=STRUCTURE, lw=0.55, alpha=0.45)
        structure_points.extend(segment)

    # Complete point populations.
    ax3d.scatter(
        reflection_global[:, 0], reflection_global[:, 1], reflection_global[:, 2],
        s=1.15, color=REFLECTION, alpha=0.34, depthshade=False,
        label=f"全部镜面反射点 {len(reflection_global):,}",
    )
    ax3d.scatter(
        blocked_in_global[:, 0], blocked_in_global[:, 1], blocked_in_global[:, 2],
        s=5.0, marker="x", linewidths=0.35, color=BLOCKED_IN, alpha=0.48,
        depthshade=False, label=f"全部入射遮挡端点 {len(blocked_in_global):,}",
    )
    ax3d.scatter(
        blocked_out_global[:, 0], blocked_out_global[:, 1], blocked_out_global[:, 2],
        s=8.0, marker="x", linewidths=0.45, color=BLOCKED_OUT, alpha=0.75,
        depthshade=False, label=f"全部反射遮挡端点 {len(blocked_out_global):,}",
    )

    # Only ray segments are sampled, exactly as in the HTML page.
    for ray in case["rays"]:
        incoming = np.asarray([ray["input"], ray["mirror"]], dtype=float)
        if ray["blocked_incoming"]:
            ax3d.plot(*incoming.T, color=BLOCKED_IN, lw=0.8, alpha=0.52, ls="--")
            continue
        ax3d.plot(*incoming.T, color=INCOMING, lw=0.55, alpha=0.34)
        reflected = np.asarray([ray["mirror"], ray["surface"]], dtype=float)
        ax3d.plot(
            *reflected.T,
            color=BLOCKED_OUT if ray["blocked_reflected"] else REFLECTION,
            lw=0.8,
            alpha=0.58,
            ls="--" if ray["blocked_reflected"] else "-",
        )

    # Horizontal ground reference and program directions.
    dish = expand_local([0.0, 0.0, -16.0], basis)
    ground_z = dish[2] - 1.0
    gx, gy = np.meshgrid(np.linspace(dish[0] - 9, dish[0] + 9, 2), np.linspace(-9, 9, 2))
    ax3d.plot_surface(gx, gy, np.full_like(gx, ground_z), color="#7aa58b", alpha=0.08, shade=False)
    pointing = np.asarray(basis["z"], dtype=float)
    ax3d.quiver(dish[0], dish[1], ground_z, *pointing, length=10, color=INK, lw=1.4, arrow_length_ratio=0.1)
    ax3d.text(*(np.asarray([dish[0], dish[1], ground_z]) + pointing * 10.5), "望远镜光轴 el=70°", color=INK, fontsize=9)
    ray_direction = np.asarray(case["rays"][0]["input_dir"], dtype=float)
    source_anchor = np.asarray(case["rays"][0]["input"], dtype=float) - ray_direction * 3.0
    ax3d.quiver(*source_anchor, *ray_direction, length=5, color=INCOMING, lw=1.4, arrow_length_ratio=0.12)
    ax3d.text(*source_anchor, "71°光源传播方向", color=INCOMING, fontsize=9)

    all_extent = np.vstack([
        np.asarray(structure_points), reflection_global, blocked_in_global,
        blocked_out_global, np.asarray([ray["input"] for ray in case["rays"]]),
    ])
    set_equal_3d(ax3d, all_extent)
    ax3d.view_init(elev=19, azim=-62)
    ax3d.set_xlabel("通用全局 x / m", color=INK)
    ax3d.set_ylabel("通用全局 y / m", color=INK)
    ax3d.set_zlabel("通用全局 z / m", color=INK)
    ax3d.tick_params(colors=MUTED, labelsize=8)
    ax3d.set_title("三维望远镜与光路（点全量、线抽样）", loc="left", color=INK, pad=12)
    ax3d.legend(loc="upper left", fontsize=8, frameon=False)
    ax3d.xaxis.pane.set_alpha(0.0)
    ax3d.yaxis.pane.set_alpha(0.0)
    ax3d.zaxis.pane.set_alpha(0.0)

    ax_mirror.scatter(reflection_local[:, 0], reflection_local[:, 1], s=1.0, color=REFLECTION, alpha=0.38, rasterized=True, label="真实反射点")
    ax_mirror.scatter(blocked_in_local[:, 0], blocked_in_local[:, 1], s=5.0, marker="x", linewidths=0.35, color=BLOCKED_IN, alpha=0.55, rasterized=True, label="入射遮挡理论镜面端点")
    ax_mirror.set_aspect("equal", adjustable="box")
    configure_axis(ax_mirror, "镜面 local x / m", "镜面 local y / m")
    ax_mirror.set_title("镜面正视图：完整落点", loc="left", color=INK)
    ax_mirror.legend(loc="upper right", fontsize=8, frameon=False)

    ax_output.scatter(output_uv[:, 0], output_uv[:, 1], s=1.0, color="#73808b", alpha=0.25, rasterized=True, label="无遮挡 output u/v")
    ax_output.scatter(blocked_out_local[:, 0], blocked_out_local[:, 1], s=10.0, marker="x", linewidths=0.55, color=BLOCKED_OUT, alpha=0.8, rasterized=True, label="反射遮挡理论输出端点")
    ax_output.axhline(0.0, color=GRID, lw=0.7)
    ax_output.axvline(0.0, color=GRID, lw=0.7)
    output_and_blocked = np.vstack([output_uv, blocked_out_local[:, :2]])
    output_center = (output_and_blocked.min(axis=0) + output_and_blocked.max(axis=0)) / 2.0
    output_half = max(output_and_blocked.max(axis=0) - output_and_blocked.min(axis=0)) * 0.58
    ax_output.set_xlim(output_center[0] - output_half, output_center[0] + output_half)
    ax_output.set_ylim(output_center[1] - output_half, output_center[1] + output_half)
    ax_output.set_aspect("equal", adjustable="box")
    configure_axis(ax_output, "原始 u / m", "原始 v / m")
    ax_output.set_title("输出面：完整无遮挡输出与反射遮挡端点", loc="left", color=INK)
    ax_output.legend(loc="upper right", fontsize=8, frameon=False)

    summary = case["summary"]
    fig.suptitle("71°平行光（望远镜 el=70°）真实 C++ 光学输出", x=0.055, ha="left", y=0.985, color=INK, fontsize=15)
    fig.text(
        0.055, 0.947,
        f"输入 30,000｜镜面反射 {summary['hit_mirror']:,}｜无遮挡输出 {summary['hit_output_plane']:,}｜入射遮挡 {summary['blocked_incoming']:,}｜反射遮挡 {summary['blocked_reflected']:,}｜仅48条光路线抽样",
        ha="left", color=MUTED, fontsize=9,
    )
    fig.text(0.055, 0.02, "数据：coordinate_parallel_sky_up 的完整 hits CSV；红/紫叉是程序保存的理论端点，不是遮挡 primitive 精确交点。", color=MUTED, fontsize=8)
    fig.subplots_adjust(left=0.055, right=0.98, top=0.91, bottom=0.08)
    fig.savefig(output, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def pixel_patches(camera_geometry: list[dict[str, str]], signal: dict[int, dict], max_count: int):
    patches = []
    colors = []
    for row in camera_geometry:
        u = float(row["x_m"])
        v = float(row["y_m"])
        size = float(row["size_m"])
        patches.append(Rectangle((u - size / 2, v - size / 2), size, size))
        count = int(signal.get(int(row["id"]), {}).get("photon_count", 0))
        colors.append(math.sqrt(count / max_count) if count > 0 else 0.0)
    return patches, np.asarray(colors)


def draw_camera_panel(ax, camera_geometry, signal, output_uv, camera_uv, limits, title):
    max_count = max(int(row["photon_count"]) for row in signal.values())
    patches, values = pixel_patches(camera_geometry, signal, max_count)
    collection = PatchCollection(
        patches, cmap=mpl.colors.LinearSegmentedColormap.from_list("camera", ["#edf1f4", "#f3c86a", "#b66d08"]),
        edgecolor="#b6c1ca", linewidth=0.22,
    )
    collection.set_array(values)
    collection.set_clim(0.0, 1.0)
    ax.add_collection(collection)
    ax.scatter(output_uv[:, 0], output_uv[:, 1], s=1.0, color="#54616b", alpha=0.22, rasterized=True, label="全部 output u/v")
    ax.scatter(camera_uv[:, 0], camera_uv[:, 1], s=1.2, color="#16222c", alpha=0.32, rasterized=True, label="全部 camera-hit u/v")
    ax.axhline(0.0, color=GRID, lw=0.7)
    ax.axvline(0.0, color=GRID, lw=0.7)
    ax.set_xlim(*limits[0])
    ax.set_ylim(*limits[1])
    ax.set_aspect("equal", adjustable="box")
    configure_axis(ax, "原始 u / m", "原始 v / m")
    ax.set_title(title, loc="left", color=INK)
    ax.legend(loc="upper right", fontsize=8, frameon=False, markerscale=4)
    return collection


def plot_camera(case: dict, camera_geometry: list[dict[str, str]], output: Path) -> None:
    output_uv = np.asarray(case["full_output_uv_m"], dtype=float)
    camera_uv = np.asarray(case["full_camera_hit_uv_m"], dtype=float)
    signal = {int(row["pixel_id"]): row for row in case["camera_signal"]}
    geometry_u = np.asarray([float(row["x_m"]) for row in camera_geometry])
    geometry_v = np.asarray([float(row["y_m"]) for row in camera_geometry])
    pixel_size = max(float(row["size_m"]) for row in camera_geometry)
    physical_limits = (
        (geometry_u.min() - pixel_size, geometry_u.max() + pixel_size),
        (geometry_v.min() - pixel_size, geometry_v.max() + pixel_size),
    )
    low = output_uv.min(axis=0)
    high = output_uv.max(axis=0)
    center = (low + high) / 2.0
    half = max(high - low) * 0.58
    half = max(half, pixel_size * 2.5)
    spot_limits = ((center[0] - half, center[0] + half), (center[1] - half, center[1] + half))

    fig = plt.figure(figsize=(14.2, 6.6), facecolor="white")
    grid = fig.add_gridspec(1, 3, width_ratios=(1.0, 1.0, 0.035), wspace=0.28)
    axes = [fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])]
    colorbar_axis = fig.add_subplot(grid[0, 2])
    collection = draw_camera_panel(
        axes[0], camera_geometry, signal, output_uv, camera_uv,
        physical_limits, "完整物理相机范围（1616像素）",
    )
    draw_camera_panel(
        axes[1], camera_geometry, signal, output_uv, camera_uv,
        spot_limits, "同一原始 u/v 光斑放大",
    )
    summary = case["summary"]
    centroid = case["camera_summary"]["output_uv_centroid_m"]
    fig.suptitle("71°平行光相机图：LACT_sim 原始 u/v", x=0.055, ha="left", y=0.985, color=INK, fontsize=15)
    fig.text(
        0.055, 0.94,
        f"全部 output={summary['hit_output_plane']:,}｜camera-hit={summary['hit_camera']:,}｜非零像素={summary['unique_hit_pixels']}｜质心=({centroid[0]:.6f}, {centroid[1]:.6f}) m｜未做 pyLAST 翻转",
        ha="left", color=MUTED, fontsize=9,
    )
    colorbar = fig.colorbar(collection, cax=colorbar_axis)
    colorbar.set_label("像素 photon_count 的平方根归一化", color=INK)
    colorbar.ax.tick_params(colors=MUTED, labelsize=8)
    fig.text(0.055, 0.025, "灰点：完整 hits CSV 的 u_m/v_m；黑点：同一 CSV 的 camera_x_m/camera_y_m；像素颜色：完整 camera CSV。", color=MUTED, fontsize=8)
    fig.subplots_adjust(left=0.055, right=0.96, top=0.88, bottom=0.11)
    fig.savefig(output, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", default="parallel_sky_up")
    parser.add_argument("--parallel-data", default="docs/assets/data/coordinate-parallel-cases.json")
    parser.add_argument("--ideal-mirror", default="configs/mirror_1229_facets.csv")
    parser.add_argument("--deformation-series", default="configs/mirror_1229_elevation_series.csv")
    parser.add_argument("--obstruction-primitives", default="configs/obstructions/raytrace_final_structure_primitives.csv")
    parser.add_argument("--camera-geometry", default="configs/cameras/new_camera_pixels.csv")
    parser.add_argument("--output-dir", default="docs/assets/plots")
    args = parser.parse_args()

    choose_font()
    payload = json.loads(Path(args.parallel_data).read_text(encoding="utf-8"))
    case = next((row for row in payload["four_direction_cases"] if row["id"] == args.case), None)
    if case is None:
        raise ValueError(f"parallel case not found: {args.case}")
    if case["validation"]["status"] != "passed":
        raise ValueError("refusing to plot an unvalidated parallel-light case")

    ideal, _, primitives = read_geometry(
        Path(args.ideal_mirror), Path(args.deformation_series), Path(args.obstruction_primitives)
    )
    camera_geometry = read_csv(Path(args.camera_geometry))
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    optics_path = output_dir / f"{args.case}-optics.png"
    camera_path = output_dir / f"{args.case}-camera.png"
    plot_optics(case, ideal, primitives, optics_path)
    plot_camera(case, camera_geometry, camera_path)
    print(f"Saved {optics_path}")
    print(f"Saved {camera_path}")


if __name__ == "__main__":
    main()
