#!/usr/bin/env python3
"""Compare CORSIKA mirror-hit outcomes for two EventIO local-origin shifts."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D

from config_io import expand_component_config, load_facets_from_config
from plot_optical_layout_3d import aperture_polygon


STATUS_STYLE = {
    "reflected_to_output": ("#2166ac", ".", "Reflected to output plane"),
    "blocked_incoming": ("#e66101", "x", "Blocked before reflection"),
    "blocked_reflected": ("#b2182b", "+", "Blocked after reflection"),
    "reflected_missed_output": ("#878787", ".", "Reflected, missed output plane"),
}


def local_facet_outlines(config_path: Path) -> list[np.ndarray]:
    config, _ = expand_component_config(config_path)
    facets = load_facets_from_config(config_path, config)
    outlines = []
    for facet in facets:
        polygon = aperture_polygon(facet)[:, :2] * 1000.0
        outlines.append(np.vstack([polygon, polygon[0]]))
    return outlines


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("minus16_csv")
    parser.add_argument("zero_csv")
    parser.add_argument("--config", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--dpi", type=int, default=260)
    parser.add_argument("--point-size", type=float, default=1.0)
    parser.add_argument("--event-id", type=int)
    parser.add_argument("--telescope-id", type=int)
    return parser.parse_args()


def load_hits(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path)
    required = {
        "mirror_x_m",
        "mirror_y_m",
        "mirror_z_m",
        "status",
    }
    missing = required - set(frame.columns)
    if missing:
        raise SystemExit(f"{path} is missing: {', '.join(sorted(missing))}")
    unknown = set(frame["status"].unique()) - set(STATUS_STYLE)
    if unknown:
        raise SystemExit(f"{path} has unknown status values: {sorted(unknown)}")
    # The tracer and imported mirror layout both use telescope-local
    # coordinates here. Applying the telescope pointing transform again would
    # rotate and translate the hits a second time.
    frame["display_x_mm"] = frame["mirror_x_m"].to_numpy(dtype=float) * 1000.0
    frame["display_y_mm"] = frame["mirror_y_m"].to_numpy(dtype=float) * 1000.0
    return frame


def draw_panel(
    ax: plt.Axes,
    frame: pd.DataFrame,
    title: str,
    outlines: list[np.ndarray],
    point_size: float,
) -> None:
    # Draw the numerous successful/reflected points first. The two obstruction
    # stages then remain legible on top with distinct colours and marker shapes.
    for status in (
        "reflected_to_output",
        "reflected_missed_output",
        "blocked_incoming",
        "blocked_reflected",
    ):
        selected = frame[frame["status"] == status]
        if selected.empty:
            continue
        colour, marker, _ = STATUS_STYLE[status]
        is_obstruction = status.startswith("blocked_")
        ax.scatter(
            selected["display_x_mm"],
            selected["display_y_mm"],
            s=point_size * (4.0 if is_obstruction else 1.0),
            c=colour,
            marker=marker,
            alpha=0.55 if is_obstruction else 0.24,
            linewidths=0.28 if is_obstruction else 0,
            rasterized=True,
            zorder=4 if is_obstruction else 2,
        )

    ax.add_collection(
        LineCollection(outlines, colors="black", linewidths=0.38, alpha=0.55, zorder=5)
    )
    counts = frame["status"].value_counts()
    total = len(frame)
    note_lines = [f"Mirror hits: {total:,}"]
    for status in (
        "reflected_to_output",
        "blocked_incoming",
        "blocked_reflected",
        "reflected_missed_output",
    ):
        count = int(counts.get(status, 0))
        note_lines.append(f"{STATUS_STYLE[status][2]}: {count:,} ({count / total:.1%})")
    ax.text(
        0.02,
        0.02,
        "\n".join(note_lines),
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=7.6,
        bbox=dict(
            facecolor="white",
            edgecolor="0.7",
            alpha=0.92,
            boxstyle="round,pad=0.28",
        ),
        zorder=8,
    )
    ax.set_title(title)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(alpha=0.12, linewidth=0.5)
    ax.set_xlabel("Mirror-local x [mm]")


def main() -> None:
    args = parse_args()
    config = Path(args.config).resolve()
    minus16 = load_hits(Path(args.minus16_csv))
    zero = load_hits(Path(args.zero_csv))
    if args.event_id is not None:
        minus16 = minus16[minus16["event_id"] == args.event_id].copy()
        zero = zero[zero["event_id"] == args.event_id].copy()
    if args.telescope_id is not None:
        minus16 = minus16[minus16["telescope_id"] == args.telescope_id].copy()
        zero = zero[zero["telescope_id"] == args.telescope_id].copy()
    if minus16.empty or zero.empty:
        raise SystemExit("the requested event/telescope selection is empty")
    outlines = local_facet_outlines(config)

    all_x = np.concatenate(
        [minus16["display_x_mm"].to_numpy(), zero["display_x_mm"].to_numpy()]
    )
    all_y = np.concatenate(
        [minus16["display_y_mm"].to_numpy(), zero["display_y_mm"].to_numpy()]
    )
    max_abs = max(float(np.max(np.abs(all_x))), float(np.max(np.abs(all_y))))
    limit = max_abs * 1.025

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(13.2, 6.25),
        dpi=args.dpi,
        sharex=True,
        sharey=True,
    )
    draw_panel(
        axes[0],
        minus16,
        "EventIO local optical-origin shift = -16 m (default)",
        outlines,
        args.point_size,
    )
    draw_panel(
        axes[1],
        zero,
        "EventIO local optical-origin shift = 0 m (diagnostic)",
        outlines,
        args.point_size,
    )
    axes[0].set_ylabel("Mirror-local y [mm]")
    for ax in axes:
        ax.set_xlim(-limit, limit)
        ax.set_ylim(-limit, limit)

    legend_handles = []
    for status in (
        "reflected_to_output",
        "blocked_incoming",
        "blocked_reflected",
        "reflected_missed_output",
    ):
        colour, marker, label = STATUS_STYLE[status]
        legend_handles.append(
            Line2D(
                [0],
                [0],
                marker=marker,
                linestyle="none",
                color=colour,
                markersize=6,
                markeredgewidth=1.0,
                label=label,
            )
        )
    fig.legend(
        handles=legend_handles,
        loc="upper center",
        ncol=4,
        frameon=True,
        framealpha=0.97,
        bbox_to_anchor=(0.5, 0.985),
    )
    selection = []
    if args.event_id is not None:
        selection.append(f"event {args.event_id}")
    if args.telescope_id is not None:
        selection.append(f"telescope {args.telescope_id}")
    suffix = f" ({', '.join(selection)})" if selection else ""
    fig.suptitle(
        "CORSIKA mirror-plane outcomes for two local optical-origin shifts" + suffix,
        y=1.035,
        fontsize=14,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved CORSIKA mirror z comparison = {output}")


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    main()
