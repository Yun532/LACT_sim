#!/usr/bin/env python3
"""Plot telescope array layout from a LACT_sim HDF5 trace file."""

import argparse
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm, Normalize
from matplotlib.ticker import MaxNLocator


def decode_attr(value):
    if isinstance(value, bytes):
        return value.decode()
    return str(value)


def event_table(h5):
    if "events/corsika" in h5:
        return h5["events/corsika"][:]
    if "events/table" in h5:
        return h5["events/table"][:]
    raise SystemExit("This HDF5 file has no events table.")


def event_id_from_shower_id(h5, shower_event_id, array_id):
    rows = event_table(h5)
    names = rows.dtype.names or ()
    if "shower_event_id" not in names or "array_id" not in names:
        raise SystemExit(
            "--shower-event-id requires events/corsika metadata. "
            "Use --event-id for this file."
        )
    matches = rows[
        (rows["shower_event_id"] == shower_event_id) & (rows["array_id"] == array_id)
    ]
    if len(matches) == 0:
        raise SystemExit(
            f"No event found for shower_event_id={shower_event_id} and array_id={array_id}."
        )
    if len(matches) > 1:
        raise SystemExit("Internal error: multiple event rows match shower/array selection.")
    return int(matches[0]["event_id"])


def event_id_from_shower_number(h5, shower_event_number, array_id):
    rows = event_table(h5)
    if "shower_event_id" in (rows.dtype.names or ()):
        shower_ids = []
        seen = set()
        for row in rows:
            shower = int(row["shower_event_id"])
            if shower not in seen:
                seen.add(shower)
                shower_ids.append(shower)
        if shower_event_number < 1 or shower_event_number > len(shower_ids):
            raise SystemExit(
                f"--shower-event-number must be 1..{len(shower_ids)} for this file."
            )
        shower = shower_ids[shower_event_number - 1]
        return event_id_from_shower_id(h5, shower, array_id)

    if shower_event_number < 1 or shower_event_number > len(rows):
        raise SystemExit(f"--shower-event-number must be 1..{len(rows)} for this file.")
    return int(rows[shower_event_number - 1]["event_id"])


def resolve_event_id(h5, event_id=None, shower_event_id=None, shower_event_number=None, array_id=0):
    choices = [
        event_id is not None,
        shower_event_id is not None,
        shower_event_number is not None,
    ]
    if sum(choices) > 1:
        raise SystemExit(
            "Use only one event selector: --event-id, --shower-event-id, "
            "or --shower-event-number."
        )
    if event_id is not None:
        return int(event_id)
    if shower_event_id is not None:
        return event_id_from_shower_id(h5, shower_event_id, array_id)
    if shower_event_number is not None:
        return event_id_from_shower_number(h5, shower_event_number, array_id)
    return None


def get_telescope_positions(telescopes):
    names = telescopes.dtype.names or ()
    if "array_x_north_m" in names and "array_y_west_m" in names:
        north = telescopes["array_x_north_m"].astype(float)
        west = telescopes["array_y_west_m"].astype(float)
    else:
        north = telescopes["x_m"].astype(float)
        west = telescopes["y_m"].astype(float)
    east = -west
    return east, north


def image_totals(h5, event_id, quantity):
    index = h5["images/index"][:]
    rows = index[index["event_id"] == event_id]
    totals = {}
    if len(rows) == 0:
        raise SystemExit(f"No images found for event_id={event_id}.")

    dense_path = f"images/dense/{quantity}"
    if dense_path in h5:
        dense = h5[dense_path]
        for row in rows:
            image_index = int(row["image_index"])
            totals[int(row["telescope_id"])] = float(np.sum(dense[image_index, :]))
        return totals

    if quantity not in ("pe", "signal", "photon_count"):
        raise SystemExit(f"{quantity} is not available in this HDF5 file.")
    for row in rows:
        if quantity == "pe":
            totals[int(row["telescope_id"])] = float(row["total_pe"])
        elif quantity == "signal":
            totals[int(row["telescope_id"])] = float(row["total_signal"])
        else:
            start = int(row["start"])
            count = int(row["count"])
            sparse = h5["images/sparse/pixels"][start:start + count]
            totals[int(row["telescope_id"])] = float(np.sum(sparse["photon_count"]))
    return totals


def get_event_metadata(h5, event_id):
    if "events/corsika" not in h5:
        return None
    rows = h5["events/corsika"][:]
    matches = rows[rows["event_id"] == event_id]
    if len(matches) == 0:
        return None
    return matches[0]


def nearby_core_rows(h5, n_cores, center_east, center_north):
    if n_cores <= 0:
        return []
    if "events/corsika_showers" in h5:
        rows = h5["events/corsika_showers"][:]
    elif "events/corsika" in h5:
        raw = h5["events/corsika"][:]
        selected = []
        seen = set()
        for row in raw:
            shower = int(row["shower_event_id"])
            if shower in seen:
                continue
            seen.add(shower)
            selected.append(row)
        rows = np.asarray(selected, dtype=raw.dtype)
    else:
        return []
    if len(rows) == 0:
        return []

    east = -rows["core_y_west_m"].astype(float)
    north = rows["core_x_north_m"].astype(float)
    distance2 = (east - center_east) ** 2 + (north - center_north) ** 2
    return [rows[i] for i in np.argsort(distance2)[:n_cores]]


def add_compass(ax, x0, y0, length):
    ax.annotate(
        "",
        xy=(x0, y0 + length),
        xytext=(x0, y0),
        arrowprops=dict(arrowstyle="-|>", lw=1.2, color="0.08"),
        zorder=8,
    )
    ax.text(x0, y0 + length * 1.12, "N", ha="center", va="bottom", fontsize=11, weight="bold")
    ax.text(x0, y0 - length * 0.18, "S", ha="center", va="top", fontsize=9, color="0.35")
    ax.annotate(
        "",
        xy=(x0 + length, y0),
        xytext=(x0, y0),
        arrowprops=dict(arrowstyle="-|>", lw=1.0, color="0.35"),
        zorder=8,
    )
    ax.text(x0 + length * 1.12, y0, "E", ha="left", va="center", fontsize=9, color="0.35")
    ax.text(x0 - length * 0.18, y0, "W", ha="right", va="center", fontsize=9, color="0.35")


def add_arrival_arrow(ax, core_east, core_north, az_deg, span, mode):
    az = np.deg2rad(az_deg)
    if mode == "incoming":
        dx = -np.sin(az)
        dy = -np.cos(az)
    else:
        dx = np.sin(az)
        dy = np.cos(az)
    length = 0.20 * span
    ax.annotate(
        "",
        xy=(core_east + dx * length, core_north + dy * length),
        xytext=(core_east, core_north),
        arrowprops=dict(arrowstyle="-|>", lw=1.5, color="#b2182b"),
    )


def format_event_info(event_id, core_x, core_y, arrival_az):
    lines = []
    if event_id is not None:
        lines.append(f"event_id = {event_id}")
    if core_x is not None and core_y is not None:
        lines.append(f"core: N = {core_x:.1f} m, E = {-core_y:.1f} m")
    if arrival_az is not None:
        lines.append(f"arrival azimuth = {arrival_az:.2f} deg")
    return "\n".join(lines)


def total_quantity_label(quantity):
    if quantity == "pe":
        return "Total p.e."
    if quantity == "photon_count":
        return "Total photon count"
    if quantity == "cherenkov_pe":
        return "Total Cherenkov p.e."
    if quantity == "nsb_pe":
        return "Total NSB p.e."
    return f"Total {quantity}"


def event_title(event_id, event_meta):
    if event_id is None:
        return "LACT telescope array layout"
    if event_meta is not None:
        return (
            f"LACT array event_id={event_id} "
            f"(shower={int(event_meta['shower_event_id'])}, array={int(event_meta['array_id'])})"
        )
    return f"LACT array event_id={event_id}"


def main():
    parser = argparse.ArgumentParser(
        description="Plot telescope array layout and optional event-level p.e. totals."
    )
    parser.add_argument("h5", help="HDF5 file from run_corsika_trace")
    parser.add_argument("--event-id", type=int, default=None)
    parser.add_argument(
        "--shower-event-number",
        "--event-number",
        dest="shower_event_number",
        type=int,
        default=None,
        help="1-based original CORSIKA shower-event order in the HDF5 file.",
    )
    parser.add_argument(
        "--shower-event-id",
        type=int,
        default=None,
        help="Original CORSIKA shower event id. Use with --array-id to select one output event.",
    )
    parser.add_argument(
        "--array-id",
        type=int,
        default=0,
        help=(
            "CORSIKA CSCAT/MC_TELOFF array-use index for the selected shower event; "
            "not the telescope ID."
        ),
    )
    parser.add_argument(
        "--quantity",
        choices=("pe", "signal", "photon_count", "cherenkov_pe", "nsb_pe"),
        default="pe",
    )
    parser.add_argument("--core-x-m", type=float, default=None, help="Core north coordinate [m].")
    parser.add_argument("--core-y-m", type=float, default=None, help="Core west coordinate [m].")
    parser.add_argument("--arrival-az-deg", type=float, default=None, help="Azimuth North-to-East [deg].")
    parser.add_argument(
        "--direction-mode",
        choices=("incoming", "source"),
        default="incoming",
        help="Draw horizontal projection of incoming shower direction or source direction.",
    )
    parser.add_argument("--no-labels", action="store_true", help="Do not label telescope IDs.")
    parser.add_argument("--label-prefix", default="T", help="Telescope label prefix.")
    parser.add_argument(
        "--cmap",
        default="inferno",
        help="Matplotlib colormap for event quantities. Default: inferno.",
    )
    parser.add_argument(
        "--vmax-percentile",
        type=float,
        default=98.0,
        help="Percentile used for the colorbar upper limit. Default: 98.",
    )
    parser.add_argument(
        "--log-color",
        action="store_true",
        help="Use logarithmic color normalization for event quantities.",
    )
    parser.add_argument(
        "--linear-color",
        action="store_true",
        help="Use linear color normalization. By default, p.e. maps use log color.",
    )
    parser.add_argument(
        "--show-nearby-cores",
        type=int,
        default=0,
        help="Mark N CORSIKA shower cores closest to the array center, if metadata are available.",
    )
    parser.add_argument(
        "--array-only-limits",
        action="store_true",
        help="Keep axis limits around telescopes only, even when cores are outside the array.",
    )
    parser.add_argument("--title", default=None)
    parser.add_argument("--output", default="array_layout.png")
    parser.add_argument("--dpi", type=int, default=300)
    args = parser.parse_args()

    with h5py.File(args.h5, "r") as h5:
        if "telescopes/table" not in h5:
            raise SystemExit("This HDF5 file has no telescopes/table dataset.")
        telescopes = h5["telescopes/table"][:]
        tel_ids = telescopes["telescope_id"].astype(int)
        east, north = get_telescope_positions(telescopes)

        event_id = resolve_event_id(
            h5,
            event_id=args.event_id,
            shower_event_id=args.shower_event_id,
            shower_event_number=args.shower_event_number,
            array_id=args.array_id,
        )

        values = None
        event_meta = None
        if event_id is not None:
            totals = image_totals(h5, event_id, args.quantity)
            values = np.array([totals.get(int(tel), 0.0) for tel in tel_ids], dtype=float)
            event_meta = get_event_metadata(h5, event_id)

        core_x = args.core_x_m
        core_y = args.core_y_m
        arrival_az = args.arrival_az_deg
        if event_meta is not None:
            if core_x is None:
                core_x = float(event_meta["core_x_north_m"])
            if core_y is None:
                core_y = float(event_meta["core_y_west_m"])
            if arrival_az is None:
                arrival_az = float(event_meta["azimuth_north_to_east_deg"])

        attrs = {}
        if "metadata/coordinates" in h5:
            attrs = {k: decode_attr(v) for k, v in h5["metadata/coordinates"].attrs.items()}
        nearby_cores = nearby_core_rows(
            h5,
            args.show_nearby_cores,
            float(np.mean(east)),
            float(np.mean(north)),
        )

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.labelsize": 11,
        "axes.titlesize": 13,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "axes.linewidth": 0.8,
    })
    fig, ax = plt.subplots(figsize=(7.6, 7.0))
    span = max(float(np.ptp(east)), float(np.ptp(north)), 1.0)
    pad = 0.18 * span
    marker_size = 86

    if values is None:
        sc = ax.scatter(
            east,
            north,
            s=marker_size,
            c="white",
            edgecolor="0.08",
            linewidth=0.9,
            zorder=3,
        )
    else:
        positive = values[values > 0]
        vmax = float(np.percentile(positive, args.vmax_percentile)) if positive.size else 1.0
        vmax = max(vmax, float(positive.max()) if positive.size and positive.max() < vmax else vmax, 1.0)
        use_log_color = args.log_color or (args.quantity == "pe" and not args.linear_color)
        if use_log_color:
            zero_mask = values <= 0.0
            if np.any(zero_mask):
                ax.scatter(
                    east[zero_mask],
                    north[zero_mask],
                    s=marker_size,
                    facecolors="none",
                    edgecolors="0.25",
                    linewidth=0.9,
                    zorder=3,
                )
            if positive.size:
                vmin = max(float(np.min(positive)), 1.0e-12)
                norm = LogNorm(vmin=vmin, vmax=vmax)
                color_mask = values > 0.0
            else:
                norm = LogNorm(vmin=1.0, vmax=vmax)
                color_mask = np.zeros_like(values, dtype=bool)
        else:
            norm = Normalize(vmin=0.0, vmax=vmax)
            color_mask = np.ones_like(values, dtype=bool)
        sc = ax.scatter(
            east[color_mask],
            north[color_mask],
            s=marker_size,
            c=values[color_mask],
            cmap=args.cmap,
            norm=norm,
            edgecolor="0.08",
            linewidth=0.35,
            zorder=3,
        )
        cbar = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label(total_quantity_label(args.quantity), fontsize=11)
        cbar.ax.tick_params(labelsize=9, direction="in")
        if not use_log_color:
            cbar.ax.yaxis.set_major_locator(MaxNLocator(nbins=6))

    if not args.no_labels:
        for tel, x, y in zip(tel_ids, east, north):
            ax.annotate(
                f"{args.label_prefix}{int(tel) + 1}",
                xy=(x, y),
                xytext=(3.0, 3.0),
                textcoords="offset points",
                ha="left",
                va="bottom",
                fontsize=7.2,
                color="0.12",
            )

    extra_east = []
    extra_north = []

    if core_x is not None and core_y is not None:
        core_east = -core_y
        core_north = core_x
        extra_east.append(core_east)
        extra_north.append(core_north)
        ax.scatter(
            [core_east],
            [core_north],
            marker="*",
            s=230,
            c="#d73027",
            edgecolor="white",
            linewidth=0.8,
            zorder=6,
        )
        core_label = (
            "Core\n"
            f"N = {core_x:.1f} m\n"
            f"E = {-core_y:.1f} m"
        )
        x_offset = -12 if core_east > float(np.mean(east)) else 12
        y_offset = -12 if core_north > float(np.mean(north)) else 12
        core_xytext = (x_offset, y_offset)
        core_va = "top" if y_offset < 0 else "bottom"
        if x_offset < 0:
            core_ha = "right"
        else:
            core_ha = "left"
        ax.annotate(
            core_label,
            xy=(core_east, core_north),
            xytext=core_xytext,
            textcoords="offset points",
            ha=core_ha,
            va=core_va,
            fontsize=8.5,
            color="#7f0000",
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="#d73027", alpha=0.94),
            arrowprops=dict(arrowstyle="-", color="#d73027", lw=0.9),
            zorder=7,
        )
        if arrival_az is not None:
            add_arrival_arrow(ax, core_east, core_north, arrival_az, span, args.direction_mode)
    elif arrival_az is not None:
        add_arrival_arrow(ax, float(np.mean(east)), float(np.mean(north)), arrival_az, span, args.direction_mode)

    if nearby_cores:
        core_east_values = -np.asarray([float(row["core_y_west_m"]) for row in nearby_cores])
        core_north_values = np.asarray([float(row["core_x_north_m"]) for row in nearby_cores])
        extra_east.extend(core_east_values.tolist())
        extra_north.extend(core_north_values.tolist())
        ax.scatter(
            core_east_values,
            core_north_values,
            marker="o",
            s=58,
            facecolors="none",
            edgecolors="#2166ac",
            linewidth=1.15,
            label=f"{len(nearby_cores)} nearest cores",
            zorder=5,
        )
        for row, x, y in zip(nearby_cores, core_east_values, core_north_values):
            ax.annotate(
                f"{int(row['shower_event_id'])}",
                xy=(x, y),
                xytext=(4, 4),
                textcoords="offset points",
                fontsize=6.8,
                color="#2166ac",
                ha="left",
                va="bottom",
            )

    ax.set_aspect("equal", adjustable="box")
    if extra_east and not args.array_only_limits:
        limit_east = np.concatenate([east, np.asarray(extra_east, dtype=float)])
        limit_north = np.concatenate([north, np.asarray(extra_north, dtype=float)])
    else:
        limit_east = east
        limit_north = north
    ax.set_xlim(float(np.min(limit_east)) - pad, float(np.max(limit_east)) + pad)
    ax.set_ylim(float(np.min(limit_north)) - pad, float(np.max(limit_north)) + pad)
    ax.set_xlabel("East [m]")
    ax.set_ylabel("North [m]")
    ax.grid(True, alpha=0.25, linewidth=0.55)
    ax.tick_params(direction="in", top=True, right=True)

    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()
    axis_span = min(x_max - x_min, y_max - y_min)
    compass_length = 0.075 * axis_span
    compass_x = x_min + 0.070 * (x_max - x_min)
    compass_y = y_min + 0.085 * (y_max - y_min)
    add_compass(ax, compass_x, compass_y, compass_length)

    if args.title is not None:
        title = args.title
    else:
        title = event_title(event_id, event_meta)
    ax.set_title(title)

    if attrs:
        ax.text(
            0.01,
            0.99,
            "HDF5 positions: x=N, y=W; plotted as East vs North",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=8.2,
            color="0.25",
            bbox=dict(boxstyle="round,pad=0.24", facecolor="white", edgecolor="0.82", alpha=0.9),
        )
    info = format_event_info(event_id, core_x, core_y, arrival_az)
    if info:
        ax.text(
            0.99,
            0.01,
            info,
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=8.4,
            color="0.12",
            bbox=dict(boxstyle="round,pad=0.28", facecolor="white", edgecolor="0.55", alpha=0.94),
        )
    if nearby_cores:
        ax.legend(loc="upper right", frameon=True, framealpha=0.92)

    out = Path(args.output)
    if out.parent:
        out.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out, dpi=args.dpi, bbox_inches="tight")
    print(f"Saved {out}")
    if values is not None:
        top = sorted(zip(values, tel_ids), reverse=True)[:8]
        print("Top telescopes:", ", ".join(f"TEL {int(tel) + 1}={val:.3g}" for val, tel in top))


if __name__ == "__main__":
    main()
