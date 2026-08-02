#!/usr/bin/env python3
"""Quick-look plots for LACT_sim lact_event_root_v1 files.

This script uses uproot instead of PyROOT so it can run in a lightweight
plotting environment after the ROOT file has been produced.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import numpy as np


def require_modules() -> tuple[Any, Any, Any]:
    try:
        import awkward as ak
        import matplotlib.pyplot as plt
        import uproot
    except ImportError as exc:
        raise SystemExit(
            "Missing plotting dependency. Install with:\n"
            "  python -m pip install uproot awkward matplotlib numpy\n"
            "or:\n"
            "  mamba install -c conda-forge uproot matplotlib numpy\n"
        ) from exc
    return ak, plt, uproot


def as_numpy(ak: Any, value: Any) -> np.ndarray:
    return np.asarray(ak.to_numpy(value))


def get_scalar_column(arrays: Any, name: str, fallback: Any = None) -> Any:
    if name not in arrays.fields:
        return fallback
    return arrays[name]


def choose_observation(ak: Any, observations: Any, event_id: int | None, tel_id: int | None) -> int:
    event_ids = as_numpy(ak, observations["event_id"])
    tel_ids = as_numpy(ak, observations["telescope_id"])
    mask = np.ones(len(event_ids), dtype=bool)
    if event_id is not None:
        mask &= event_ids == event_id
    if tel_id is not None:
        mask &= tel_ids == tel_id
    matches = np.nonzero(mask)[0]
    if len(matches) == 0:
        raise SystemExit(
            f"No observation found for event_id={event_id!r}, telescope_id={tel_id!r}"
        )
    return int(matches[0])


def dense_camera_values(
    ak: Any,
    camera_pixel_ids: np.ndarray,
    sparse_pixel_ids: Any,
    sparse_values: Any,
    fill: float,
) -> np.ndarray:
    dense = np.full(len(camera_pixel_ids), fill, dtype=float)
    lookup = {int(pixel_id): i for i, pixel_id in enumerate(camera_pixel_ids)}
    pixel_ids = as_numpy(ak, sparse_pixel_ids).astype(int)
    values = as_numpy(ak, sparse_values).astype(float)
    for pixel_id, value in zip(pixel_ids, values):
        idx = lookup.get(int(pixel_id))
        if idx is not None:
            dense[idx] = value
    return dense


def scatter_camera(plt: Any, x: np.ndarray, y: np.ndarray, values: np.ndarray,
                   title: str, colorbar_label: str, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 6.2), constrained_layout=True)
    finite = np.isfinite(values)
    plot_values = np.where(finite, values, np.nan)
    sc = ax.scatter(x, y, c=plot_values, s=18, cmap="viridis", marker="h", linewidths=0)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("LACT focal-plane u [m]")
    ax.set_ylabel("LACT focal-plane v [m]")
    ax.set_title(title)
    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label(colorbar_label)
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_total_pe_hist(ak: Any, plt: Any, observations: Any, path: Path) -> None:
    total_pe = as_numpy(ak, observations["total_pe"]).astype(float)
    fig, ax = plt.subplots(figsize=(7.2, 4.8), constrained_layout=True)
    ax.hist(total_pe[np.isfinite(total_pe)], bins=40, histtype="stepfilled", alpha=0.75)
    ax.set_xlabel("total p.e.")
    ax.set_ylabel("event-telescope count")
    ax.set_title("Observation total p.e.")
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_truth(plt: Any, corsika_events: Any, path: Path) -> None:
    if corsika_events is None or len(corsika_events["event_id"]) == 0:
        return
    energy_gev = np.asarray(corsika_events["energy_gev"], dtype=float)
    altitude_deg = np.asarray(corsika_events["altitude_deg"], dtype=float)
    x_max = np.asarray(corsika_events["x_max_g_cm2"], dtype=float)

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2), constrained_layout=True)
    axes[0].hist(energy_gev[np.isfinite(energy_gev)] / 1000.0, bins=40, histtype="stepfilled")
    axes[0].set_xlabel("energy [TeV]")
    axes[0].set_ylabel("events")
    axes[0].set_title("Primary energy")
    axes[1].hist(altitude_deg[np.isfinite(altitude_deg)], bins=40, histtype="stepfilled")
    axes[1].set_xlabel("altitude [deg]")
    axes[1].set_title("Altitude")
    axes[2].hist(x_max[np.isfinite(x_max)], bins=40, histtype="stepfilled")
    axes[2].set_xlabel("Xmax [g/cm^2]")
    axes[2].set_title("Xmax")
    fig.savefig(path, dpi=160)
    plt.close(fig)


def waveform_row_index(ak: Any, waveforms: Any, event_id: int, tel_id: int) -> int | None:
    if waveforms is None:
        return None
    event_ids = as_numpy(ak, waveforms["event_id"])
    tel_ids = as_numpy(ak, waveforms["telescope_id"])
    matches = np.nonzero((event_ids == event_id) & (tel_ids == tel_id))[0]
    if len(matches) == 0:
        return None
    return int(matches[0])


def plot_waveform_sum(ak: Any, plt: Any, root_file: Any, event_id: int, tel_id: int,
                      image_pe: np.ndarray, path: Path, summary_lines: list[str]) -> None:
    if "waveforms" not in root_file:
        summary_lines.append("waveforms: missing")
        return
    waveforms = root_file["waveforms"].arrays(library="ak")
    row = waveform_row_index(ak, waveforms, event_id, tel_id)
    if row is None:
        summary_lines.append(f"waveforms: no row for event={event_id}, tel={tel_id}")
        return

    time_bin = as_numpy(ak, waveforms["time_bin"][row]).astype(int)
    sample_value = as_numpy(ak, waveforms["sample_value"][row]).astype(float)
    pixel_id = as_numpy(ak, waveforms["pixel_id"][row]).astype(int)
    n_bins = int(waveforms["n_time_bins"][row])
    summed = np.bincount(time_bin, weights=sample_value, minlength=n_bins)
    x = np.arange(n_bins, dtype=float)
    if "waveform_config" in root_file:
        cfg = root_file["waveform_config"].arrays(library="ak")
        if len(cfg["time_centers_ns"]) > 0:
            centers = as_numpy(ak, cfg["time_centers_ns"][0])
            if len(centers) == n_bins:
                x = centers

    fig, ax = plt.subplots(figsize=(7.8, 4.8), constrained_layout=True)
    ax.step(x, summed, where="mid")
    ax.set_xlabel("time [ns]" if len(x) == n_bins else "time bin")
    ax.set_ylabel("summed sample value")
    ax.set_title(f"Camera-summed waveform, event {event_id}, tel {tel_id}")
    fig.savefig(path, dpi=160)
    plt.close(fig)

    wf_image_by_pixel: dict[int, float] = {}
    for pix, value in zip(pixel_id, sample_value):
        wf_image_by_pixel[int(pix)] = wf_image_by_pixel.get(int(pix), 0.0) + float(value)
    summary_lines.append(
        f"waveforms: nonzero samples={len(sample_value)}, "
        f"sum_sample_value={float(np.sum(sample_value)):.6g}"
    )
    summary_lines.append(f"selected image total_pe={float(np.sum(image_pe)):.6g}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root_file", type=Path, help="LACT_sim lact_events.root")
    parser.add_argument("--outdir", type=Path, default=None, help="Directory for PNG outputs")
    parser.add_argument("--event-id", type=int, default=None, help="Event id to plot")
    parser.add_argument("--tel-id", type=int, default=None, help="Telescope id to plot")
    args = parser.parse_args()

    ak, plt, uproot = require_modules()
    outdir = args.outdir or args.root_file.with_suffix("").with_name(args.root_file.stem + "_plots")
    outdir.mkdir(parents=True, exist_ok=True)

    with uproot.open(args.root_file) as root_file:
        for required in ("config", "camera_pixels", "observations"):
            if required not in root_file:
                raise SystemExit(f"{args.root_file} is missing required tree '{required}'")

        config = root_file["config"].arrays(library="np")
        camera = root_file["camera_pixels"].arrays(library="np")
        observations = root_file["observations"].arrays(library="ak")
        corsika = root_file["corsika_events"].arrays(library="np") if "corsika_events" in root_file else None

        obs_idx = choose_observation(ak, observations, args.event_id, args.tel_id)
        event_id = int(observations["event_id"][obs_idx])
        tel_id = int(observations["telescope_id"][obs_idx])

        camera_pixel_ids = np.asarray(camera["pixel_id"], dtype=int)
        x = np.asarray(camera["x_m"], dtype=float)
        y = np.asarray(camera["y_m"], dtype=float)
        image_pe = dense_camera_values(
            ak,
            camera_pixel_ids,
            observations["pixel_id"][obs_idx],
            observations["image_pe"][obs_idx],
            0.0,
        )
        peak_time = dense_camera_values(
            ak,
            camera_pixel_ids,
            observations["pixel_id"][obs_idx],
            observations["image_time_peak_ns"][obs_idx],
            np.nan,
        )

        schema = str(config["schema_name"][0]) if "schema_name" in config else "unknown"
        profile = str(config["profile"][0]) if "profile" in config else "unknown"
        summary_lines = [
            f"file={args.root_file}",
            f"schema={schema}",
            f"profile={profile}",
            f"selected_event_id={event_id}",
            f"selected_telescope_id={tel_id}",
            f"camera_pixels={len(camera_pixel_ids)}",
            f"observations={len(observations['event_id'])}",
            f"selected_nonzero_pixels={len(observations['pixel_id'][obs_idx])}",
            f"selected_image_sum_pe={float(np.sum(image_pe)):.6g}",
        ]

        scatter_camera(
            plt,
            x,
            y,
            image_pe,
            f"Integrated p.e. image, event {event_id}, tel {tel_id}",
            "p.e.",
            outdir / f"event_{event_id}_tel_{tel_id}_image_pe.png",
        )
        scatter_camera(
            plt,
            x,
            y,
            peak_time,
            f"Peak time, event {event_id}, tel {tel_id}",
            "time [ns]",
            outdir / f"event_{event_id}_tel_{tel_id}_peak_time.png",
        )
        plot_waveform_sum(
            ak,
            plt,
            root_file,
            event_id,
            tel_id,
            image_pe,
            outdir / f"event_{event_id}_tel_{tel_id}_waveform_sum.png",
            summary_lines,
        )
        plot_total_pe_hist(ak, plt, observations, outdir / "total_pe_hist.png")
        plot_truth(plt, corsika, outdir / "corsika_truth_hist.png")

    summary_path = outdir / "summary.txt"
    summary_path.write_text("\n".join(summary_lines) + "\n", encoding="utf-8")
    print(f"Wrote quick-look plots to {outdir}")
    print(f"Wrote summary to {summary_path}")


if __name__ == "__main__":
    main()
