#!/usr/bin/env python3
"""Plot LACT_sim spectral inputs and their detected responses."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

try:
    from config_io import expand_component_config, resolve_workspace_path
except ModuleNotFoundError:  # Support ``python -m python.plot_nsb_spectral_rate``.
    from .config_io import expand_component_config, resolve_workspace_path


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = REPO_ROOT / "configs/nsb/spectral_rate_check_with_obstruction.cfg"
DEFAULT_OUTPUT = (
    REPO_ROOT
    / "run_logs/manual_checks/nsb_spectral/lons_instrument_cherenkov_30ns.png"
)


def read_two_col(path: str) -> tuple[np.ndarray, np.ndarray]:
    rows: list[tuple[float, float]] = []
    for line in Path(path).read_text().splitlines():
        text = line.strip()
        if not text or text.startswith("#") or text.lower().startswith("wavelength"):
            continue
        cells = text.replace(",", " ").split()
        if len(cells) < 2:
            continue
        rows.append((float(cells[0]), float(cells[1])))
    if len(rows) < 2:
        raise ValueError(f"{path} needs at least two numeric rows")
    arr = np.asarray(rows, dtype=float)
    order = np.argsort(arr[:, 0])
    return arr[order, 0], arr[order, 1]


def interp_to(x: np.ndarray, y: np.ndarray, xnew: np.ndarray) -> np.ndarray:
    return np.interp(xnew, x, y, left=0.0, right=0.0)


def input_path(value: str | Path) -> Path:
    path = Path(value)
    return path if path.is_absolute() else REPO_ROOT / path


def first_camera_pixel_size(path: Path) -> float:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        row = next(csv.DictReader(stream), None)
    if row is None or "size_m" not in row:
        raise ValueError(f"camera CSV has no size_m data: {path}")
    return float(row["size_m"])


def apply_config_defaults(args: argparse.Namespace) -> None:
    config_path = input_path(args.config)
    cfg, _ = expand_component_config(config_path)

    def configured_path(key: str) -> Path:
        value = cfg.get(key)
        if not value:
            raise ValueError(f"{key} is missing from {config_path}")
        return resolve_workspace_path(config_path, value)

    args.config = config_path
    args.spectrum = input_path(args.spectrum) if args.spectrum else configured_path(
        "nsb.spectrum_csv"
    )
    args.mirror = input_path(args.mirror) if args.mirror else configured_path(
        "efficiency.mirror_reflectivity"
    )
    args.filter = input_path(args.filter) if args.filter else configured_path(
        "efficiency.filter_transmission"
    )
    args.sipm = input_path(args.sipm) if args.sipm else configured_path("sipm.pde")
    args.atmosphere_transmission = input_path(args.atmosphere_transmission)

    if args.effective_area_m2 is None:
        args.effective_area_m2 = float(cfg["nsb.effective_area_m2"])
    if args.collector_mean_transmission is None:
        args.collector_mean_transmission = float(
            cfg["nsb.collector_mean_transmission"]
        )
    if args.focal_length_m is None:
        args.focal_length_m = float(cfg["telescope.focal_length_m"])
    if args.pixel_size_m is None:
        camera_csv = configured_path("camera.csv_path")
        args.pixel_size_m = first_camera_pixel_size(camera_csv)


def cherenkov_spectrum(
    wave: np.ndarray,
    zenith_deg: float,
    transmission_path: str,
    trapezoid,
) -> np.ndarray:
    lines = Path(transmission_path).read_text().splitlines()
    header = next(line for line in lines if line.startswith("# H2="))
    left, right = header[1:].split("H1=")
    observer_km = float(left.split("H2=")[1].strip().rstrip(","))
    heights = np.asarray([observer_km, *map(float, right.split())])
    rows = np.asarray(
        [[float(x) for x in line.split()] for line in lines if line.strip() and not line.startswith("#")]
    )
    table_wave = rows[:, 0]
    tau = np.column_stack([np.zeros(len(rows)), rows[:, 1:]])

    cosz = np.cos(np.deg2rad(zenith_deg))
    emission_height = np.linspace(observer_km + 0.001, 30.0, 1600)
    vertical_depth = 1030.0 * np.exp(-emission_height / 8.4)
    slant_depth = vertical_depth / cosz
    shower_y = np.log(1.0e4 / 0.085)  # 10 TeV gamma, Ec = 85 MeV
    shower_t = np.maximum(slant_depth / 37.0, 1.0e-6)
    shower_age = 3.0 * shower_t / (shower_t + 2.0 * shower_y)
    height_weight = (
        0.31
        / np.sqrt(shower_y)
        * np.exp(shower_t * (1.0 - 1.5 * np.log(shower_age)))
        * vertical_depth
        / (8.4 * cosz)
    )

    spectrum = []
    for wavelength in wave:
        table_wavelength = np.clip(wavelength, table_wave[0], table_wave[-1])
        j = np.clip(
            np.searchsorted(table_wave, table_wavelength), 1, len(table_wave) - 1
        )
        f = (table_wavelength - table_wave[j - 1]) / (
            table_wave[j] - table_wave[j - 1]
        )
        tau_row = np.minimum((1.0 - f) * tau[j - 1] + f * tau[j], 700.0)
        tau_height = np.interp(
            np.log10(emission_height), np.log10(heights), tau_row
        )
        transmission = np.exp(-np.minimum(tau_height / cosz, 700.0))
        spectrum.append(
            wavelength**-2
            * trapezoid(height_weight * transmission, emission_height)
        )
    spectrum = np.asarray(spectrum)
    return spectrum / spectrum.max()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(DEFAULT_CONFIG))
    parser.add_argument("--spectrum", default=None)
    parser.add_argument(
        "--mirror",
        default=None,
    )
    parser.add_argument("--filter", default=None)
    parser.add_argument("--sipm", default=None)
    parser.add_argument(
        "--atmosphere-transmission",
        default=str(
            REPO_ROOT / "configs/atmosphere/atm_trans_4400_1_10_0_0_4400.dat"
        ),
    )
    parser.add_argument("--effective-area-m2", type=float, default=None)
    parser.add_argument(
        "--collector-mean-transmission", type=float, default=None
    )
    parser.add_argument("--pixel-size-m", type=float, default=None)
    parser.add_argument("--focal-length-m", type=float, default=None)
    parser.add_argument("--window-ns", type=float, default=30.0)
    parser.add_argument("--n-samples", type=int, default=200000)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--output", default=str(DEFAULT_OUTPUT))
    parser.add_argument("--diagnostic-csv", default=None)
    parser.add_argument("--summary", default=None)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    apply_config_defaults(args)
    wave, lons = read_two_col(args.spectrum)
    wm, mirror = read_two_col(args.mirror)
    wf, filt = read_two_col(args.filter)
    ws, sipm = read_two_col(args.sipm)

    mirror_i = interp_to(wm, mirror, wave)
    filter_i = interp_to(wf, filt, wave)
    sipm_i = interp_to(ws, sipm, wave)
    total_eff = mirror_i * filter_i * sipm_i
    detected = lons * total_eff

    instrument_wave = np.arange(260.0, 1000.1, 1.0)
    response_min = max(wm[0], wf[0], ws[0])
    response_max = min(wm[-1], wf[-1], ws[-1])
    response_wave = np.unique(
        np.concatenate(
            ([response_min], instrument_wave[(instrument_wave > response_min) &
                                              (instrument_wave < response_max)],
             [response_max])
        )
    )
    total_eff_plot = (
        np.interp(response_wave, wm, mirror)
        * np.interp(response_wave, wf, filt)
        * np.interp(response_wave, ws, sipm)
    )
    detected_plot = (wave >= response_min) & (wave <= response_max)

    # NumPy 2.x exposes the trapezoidal integrator as ``trapezoid``;
    # older supported environments only provide the deprecated ``trapz``.
    trapezoid = getattr(np, "trapezoid", None)
    if trapezoid is None:
        trapezoid = np.trapz
    cherenkov_wave = instrument_wave
    cherenkov_20_plot = cherenkov_spectrum(
        cherenkov_wave, 20.0, args.atmosphere_transmission, trapezoid
    )
    cherenkov_60_plot = cherenkov_spectrum(
        cherenkov_wave, 60.0, args.atmosphere_transmission, trapezoid
    )
    cherenkov_20_detected_plot = (
        np.interp(response_wave, cherenkov_wave, cherenkov_20_plot)
        * total_eff_plot
    )
    cherenkov_60_detected_plot = (
        np.interp(response_wave, cherenkov_wave, cherenkov_60_plot)
        * total_eff_plot
    )
    cherenkov_20 = interp_to(cherenkov_wave, cherenkov_20_plot, wave)
    cherenkov_60 = interp_to(cherenkov_wave, cherenkov_60_plot, wave)
    red_integral = float(trapezoid(detected, wave))
    omega = (args.pixel_size_m / args.focal_length_m) ** 2
    if not 0.0 < args.collector_mean_transmission <= 1.0:
        raise ValueError("--collector-mean-transmission must be in (0, 1]")
    rate = (
        1.0e-9
        * red_integral
        * args.effective_area_m2
        * omega
        * args.collector_mean_transmission
    )
    mean = rate * args.window_ns

    rng = np.random.default_rng(args.seed)
    samples = rng.poisson(mean, size=args.n_samples)

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    if args.diagnostic_csv:
        diag = Path(args.diagnostic_csv)
        diag.parent.mkdir(parents=True, exist_ok=True)
        np.savetxt(
            diag,
            np.column_stack(
                [
                    wave,
                    lons,
                    mirror_i,
                    filter_i,
                    sipm_i,
                    total_eff,
                    detected,
                    cherenkov_20,
                    cherenkov_60,
                ]
            ),
            delimiter=",",
            header=(
                "wavelength_nm,lons_flux_ph_s_nm_sr_m2,mirror_reflectivity,"
                "filter_transmission,sipm_pde,total_efficiency,"
                "lons_times_eff_pe_s_nm_sr_m2,"
                "cherenkov_20deg_peak_normalized,cherenkov_60deg_peak_normalized"
            ),
            comments="",
        )
    if args.summary:
        summary = Path(args.summary)
        summary.parent.mkdir(parents=True, exist_ok=True)
        summary.write_text(
            "NSB spectral-rate summary\n"
            "=========================\n"
            f"config: {args.config}\n"
            f"spectrum: {args.spectrum}\n"
            f"atmosphere_transmission: {args.atmosphere_transmission}\n"
            f"effective_area_m2: {args.effective_area_m2:.9f}\n"
            "collector_mean_transmission: "
            f"{args.collector_mean_transmission:.9f}\n"
            f"pixel_solid_angle_sr: {omega:.12e}\n"
            f"spectral_integral_pe_s_sr_m2: {red_integral:.10e}\n"
            f"rate_pe_per_ns_per_pixel: {rate:.10e}\n"
            f"mean_pe_per_pixel_{args.window_ns:g}ns: {mean:.10e}\n"
            f"sample_mean_pe: {samples.mean():.10e}\n"
            f"sample_rms_pe: {samples.std(ddof=0):.10e}\n"
            "cherenkov_wavelength_range_nm: 260-1000\n"
            f"cherenkov_peak_nm_20deg: {cherenkov_wave[np.argmax(cherenkov_20_plot)]:.1f}\n"
            f"cherenkov_peak_nm_60deg: {cherenkov_wave[np.argmax(cherenkov_60_plot)]:.1f}\n"
        )

    plt.rcParams.update(
        {
            "font.size": 12,
            "axes.labelsize": 13,
            "legend.fontsize": 9.5,
            "axes.linewidth": 1.1,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )
    fig, ax1 = plt.subplots(figsize=(9.0, 5.9))
    ax2 = ax1.twinx()
    ax1.set_zorder(1)
    ax2.set_zorder(2)
    ax1.patch.set_alpha(0.0)
    ax2.patch.set_alpha(0.0)

    ax1.plot(wave, lons, color="black", lw=2.4, zorder=10, label="LoNS")
    ax1.plot(
        wave[detected_plot],
        detected[detected_plot],
        color="tab:red",
        lw=2.6,
        label=r"LoNS $\times$ total efficiency",
    )
    ax1.set_xlabel("Wavelength [nm]")
    ax1.set_ylabel(r"Flux [ph nm$^{-1}$ s$^{-1}$ sr$^{-1}$ m$^{-2}$]")
    ax1.set_yscale("log")
    ax1.set_xlim(260.0, 900.0)
    ax1.set_xticks([260, *range(300, 901, 100)])
    ax1.set_ylim(1.0e9, 1.0e13)

    ax2.plot(
        wm,
        mirror,
        color="tab:blue",
        lw=2.0,
        label="Mirror reflectivity",
    )
    ax2.plot(
        wf,
        filt,
        color="tab:green",
        lw=2.0,
        ls="--",
        label="Filter transmission",
    )
    ax2.plot(
        ws,
        sipm,
        color="tab:orange",
        lw=2.0,
        ls="-.",
        label="SiPM PDE",
    )
    ax2.plot(
        response_wave,
        total_eff_plot,
        color="tab:purple",
        lw=2.4,
        ls=":",
        label="Total efficiency",
    )
    ax2.plot(
        cherenkov_wave,
        cherenkov_20_plot,
        color="#4a4a4a",
        lw=2.2,
        ls=(0, (6, 2)),
        label=r"Cherenkov $20^\circ$",
    )
    ax2.plot(
        response_wave,
        cherenkov_20_detected_plot,
        color="#d76a6a",
        lw=2.2,
        ls=(0, (6, 2)),
        label=r"Cherenkov $20^\circ$ $\times$ total efficiency",
    )
    ax2.plot(
        cherenkov_wave,
        cherenkov_60_plot,
        color="#9a9a9a",
        lw=2.2,
        ls=(0, (2, 1)),
        label=r"Cherenkov $60^\circ$",
    )
    ax2.plot(
        response_wave,
        cherenkov_60_detected_plot,
        color="#efaaaa",
        lw=2.2,
        ls=(0, (2, 1)),
        label=r"Cherenkov $60^\circ$ $\times$ total efficiency",
    )
    ax2.set_ylabel("Efficiency / normalized spectrum")
    ax2.set_ylim(0.0, 1.05)

    ax1.set_title("LoNS, instrumental response, and Cherenkov spectra")
    ax1.grid(True, which="major", alpha=0.18, lw=0.6)
    lines = {line.get_label(): line for line in ax1.lines + ax2.lines}
    legend_lines = [
        lines["Mirror reflectivity"],
        lines["Filter transmission"],
        lines[r"Cherenkov $20^\circ$"],
        lines[r"Cherenkov $60^\circ$"],
        lines["LoNS"],
        lines["SiPM PDE"],
        lines["Total efficiency"],
        lines[r"Cherenkov $20^\circ$ $\times$ total efficiency"],
        lines[r"Cherenkov $60^\circ$ $\times$ total efficiency"],
        lines[r"LoNS $\times$ total efficiency"],
    ]
    ax2.legend(
        legend_lines,
        [line.get_label() for line in legend_lines],
        loc="upper right",
        bbox_to_anchor=(0.99, 0.99),
        frameon=True,
        facecolor="white",
        edgecolor="0.7",
        framealpha=0.62,
        fancybox=False,
        fontsize=8.5,
        handlelength=2.6,
        handletextpad=0.5,
        columnspacing=1.2,
        labelspacing=0.3,
        borderpad=0.45,
        markerfirst=True,
        ncol=2,
    )
    fig.tight_layout()
    fig.savefig(out, dpi=350, bbox_inches="tight")
    print(f"rate_pe_per_ns_per_pixel={rate:.10e}")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
