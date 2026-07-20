#!/usr/bin/env python3
"""Plot the spectral NSB-rate calculation and a Poisson sampling check."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spectrum", default="configs/nsb/nsb_spectrum.csv")
    parser.add_argument("--mirror", default="configs/efficiency/mirror_reflectivity.csv")
    parser.add_argument("--filter", default="configs/efficiency/filter_transmission.csv")
    parser.add_argument("--sipm", default="configs/efficiency/sipm_pde.csv")
    parser.add_argument("--effective-area-m2", type=float, default=24.576860)
    parser.add_argument("--pixel-size-m", type=float, default=0.0244)
    parser.add_argument("--focal-length-m", type=float, default=8.0)
    parser.add_argument("--window-ns", type=float, default=25.0)
    parser.add_argument("--n-samples", type=int, default=200000)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--output", default="nsb_spectral_response.png")
    parser.add_argument("--diagnostic-csv", default=None)
    parser.add_argument("--summary", default=None)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    wave, lons = read_two_col(args.spectrum)
    wm, mirror = read_two_col(args.mirror)
    wf, filt = read_two_col(args.filter)
    ws, sipm = read_two_col(args.sipm)

    mirror_i = interp_to(wm, mirror, wave)
    filter_i = interp_to(wf, filt, wave)
    sipm_i = interp_to(ws, sipm, wave)
    total_eff = mirror_i * filter_i * sipm_i
    detected = lons * total_eff

    red_integral = float(np.trapz(detected, wave))
    omega = (args.pixel_size_m / args.focal_length_m) ** 2
    rate = 1.0e-9 * red_integral * args.effective_area_m2 * omega
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
                [wave, lons, mirror_i, filter_i, sipm_i, total_eff, detected]
            ),
            delimiter=",",
            header=(
                "wavelength_nm,lons_flux_ph_s_nm_sr_m2,mirror_reflectivity,"
                "filter_transmission,sipm_pde,total_efficiency,"
                "lons_times_eff_pe_s_nm_sr_m2"
            ),
            comments="",
        )
    if args.summary:
        summary = Path(args.summary)
        summary.parent.mkdir(parents=True, exist_ok=True)
        summary.write_text(
            "NSB spectral-rate summary\n"
            "=========================\n"
            f"spectrum: {args.spectrum}\n"
            f"effective_area_m2: {args.effective_area_m2:.9f}\n"
            f"pixel_solid_angle_sr: {omega:.12e}\n"
            f"spectral_integral_pe_s_sr_m2: {red_integral:.10e}\n"
            f"rate_pe_per_ns_per_pixel: {rate:.10e}\n"
            f"mean_pe_per_pixel_{args.window_ns:g}ns: {mean:.10e}\n"
            f"sample_mean_pe: {samples.mean():.10e}\n"
            f"sample_rms_pe: {samples.std(ddof=0):.10e}\n"
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
    fig = plt.figure(figsize=(12.5, 5.6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.45, 1.0], wspace=0.28)

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = ax1.twinx()
    ax1.plot(wave, lons, color="black", lw=2.0, label="LoNS")
    ax1.plot(wave, detected, color="tab:red", lw=2.2, label="LoNS x total efficiency")
    ax1.set_yscale("log")
    ax1.set_xlim(float(wave.min()), float(wave.max()))
    ax1.set_xlabel("Wavelength [nm]")
    ax1.set_ylabel(r"Flux [ph s$^{-1}$ nm$^{-1}$ sr$^{-1}$ m$^{-2}$]")
    ax1.grid(True, which="major", alpha=0.18, lw=0.6)
    ax2.plot(wave, mirror_i, color="tab:blue", lw=1.7, label="Mirror")
    ax2.plot(wave, filter_i, color="tab:green", lw=1.7, ls="--", label="Filter")
    ax2.plot(wave, sipm_i, color="tab:orange", lw=1.7, ls="-.", label="SiPM PDE")
    ax2.plot(wave, total_eff, color="tab:purple", lw=2.0, ls=":", label="Total efficiency")
    ax2.set_ylim(0, 1.05)
    ax2.set_ylabel("Efficiency")
    lines = ax1.lines + ax2.lines
    ax1.legend(lines, [line.get_label() for line in lines], loc="upper left", ncol=2)
    ax1.text(
        0.98,
        0.04,
        f"Aeff = {args.effective_area_m2:.2f} m2\n"
        f"Omega = {omega:.2e} sr\n"
        f"rate = {rate:.4f} pe/ns/pixel",
        transform=ax1.transAxes,
        ha="right",
        va="bottom",
        bbox=dict(facecolor="white", edgecolor="0.75", alpha=0.95),
    )

    axh = fig.add_subplot(gs[0, 1])
    max_count = int(max(samples.max(), np.ceil(mean + 6.0 * np.sqrt(max(mean, 1e-12)))))
    bins = np.arange(-0.5, max_count + 1.5, 1.0)
    axh.hist(samples, bins=bins, density=True, alpha=0.55, color="tab:blue", label="simulation")
    xs = np.arange(0, max_count + 1)
    probs = np.exp(-mean) * np.ones_like(xs, dtype=float)
    for i in range(1, len(xs)):
        probs[i] = probs[i - 1] * mean / i
    axh.plot(xs, probs, "o-", color="black", ms=3.5, lw=1.2, label="Poisson theory")
    axh.set_xlabel(f"NSB p.e. / pixel / {args.window_ns:g} ns")
    axh.set_ylabel("Probability")
    axh.legend(frameon=True)
    axh.text(
        0.98,
        0.96,
        f"mean = {mean:.3f}\nsample = {samples.mean():.3f}",
        transform=axh.transAxes,
        ha="right",
        va="top",
        bbox=dict(facecolor="white", edgecolor="0.75", alpha=0.95),
    )

    fig.savefig(out, dpi=350, bbox_inches="tight")
    print(f"rate_pe_per_ns_per_pixel={rate:.10e}")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
