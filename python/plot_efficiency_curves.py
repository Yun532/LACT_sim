#!/usr/bin/env python3
"""Plot and validate wavelength-dependent efficiency curves."""

from __future__ import annotations

import argparse
import bisect
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


@dataclass
class Curve:
    name: str
    path: Path | None
    wavelengths: np.ndarray
    values: np.ndarray
    raw_wavelengths: np.ndarray
    raw_values: np.ndarray
    duplicate_groups: int = 0
    duplicate_rows: int = 0
    out_of_range_rows: int = 0

    @property
    def enabled(self) -> bool:
        return self.path is not None

    def evaluate(self, grid: np.ndarray) -> np.ndarray:
        if not self.enabled:
            return np.ones_like(grid, dtype=float)
        return np.array([self.evaluate_one(float(wl)) for wl in grid], dtype=float)

    def evaluate_one(self, wavelength_nm: float) -> float:
        if self.wavelengths.size == 0:
            return 1.0
        xs = self.wavelengths.tolist()
        if wavelength_nm < xs[0] or wavelength_nm > xs[-1]:
            return 0.0
        idx = bisect.bisect_left(xs, wavelength_nm)
        if idx == 0:
            return float(self.values[0])
        if idx == len(xs):
            return float(self.values[-1])
        if xs[idx] == wavelength_nm:
            return float(self.values[idx])
        x0 = float(self.wavelengths[idx - 1])
        x1 = float(self.wavelengths[idx])
        y0 = float(self.values[idx - 1])
        y1 = float(self.values[idx])
        return y0 + (wavelength_nm - x0) / (x1 - x0) * (y1 - y0)

    def raw_sorted(self) -> tuple[np.ndarray, np.ndarray]:
        if not self.enabled:
            return np.array([], dtype=float), np.array([], dtype=float)
        order = np.argsort(self.raw_wavelengths, kind="stable")
        return self.raw_wavelengths[order], self.raw_values[order]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate and plot mirror/filter/SiPM/atmosphere efficiency curves."
    )
    parser.add_argument("--mirror-csv", default="configs/efficiency/mirror_reflectivity.csv")
    parser.add_argument("--filter-csv", default="configs/efficiency/filter_transmission.csv")
    parser.add_argument("--sipm-csv", default="configs/efficiency/sipm_pde.csv")
    parser.add_argument("--atmosphere-csv", default="none")
    parser.add_argument("--constant-scale", type=float, default=1.0)
    parser.add_argument("--wavelength-min", type=float, default=200.0)
    parser.add_argument("--wavelength-max", type=float, default=700.0)
    parser.add_argument("--step", type=float, default=1.0)
    parser.add_argument("--output-dir", default="run_logs/official_tests/efficiency_curves")
    parser.add_argument("--dpi", type=int, default=350)
    return parser.parse_args()


def normalize_path(value: str) -> Path | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text or text.lower() in {"none", "off", "false", "1", "ideal"}:
        return None
    return Path(text)


def read_curve(name: str, value: str) -> Curve:
    path = normalize_path(value)
    if path is None:
        return Curve(
            name=name,
            path=None,
            wavelengths=np.array([], dtype=float),
            values=np.array([], dtype=float),
            raw_wavelengths=np.array([], dtype=float),
            raw_values=np.array([], dtype=float),
        )
    if not path.exists():
        raise FileNotFoundError(path)

    rows: list[tuple[float, float]] = []
    with path.open() as handle:
        for line in handle:
            line = line.split("#", 1)[0].strip()
            if not line:
                continue
            parts = line.replace(",", " ").split()
            if len(parts) < 2:
                continue
            try:
                rows.append((float(parts[0]), float(parts[1])))
            except ValueError:
                continue
    if not rows:
        raise ValueError(f"{path} has no numeric rows")

    raw_wl = np.array([r[0] for r in rows], dtype=float)
    raw_val = np.array([r[1] for r in rows], dtype=float)
    out_of_range = int(np.count_nonzero((raw_val < 0.0) | (raw_val > 1.0)))

    rows.sort(key=lambda item: item[0])
    merged_wl: list[float] = []
    merged_val: list[float] = []
    duplicate_groups = 0
    duplicate_rows = 0
    i = 0
    while i < len(rows):
        wl = rows[i][0]
        vals = []
        while i < len(rows) and rows[i][0] == wl:
            vals.append(rows[i][1])
            i += 1
        if len(vals) > 1:
            duplicate_groups += 1
            duplicate_rows += len(vals)
        merged_wl.append(wl)
        merged_val.append(float(np.mean(vals)))

    return Curve(
        name=name,
        path=path,
        wavelengths=np.array(merged_wl, dtype=float),
        values=np.array(merged_val, dtype=float),
        raw_wavelengths=raw_wl,
        raw_values=raw_val,
        duplicate_groups=duplicate_groups,
        duplicate_rows=duplicate_rows,
        out_of_range_rows=out_of_range,
    )


def plot_single(curve: Curve, grid: np.ndarray, output_dir: Path, dpi: int) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    if curve.enabled:
        raw_wl, raw_val = curve.raw_sorted()
        ax.plot(
            raw_wl,
            raw_val,
            color="#ff7f0e",
            lw=1.4,
            ls="--",
            alpha=0.85,
            label="input table line",
        )
        ax.scatter(
            curve.raw_wavelengths,
            curve.raw_values,
            s=11,
            color="#888888",
            alpha=0.45,
            linewidths=0,
            label="input table points",
        )
        ax.plot(grid, curve.evaluate(grid), color="#1f77b4", lw=2.0, label="used curve")
    else:
        ax.plot(grid, np.ones_like(grid), color="#1f77b4", lw=2.0, label="disabled -> 1")
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("Efficiency")
    ax.set_ylim(-0.04, 1.08)
    ax.set_title(curve.name)
    ax.grid(True, color="#d9d9d9", linewidth=0.7, alpha=0.8)
    ax.legend(frameon=True, framealpha=1.0, edgecolor="black")
    fig.tight_layout()
    fig.savefig(output_dir / f"{curve.name}.png", dpi=dpi)
    plt.close(fig)


def display_name(name: str) -> str:
    return {
        "mirror_reflectivity": "Mirror reflectivity",
        "filter_transmission": "Filter transmission",
        "sipm_pde": "SiPM PDE",
        "atmosphere_transmission": "Atmosphere transmission",
    }.get(name, name)


def plot_summary(
    curves: list[Curve],
    grid: np.ndarray,
    total: np.ndarray,
    constant_scale: float,
    output_dir: Path,
    dpi: int,
) -> None:
    del constant_scale
    plt.rcParams.update(
        {
            "font.size": 11,
            "axes.labelsize": 12,
            "axes.titlesize": 13,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "legend.fontsize": 9,
        }
    )
    fig, ax = plt.subplots(figsize=(11.2, 5.8))
    colors = {
        "mirror_reflectivity": "#0072B2",
        "filter_transmission": "#009E73",
        "sipm_pde": "#D55E00",
        "atmosphere_transmission": "#CC79A7",
    }

    for curve in curves:
        color = colors.get(curve.name, "#666666")
        label = display_name(curve.name)
        if curve.enabled:
            raw_wl, raw_val = curve.raw_sorted()
            ax.plot(
                raw_wl,
                raw_val,
                color=color,
                lw=1.05,
                ls=(0, (3.0, 2.0)),
                alpha=0.52,
                label=f"{label} input",
            )
        ax.plot(
            grid,
            curve.evaluate(grid),
            color=color,
            lw=2.0,
            alpha=0.95,
            label=f"{label} used",
        )

    ax.plot(grid, total, color="black", lw=2.8, label="Total used")
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("Efficiency")
    ax.set_xlim(float(grid[0]), float(grid[-1]))
    ax.set_ylim(-0.02, 1.04)
    ax.set_title("Wavelength-dependent efficiency validation")
    ax.grid(True, which="major", color="#d0d0d0", linewidth=0.65, alpha=0.75)
    ax.minorticks_on()
    ax.grid(True, which="minor", color="#eeeeee", linewidth=0.45, alpha=0.65)
    for spine in ax.spines.values():
        spine.set_linewidth(0.9)

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(
        handles,
        labels,
        loc="center left",
        bbox_to_anchor=(1.015, 0.5),
        frameon=False,
        borderaxespad=0.0,
        handlelength=2.8,
    )
    fig.tight_layout(rect=(0.0, 0.0, 0.80, 1.0))
    fig.savefig(output_dir / "efficiency_summary.png", dpi=dpi)
    fig.savefig(output_dir / "efficiency_summary.pdf")
    plt.close(fig)


def plot_total(
    curves: list[Curve],
    grid: np.ndarray,
    total: np.ndarray,
    constant_scale: float,
    output_dir: Path,
    dpi: int,
) -> None:
    fig, ax = plt.subplots(figsize=(7.8, 4.8))
    colors = {
        "mirror_reflectivity": "#1f77b4",
        "filter_transmission": "#2ca02c",
        "sipm_pde": "#d62728",
        "atmosphere_transmission": "#9467bd",
    }
    for curve in curves:
        ax.plot(grid, curve.evaluate(grid), lw=1.6, color=colors.get(curve.name), label=curve.name)
    ax.plot(grid, total, color="black", lw=2.4, label="total")

    union_wavelengths = sorted(
        {
            float(wl)
            for curve in curves
            if curve.enabled
            for wl in curve.wavelengths
            if grid[0] <= wl <= grid[-1]
        }
    )
    if union_wavelengths:
        union_grid = np.array(union_wavelengths, dtype=float)
        total_at_input_wavelengths = constant_scale * np.prod(
            [curve.evaluate(union_grid) for curve in curves],
            axis=0,
        )
        ax.scatter(
            union_grid,
            total_at_input_wavelengths,
            s=9,
            color="#666666",
            alpha=0.42,
            linewidths=0,
            label="total at input wavelengths",
        )
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("Efficiency")
    ax.set_ylim(-0.04, 1.08)
    ax.set_title("Total wavelength efficiency")
    ax.grid(True, color="#d9d9d9", linewidth=0.7, alpha=0.8)
    ax.legend(frameon=True, framealpha=1.0, edgecolor="black", ncol=2)
    fig.tight_layout()
    fig.savefig(output_dir / "total_efficiency.png", dpi=dpi)
    plt.close(fig)


def write_samples(curves: list[Curve], grid: np.ndarray, total: np.ndarray, output_dir: Path) -> None:
    columns = ["wavelength_nm"] + [curve.name for curve in curves] + ["total_efficiency"]
    with (output_dir / "efficiency_curve_samples.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(columns)
        values = [curve.evaluate(grid) for curve in curves]
        for i, wl in enumerate(grid):
            writer.writerow([f"{wl:.6g}"] + [f"{v[i]:.8g}" for v in values] + [f"{total[i]:.8g}"])


def write_report(curves: list[Curve], grid: np.ndarray, total: np.ndarray, args: argparse.Namespace, output_dir: Path) -> None:
    lines = [
        "Efficiency curve validation report",
        "==================================",
        "",
        f"constant_scale: {args.constant_scale:g}",
        f"wavelength grid: {grid[0]:g} to {grid[-1]:g} nm, step {args.step:g} nm",
        "",
    ]
    for curve in curves:
        if not curve.enabled:
            lines.append(f"[{curve.name}] disabled -> factor 1")
            lines.append("")
            continue
        lines.extend(
            [
                f"[{curve.name}]",
                f"path: {curve.path}",
                f"raw rows: {curve.raw_wavelengths.size}",
                f"unique wavelengths after merge: {curve.wavelengths.size}",
                f"duplicate wavelength groups: {curve.duplicate_groups}",
                f"duplicate rows included in groups: {curve.duplicate_rows}",
                f"rows outside [0, 1]: {curve.out_of_range_rows}",
                f"wavelength range: {curve.wavelengths[0]:g} to {curve.wavelengths[-1]:g} nm",
                f"efficiency range after merge: {np.min(curve.values):.6g} to {np.max(curve.values):.6g}",
                "",
            ]
        )
    lines.extend(
        [
            "[total]",
            f"min total on grid: {np.min(total):.6g}",
            f"max total on grid: {np.max(total):.6g}",
            f"mean total on grid: {np.mean(total):.6g}",
            "",
            "Notes:",
            "- Duplicate wavelengths are merged by arithmetic mean, matching the C++ curve loader.",
            "- A curve evaluates to 0 outside its wavelength table range.",
            "- A disabled component evaluates to 1 and is still shown in the combined plot.",
        ]
    )
    (output_dir / "efficiency_curve_report.txt").write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    grid = np.arange(args.wavelength_min, args.wavelength_max + 0.5 * args.step, args.step)
    curves = [
        read_curve("mirror_reflectivity", args.mirror_csv),
        read_curve("filter_transmission", args.filter_csv),
        read_curve("sipm_pde", args.sipm_csv),
        read_curve("atmosphere_transmission", args.atmosphere_csv),
    ]
    component_values = [curve.evaluate(grid) for curve in curves]
    total = args.constant_scale * np.prod(component_values, axis=0)

    for curve in curves:
        plot_single(curve, grid, output_dir, args.dpi)
    plot_summary(curves, grid, total, args.constant_scale, output_dir, args.dpi)
    plot_total(curves, grid, total, args.constant_scale, output_dir, args.dpi)
    write_samples(curves, grid, total, output_dir)
    write_report(curves, grid, total, args, output_dir)

    print(f"wrote efficiency validation outputs to {output_dir}")
    print(f"total efficiency max={np.max(total):.6g}, mean={np.mean(total):.6g}")
    for curve in curves:
        if curve.enabled and curve.duplicate_groups:
            print(
                f"warning: {curve.name} has {curve.duplicate_groups} duplicate wavelength groups; "
                "merged by arithmetic mean"
            )


if __name__ == "__main__":
    main()
