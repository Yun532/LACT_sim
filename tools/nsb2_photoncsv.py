#!/usr/bin/env python3
"""Generate expectation or time-resolved LACT PhotonCsv input from nsb2."""

from __future__ import annotations

import argparse
import csv
import importlib.metadata
import importlib.util
import json
import math
from pathlib import Path

import astropy.units as u
import numpy as np
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, SkyOffsetFrame
from astropy.time import Time
from astropy.utils import iers

iers.conf.auto_download = False
iers.conf.auto_max_age = None

CRAB = SkyCoord(ra=83.633083 * u.deg, dec=22.0145 * u.deg, frame="icrs")
LACT = EarthLocation(
    lat=29.3576667 * u.deg,
    lon=100.1387778 * u.deg,
    height=4410 * u.m,
)


def atmosphere():
    from nsb2.atmosphere import SingleScatteringAtmosphere

    def airmass(zenith):
        return 1 / (
            np.cos(zenith)
            + 0.50572 * (96.07995 - np.rad2deg(zenith)) ** (-1.6364)
        )

    def tau_rayleigh(wavelength):
        return 0.00879 * wavelength.to_value(u.micron) ** -4.09 * np.exp(-8 / 8)

    def tau_mie(wavelength):
        return 0.5 * (wavelength.to_value(u.nm) / 380) ** -1 * np.exp(-1.8 / 1.54)

    def tau_absorption(wavelength):
        return np.zeros(wavelength.shape)

    return SingleScatteringAtmosphere(
        airmass, tau_rayleigh, tau_mie, tau_absorption, 0.8
    )


def angular_cells(radius_deg: float, side: int) -> tuple[dict, np.ndarray]:
    edges = np.deg2rad(np.linspace(-radius_deg, radius_deg, side + 1))
    x0, y0 = np.meshgrid(edges[:-1], edges[:-1])
    x1, y1 = np.meshgrid(edges[1:], edges[1:])
    bounds = np.column_stack((x0.ravel(), x1.ravel(), y0.ravel(), y1.ravel()))
    response = {
        "x": bounds[:, :2],
        "y": bounds[:, 2:],
        "values": np.ones((side * side, 2, 2)),
    }
    return response, bounds


def make_bandpass(low_nm: float, high_nm: float):
    from nsb2.core.spectral import Bandpass

    return Bandpass(
        np.linspace(low_nm, high_nm, 5) * u.nm,
        np.ones(5),
    )


def observation(obstime: str) -> SkyOffsetFrame:
    horizon = AltAz(
        obstime=Time(obstime, format="isot", scale="utc"),
        location=LACT,
        pressure=600 * u.hPa,
        temperature=-5 * u.deg_C,
        relative_humidity=0.3,
    )
    return SkyOffsetFrame(origin=CRAB.transform_to(horizon), rotation=0 * u.deg)


def build_sources():
    from nsb2.emitter import airglow, moon, stars, zodiacal

    print("Loading nsb2 emitters...", flush=True)
    glow = airglow.from_eso_skycalc(87 * u.km, 100)
    zodi = zodiacal.from_leinert1998()
    bright = stars.from_gaia_suppl_catalog()
    bright.build_balltree()
    gaia = stars.from_gaia_dr3_catalog()
    gaia.build_balltree()
    gaia_scatter = gaia.to_map(2**4)
    faint = stars.from_gaia_dr3_map()
    moon_model = moon.from_noll2013()
    return glow, zodi, bright, gaia, gaia_scatter, faint, moon_model


def predict_cube(args, response, obs):
    from nsb2.core.instrument import EffectiveApertureInstrument
    from nsb2.core.lightpath import DirectPath, ScatteredPath
    from nsb2.core.pipeline import Pipeline
    from nsb2.core.solver import LUTDirectSolver, LUTScatteredSolver

    glow, zodi, bright, gaia, gaia_scatter, faint, moon_model = build_sources()
    edges = np.linspace(
        args.wavelength_min_nm,
        args.wavelength_max_nm,
        args.bands + 1,
    )
    rates = []
    component_names = None
    atmo = atmosphere()
    for index, (low, high) in enumerate(zip(edges[:-1], edges[1:]), 1):
        print(f"nsb2 band {index}/{args.bands}: {low:.0f}-{high:.0f} nm", flush=True)
        instrument = EffectiveApertureInstrument(
            response, make_bandpass(low, high)
        )
        diffuse = Pipeline(
            instrument,
            atmo,
            [glow, zodi, moon_model],
            [DirectPath(), ScatteredPath()],
        )
        direct_stars = Pipeline(
            instrument,
            atmo,
            [bright, gaia, faint],
            DirectPath(LUTDirectSolver()),
        )
        scattered_stars = Pipeline(
            instrument,
            atmo,
            [bright, gaia_scatter, faint],
            ScatteredPath(LUTScatteredSolver()),
        )
        model = direct_stars + scattered_stars + diffuse
        model.compile(
            extinction_z_bins=30,
            scattering_z_bins=10,
            scattering_theta_bins=10,
        )
        prediction = model.predict(obs)
        names = [f"{item.source_name}|{item.path_name}" for item in prediction]
        if component_names is None:
            component_names = names
        elif names != component_names:
            raise RuntimeError("nsb2 components changed between wavelength bands")
        rates.append(
            np.stack([item.rates.to_value(u.Hz)[:, 1] for item in prediction])
        )
    return edges, component_names, np.stack(rates)


def component_mask(names, include, exclude, exclude_stars):
    lowered = [name.lower() for name in names]
    mask = np.ones(len(names), dtype=bool)
    if include:
        tokens = [token.lower() for token in include]
        mask = np.array([any(token in name for token in tokens) for name in lowered])
    for token in (item.lower() for item in exclude):
        mask &= np.array([token not in name for name in lowered])
    if exclude_stars:
        mask &= np.array([
            not any(token in name for token in ("star", "gaia"))
            for name in lowered
        ])
    if not mask.any():
        raise ValueError("component selection removed every nsb2 component")
    return mask


def sample_geometry(rng, count, bounds, pupil_radius_m, wave_low, wave_high):
    lon0, lon1, lat0, lat1 = bounds
    lon = rng.uniform(lon0, lon1, count)
    lat = np.arcsin(rng.uniform(np.sin(lat0), np.sin(lat1), count))
    cos_lat = np.cos(lat)
    direction = np.column_stack((
        cos_lat * np.sin(lon),
        np.sin(lat),
        -cos_lat * np.cos(lon),
    ))
    radius = pupil_radius_m * np.sqrt(rng.random(count))
    angle = rng.uniform(0.0, 2.0 * np.pi, count)
    wavelength = rng.uniform(wave_low, wave_high, count)
    return radius * np.cos(angle), radius * np.sin(angle), direction, wavelength


def write_photons(
    path,
    args,
    band_edges,
    angular_bounds,
    selected_rates_hz,
    generation_start_ns,
    generation_end_ns,
):
    rng = np.random.default_rng(args.seed)
    area_m2 = math.pi * args.pupil_radius_m**2
    duration_ns = (
        args.exposure_ns
        if args.mode == "expectation"
        else generation_end_ns - generation_start_ns
    )
    expected = selected_rates_hz * area_m2 * duration_ns * 1e-9
    flat_expected = np.clip(expected.ravel(), 0.0, None)
    total_expected = float(flat_expected.sum())
    if not np.isfinite(total_expected) or total_expected <= 0:
        raise RuntimeError("nsb2 predicted no photons for the selected components")

    if args.mode == "expectation":
        counts = rng.multinomial(args.rays, flat_expected / total_expected)
        multiplicity = total_expected / args.rays
        realized_rows = int(args.rays)
    else:
        counts = rng.poisson(flat_expected)
        multiplicity = 1.0
        realized_rows = int(counts.sum())
        if args.max_timed_photons > 0 and realized_rows > args.max_timed_photons:
            raise RuntimeError(
                f"timed PhotonCsv needs {realized_rows} rows; increase "
                "--max-timed-photons or shorten the generated window"
            )

    header = [
        "x_m", "y_m", "z_m", "dir_x", "dir_y", "dir_z",
    ]
    if args.mode == "timed":
        header.append("time_ns")
    header += [
        "wavelength_nm", "weight", "multiplicity", "event_id",
        "telescope_id", "origin",
    ]
    n_cells = len(angular_bounds)
    written = 0
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(header)
        for flat_index, total in enumerate(counts):
            if total == 0:
                continue
            band_index, cell_index = divmod(flat_index, n_cells)
            remaining = int(total)
            while remaining:
                count = min(remaining, args.batch_rows)
                x, y, direction, wavelength = sample_geometry(
                    rng,
                    count,
                    angular_bounds[cell_index],
                    args.pupil_radius_m,
                    band_edges[band_index],
                    band_edges[band_index + 1],
                )
                times = (
                    rng.uniform(generation_start_ns, generation_end_ns, count)
                    if args.mode == "timed" else None
                )
                for row in range(count):
                    values = [x[row], y[row], 1.0, *direction[row]]
                    if times is not None:
                        values.append(times[row])
                    values += [
                        wavelength[row], 1.0, multiplicity,
                        args.event_id, args.telescope_id, "nsb",
                    ]
                    writer.writerow(values)
                written += count
                remaining -= count
    if written != realized_rows:
        raise RuntimeError("internal PhotonCsv row-count mismatch")
    return expected, total_expected, realized_rows, area_m2


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--mode", choices=("expectation", "timed"), required=True)
    parser.add_argument("--rate-cube", type=Path)
    parser.add_argument("--obstime", default="2026-12-05T19:20:00")
    parser.add_argument("--radius-deg", type=float, default=4.5)
    parser.add_argument("--grid-side", type=int, default=41)
    parser.add_argument("--wavelength-min-nm", type=float, default=300.0)
    parser.add_argument("--wavelength-max-nm", type=float, default=900.0)
    parser.add_argument("--bands", type=int, default=12)
    parser.add_argument("--pupil-radius-m", type=float, default=4.0)
    parser.add_argument("--exposure-ns", type=float, default=1e9)
    parser.add_argument("--integration-start-ns", type=float, default=0.0)
    parser.add_argument("--integration-end-ns", type=float, default=2000.0)
    parser.add_argument("--guard-pre-ns", type=float, default=250.0)
    parser.add_argument("--guard-post-ns", type=float, default=250.0)
    parser.add_argument("--rays", type=int, default=1_000_000)
    parser.add_argument("--max-timed-photons", type=int, default=50_000_000)
    parser.add_argument("--batch-rows", type=int, default=100_000)
    parser.add_argument("--include-component", action="append", default=[])
    parser.add_argument("--exclude-component", action="append", default=[])
    parser.add_argument("--exclude-stars", action="store_true")
    parser.add_argument("--event-id", type=int, default=0)
    parser.add_argument("--telescope-id", type=int, default=0)
    parser.add_argument("--seed", type=int, default=20260825)
    args = parser.parse_args()
    if args.mode == "expectation" and args.exposure_ns <= 0:
        parser.error("--exposure-ns must be > 0")
    if args.mode == "timed" and not (
        args.integration_end_ns > args.integration_start_ns
        and args.guard_pre_ns >= 0
        and args.guard_post_ns >= 0
    ):
        parser.error("invalid integration window or guard")
    if args.rays <= 0 or args.batch_rows <= 0:
        parser.error("--rays and --batch-rows must be > 0")
    return args


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    obs = observation(args.obstime)
    if args.rate_cube:
        cached = np.load(args.rate_cube, allow_pickle=False)
        band_edges = cached["band_edges_nm"]
        bounds = cached["angular_bounds_rad"]
        component_names = cached["component_names"].astype(str).tolist()
        component_rates = cached["component_rates_hz"]
    else:
        response, bounds = angular_cells(args.radius_deg, args.grid_side)
        band_edges, component_names, component_rates = predict_cube(
            args, response, obs
        )
    selected = component_mask(
        component_names,
        args.include_component,
        args.exclude_component,
        args.exclude_stars,
    )
    selected_rates = component_rates[:, selected, :].sum(axis=1)
    generated_start = (
        args.integration_start_ns - args.guard_pre_ns
        if args.mode == "timed" else 0.0
    )
    generated_end = (
        args.integration_end_ns + args.guard_post_ns
        if args.mode == "timed" else args.exposure_ns
    )

    photon_csv = args.output_dir / f"nsb2_crab_{args.mode}_photons.csv"
    expected, total_expected, rows, area_m2 = write_photons(
        photon_csv,
        args,
        band_edges,
        bounds,
        selected_rates,
        generated_start,
        generated_end,
    )
    rate_cube = args.output_dir / "nsb2_rate_cube.npz"
    np.savez_compressed(
        rate_cube,
        band_edges_nm=band_edges,
        angular_bounds_rad=bounds,
        component_names=np.asarray(component_names),
        component_rates_hz=component_rates,
        selected_component_mask=selected,
        selected_rates_hz=selected_rates,
        expected_input_photons=expected,
    )
    origin_altaz = obs.origin
    manifest = {
        "model": "nsb2 full sky model",
        "nsb2_version": (
            importlib.metadata.version("nsb2")
            if importlib.util.find_spec("nsb2") is not None
            else "rate-cube-only"
        ),
        "mode": args.mode,
        "obstime_utc": args.obstime,
        "site": {"lat_deg": 29.3576667, "lon_deg": 100.1387778, "height_m": 4410},
        "pointing": {
            "target": "Crab",
            "ra_deg": CRAB.ra.deg,
            "dec_deg": CRAB.dec.deg,
            "alt_deg": origin_altaz.alt.deg,
            "az_deg": origin_altaz.az.deg,
        },
        "field_radius_deg": args.radius_deg,
        "wavelength_edges_nm": band_edges.tolist(),
        "launch_area_m2": area_m2,
        "integration_window_ns": (
            [args.integration_start_ns, args.integration_end_ns]
            if args.mode == "timed" else [0.0, args.exposure_ns]
        ),
        "generated_time_window_ns": [generated_start, generated_end],
        "selected_components": [
            name for name, keep in zip(component_names, selected) if keep
        ],
        "sample_rows": rows,
        "expected_input_photons": total_expected,
        "photon_csv": str(photon_csv),
        "rate_cube": str(rate_cube),
        "notes": [
            "nsb2 already applies atmosphere; use atmosphere=ideal in LACT_sim.",
            "Timed mode uses independent Poisson aperture photons with multiplicity=1.",
        ],
    }
    manifest_path = args.output_dir / f"nsb2_crab_{args.mode}_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2), flush=True)


if __name__ == "__main__":
    main()
