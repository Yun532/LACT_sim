#!/usr/bin/env python3
"""Prepare compact, provenance-rich parallel-light cases for the coordinate page."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from pathlib import Path


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def number(row: dict[str, str], key: str, default: float = 0.0) -> float:
    value = row.get(key, "")
    return default if value in ("", None) else float(value)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def combined_sha256(paths: list[Path]) -> str:
    digest = hashlib.sha256()
    for path in sorted(paths, key=lambda item: item.as_posix()):
        digest.update(path.name.encode("utf-8"))
        digest.update(b"\0")
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b''):
                digest.update(block)
    return digest.hexdigest()


def uniform_sample(rows: list[dict[str, str]], target: int) -> list[dict[str, str]]:
    if len(rows) <= target:
        return rows
    indices = sorted({round(i * (len(rows) - 1) / (target - 1)) for i in range(target)})
    return [rows[index] for index in indices]


def generic_basis(az_deg: float, el_deg: float) -> dict[str, list[float]]:
    az = math.radians(az_deg)
    el = math.radians(el_deg)
    z = [math.cos(el) * math.cos(az), math.cos(el) * math.sin(az), math.sin(el)]
    x = [-math.sin(az), math.cos(az), 0.0]
    y = [
        z[1] * x[2] - z[2] * x[1],
        z[2] * x[0] - z[0] * x[2],
        z[0] * x[1] - z[1] * x[0],
    ]
    return {"x": x, "y": y, "z": z}


def expand_local(value: list[float], basis: dict[str, list[float]]) -> list[float]:
    return [sum(value[j] * basis[axis][i] for j, axis in enumerate(("x", "y", "z")))
            for i in range(3)]


def vec(row: dict[str, str], prefix: str, suffix: str = '') -> list[float]:
    return [number(row, prefix + axis + suffix) for axis in 'xyz']


def parse_log(path: Path) -> dict[str, float | int | str]:
    text = path.read_text(encoding="utf-8", errors="replace")
    values: dict[str, float | int | str] = {}
    wanted = {
        "total_photons", "blocked_by_obstruction", "blocked_incoming",
        "blocked_reflected", "hit_mirror_before_obstruction",
        "hit_output_before_obstruction", "hit_mirror", "hit_output_plane", "hit_camera",
        "accepted_camera", "unique_hit_pixels", "weighted_spot_rms_mm",
    }
    for line in text.splitlines():
        match = re.fullmatch(r"([a-z_]+)=([^ ]+)", line.strip())
        if not match or match.group(1) not in wanted:
            continue
        key, raw = match.groups()
        try:
            value = float(raw)
            values[key] = int(value) if value.is_integer() else value
        except ValueError:
            values[key] = raw
    return values


def parse_runtime_config(path: Path) -> dict:
    text = path.read_text(encoding="utf-8", errors="replace")
    output: dict = {}
    for key in ("pointing_az_deg", "pointing_el_deg"):
        match = re.search(rf"^\s*{key}\s*:\s*([-+0-9.eE]+)\s*$", text, re.MULTILINE)
        if match:
            output[key] = float(match.group(1))
    direction = re.search(
        r"^\s*beam_direction\s*:\s*\(([^)]+)\)\s*$", text, re.MULTILINE
    )
    if direction:
        output["beam_direction"] = [float(value.strip()) for value in direction.group(1).split(",")]
    mark_only = re.search(r"^\s*mark_only\s*:\s*(true|false)\s*$", text, re.MULTILINE)
    if mark_only:
        output["obstruction_mark_only"] = mark_only.group(1) == "true"
    return output


def compact(value):
    if isinstance(value, float):
        return round(value, 10)
    if isinstance(value, list):
        return [compact(item) for item in value]
    if isinstance(value, dict):
        return {key: compact(item) for key, item in value.items()}
    return value


def load_case(case_id: str, title: str, elevation: float, hit_path: Path,
              camera_path: Path, log_path: Path, config_path: Path,
              sample_rays: int, sample_output: int,
              beam_theta_deg: float = 0.0, beam_phi_deg: float = 0.0,
              source_sky: dict | None = None,
              beam_direction_local: list[float] | None = None,
              include_blocked: bool = False) -> dict:
    azimuth = 0.0
    basis = generic_basis(azimuth, elevation)
    hit_rows = [row for row in read_csv(hit_path) if int(number(row, 'hit_surface')) == 1]
    clear_rows = [row for row in hit_rows if not int(number(row, 'obstruction_blocked'))]
    incoming_blocked_rows = [
        row for row in hit_rows if int(number(row, 'obstruction_blocked_incoming'))
    ]
    reflected_blocked_rows = [
        row for row in hit_rows
        if int(number(row, 'obstruction_blocked_reflected'))
        and not int(number(row, 'obstruction_blocked_incoming'))
    ]
    reflection_rows = [
        row for row in hit_rows
        if not int(number(row, 'obstruction_blocked_incoming'))
    ]
    if include_blocked:
        clear_target = max(1, sample_rays // 2)
        blocked_target = max(1, (sample_rays - clear_target) // 2)
        ray_rows = (
            uniform_sample(clear_rows, clear_target)
            + uniform_sample(incoming_blocked_rows, blocked_target)
            + uniform_sample(reflected_blocked_rows, sample_rays - clear_target - blocked_target)
        )
    else:
        ray_rows = uniform_sample(clear_rows, sample_rays)
    output_rows = uniform_sample(clear_rows, sample_output)
    camera_hit_rows = [row for row in clear_rows if int(number(row, "hit_camera"))]
    rays = []
    for index, row in enumerate(ray_rows):
        input_local = vec(row, 'input_local_', '_m')
        input_dir_local = vec(row, 'input_local_dir_')
        rays.append({
            "photon_index": index,
            "mirror_id": int(number(row, "mirror_id", -1)),
            "input_local": input_local,
            "input_dir_local": input_dir_local,
            "input": expand_local(input_local, basis),
            "input_dir": expand_local(input_dir_local, basis),
            "mirror": vec(row, "mirror_"),
            "surface": vec(row, "surface_"),
            "output_dir": vec(row, "dir_"),
            "u_m": number(row, "u_m"),
            "v_m": number(row, "v_m"),
            "time_ns": number(row, "time_ns"),
            "obstruction_blocked": bool(int(number(row, "obstruction_blocked"))),
            "blocked_incoming": bool(int(number(row, "obstruction_blocked_incoming"))),
            "blocked_reflected": bool(int(number(row, "obstruction_blocked_reflected"))),
        })
    output_points = [{
        "u_m": number(row, "u_m"), "v_m": number(row, "v_m"),
        "time_ns": number(row, "time_ns"),
        "signal_weight": number(row, "weight") * number(row, "relative_efficiency"),
    } for row in output_rows]
    # The 3D ray view is deliberately sampled, but the camera comparison must not be.
    # Keep every physical output-plane u/v and every accepted camera hit from the C++ CSV.
    full_output_uv_m = [
        [number(row, "u_m"), number(row, "v_m")] for row in clear_rows
    ]
    full_camera_hit_uv_m = [
        [number(row, "camera_x_m"), number(row, "camera_y_m")]
        for row in camera_hit_rows
    ]
    # Keep every point needed for the mirror illumination/obstruction diagnostic.
    # Ray segments remain sampled.  For blocked rows the current CSV does not
    # contain the primitive intersection, so preserve the program's theoretical
    # mirror/output endpoint and label it explicitly in the page.
    full_reflection_points_global_m = [
        vec(row, "mirror_") for row in reflection_rows
    ] if include_blocked else []
    full_incoming_blocked_mirror_endpoints_global_m = [
        vec(row, "mirror_") for row in incoming_blocked_rows
    ] if include_blocked else []
    full_reflected_blocked_surface_endpoints_global_m = [
        vec(row, "surface_") for row in reflected_blocked_rows
    ] if include_blocked else []
    camera_signal = [{
        "pixel_id": int(row["pixel_id"]),
        "photon_count": int(row["photon_count"]),
        "pe": number(row, "pe"),
        "time_mean_ns": number(row, "time_mean_ns"),
    } for row in read_csv(camera_path)]
    output_weight = sum(
        max(number(row, "weight") * number(row, "relative_efficiency"), 0.0)
        for row in clear_rows
    )
    output_uv_centroid_m = [0.0, 0.0]
    if output_weight:
        output_uv_centroid_m = [
            sum(
                number(row, "u_m")
                * max(number(row, "weight") * number(row, "relative_efficiency"), 0.0)
                for row in clear_rows
            ) / output_weight,
            sum(
                number(row, "v_m")
                * max(number(row, "weight") * number(row, "relative_efficiency"), 0.0)
                for row in clear_rows
            ) / output_weight,
        ]
    summary = parse_log(log_path)
    runtime_config = parse_runtime_config(log_path)
    runtime_azimuth = float(runtime_config.get("pointing_az_deg", math.nan))
    runtime_elevation = float(runtime_config.get("pointing_el_deg", math.nan))
    if not math.isfinite(runtime_azimuth) or abs(runtime_azimuth - azimuth) > 1e-9:
        raise ValueError(f'{case_id}: runtime pointing azimuth differs from requested value')
    if not math.isfinite(runtime_elevation) or abs(runtime_elevation - elevation) > 1e-9:
        raise ValueError(f'{case_id}: runtime pointing elevation differs from requested value')
    if beam_direction_local is not None:
        logged_direction = runtime_config.get("beam_direction", [])
        if len(logged_direction) != 3 or any(
            abs(actual - expected) > 1.0e-6
            for actual, expected in zip(logged_direction, beam_direction_local)
        ):
            raise ValueError(f'{case_id}: runtime beam direction differs from requested value')
    if include_blocked and runtime_config.get("obstruction_mark_only") is not True:
        raise ValueError(f'{case_id}: runtime did not enable obstruction mark-only diagnostics')
    if int(summary.get('hit_output_plane', -1)) != len(clear_rows):
        raise ValueError(f'{case_id}: physical output count does not match clear CSV rows')
    if include_blocked and int(summary.get('hit_output_before_obstruction', -1)) != len(hit_rows):
        raise ValueError(f'{case_id}: pre-obstruction output count does not match mark-only CSV rows')
    plane = expand_local([0.0, 0.0, -8.0], basis)
    max_output_coordinate_error = 0.0
    max_camera_coordinate_error = 0.0
    for row in clear_rows:
        surface = vec(row, 'surface_')
        offset = [surface[i] - plane[i] for i in range(3)]
        expected_u = sum(offset[i] * basis['x'][i] for i in range(3))
        expected_v = sum(offset[i] * basis['y'][i] for i in range(3))
        max_output_coordinate_error = max(
            max_output_coordinate_error,
            abs(expected_u - number(row, 'u_m')),
            abs(expected_v - number(row, 'v_m')),
        )
        max_camera_coordinate_error = max(
            max_camera_coordinate_error,
            abs(number(row, 'camera_x_m') - number(row, 'u_m')),
            abs(number(row, 'camera_y_m') - number(row, 'v_m')),
        )
    camera_photons = sum(row['photon_count'] for row in camera_signal)
    camera_counts_from_hits: dict[int, int] = {}
    for row in camera_hit_rows:
        pixel_id = int(number(row, "pixel_id", -1))
        camera_counts_from_hits[pixel_id] = camera_counts_from_hits.get(pixel_id, 0) + 1
    camera_counts_from_csv = {
        row["pixel_id"]: row["photon_count"] for row in camera_signal
    }
    if max_output_coordinate_error >= 1.0e-7:
        raise ValueError(f'{case_id}: raw surface does not reproduce raw u/v')
    if max_camera_coordinate_error >= 1.0e-10:
        raise ValueError(f'{case_id}: camera x/y do not equal raw u/v')
    if camera_photons != int(summary.get('hit_camera', -1)):
        raise ValueError(f'{case_id}: pixel photon sum does not match hit_camera')
    if len(camera_signal) != int(summary.get('unique_hit_pixels', -1)):
        raise ValueError(f'{case_id}: non-zero pixel rows do not match unique_hit_pixels')
    if camera_counts_from_hits != camera_counts_from_csv:
        raise ValueError(f'{case_id}: per-pixel hit rows do not match camera CSV photon_count')
    if include_blocked and (
        len(full_reflection_points_global_m) != int(summary.get('hit_mirror', -1))
        or len(full_incoming_blocked_mirror_endpoints_global_m)
        != int(summary.get('blocked_incoming', -1))
        or len(full_reflected_blocked_surface_endpoints_global_m)
        != int(summary.get('blocked_reflected', -1))
    ):
        raise ValueError(f'{case_id}: complete mirror/obstruction point arrays differ from C++ counts')
    mirror_stats: dict[int, dict[str, int]] = {}
    for row in hit_rows:
        mirror_id = int(number(row, "mirror_id", -1))
        stat = mirror_stats.setdefault(mirror_id, {
            "mirror_id": mirror_id,
            "theoretical_hits": 0,
            "clear_reflections": 0,
            "blocked_incoming": 0,
            "blocked_reflected": 0,
        })
        stat["theoretical_hits"] += 1
        if int(number(row, "obstruction_blocked_incoming")):
            stat["blocked_incoming"] += 1
        elif int(number(row, "obstruction_blocked_reflected")):
            stat["blocked_reflected"] += 1
        else:
            stat["clear_reflections"] += 1
    return compact({
        "id": case_id,
        "title": title,
        "pointing": {"az_deg": azimuth, "el_deg": elevation},
        "beam": {
            "theta_deg": beam_theta_deg,
            "phi_deg": beam_phi_deg,
            "direction_local": beam_direction_local,
        },
        "source_sky": source_sky,
        "coordinate_frame": "generic_trace_global",
        "basis": basis,
        "rays": rays,
        "output_points": output_points,
        "full_output_uv_m": full_output_uv_m,
        "full_camera_hit_uv_m": full_camera_hit_uv_m,
        "full_reflection_points_global_m": full_reflection_points_global_m,
        "full_incoming_blocked_mirror_endpoints_global_m": (
            full_incoming_blocked_mirror_endpoints_global_m
        ),
        "full_reflected_blocked_surface_endpoints_global_m": (
            full_reflected_blocked_surface_endpoints_global_m
        ),
        "camera_signal": camera_signal,
        "mirror_stats": [mirror_stats[key] for key in sorted(mirror_stats)],
        "ray_sample_counts": {
            "clear": sum(not ray["obstruction_blocked"] for ray in rays),
            "blocked_incoming": sum(ray["blocked_incoming"] for ray in rays),
            "blocked_reflected": sum(
                ray["blocked_reflected"] and not ray["blocked_incoming"] for ray in rays
            ),
        },
        "camera_summary": {
            "hit_camera": summary.get("hit_camera", 0),
            "accepted_camera": summary.get("accepted_camera", 0),
            "unique_hit_pixels": summary.get("unique_hit_pixels", 0),
            "output_uv_centroid_m": output_uv_centroid_m,
        },
        "summary": summary,
        "runtime_config": runtime_config,
        "validation": {
            "status": "passed",
            "max_surface_to_uv_error_m": max_output_coordinate_error,
            "max_camera_xy_to_uv_error_m": max_camera_coordinate_error,
            "camera_pixel_photon_sum": camera_photons,
            "full_output_uv_rows": len(full_output_uv_m),
            "full_camera_hit_uv_rows": len(full_camera_hit_uv_m),
            "full_reflection_point_rows": len(full_reflection_points_global_m),
            "full_incoming_blocked_endpoint_rows": len(
                full_incoming_blocked_mirror_endpoints_global_m
            ),
            "full_reflected_blocked_endpoint_rows": len(
                full_reflected_blocked_surface_endpoints_global_m
            ),
            "checks": [
                "every hit_surface row is counted by the C++ summary",
                "mark-only rows preserve incoming/reflected obstruction flags when enabled",
                "raw global surface projected on the logged generic basis reproduces raw u/v",
                "camera_x_m/y_m equal u_m/v_m",
                "pixel photon sum equals hit_camera and row count equals unique_hit_pixels",
                "per-pixel hit rows exactly equal the complete camera CSV photon_count",
                "all clear mirror points and all blocked diagnostic endpoints are embedded",
            ],
        },
        "provenance": {
            "run_binary": "build/run_optical_sim",
            "run_config": config_path.as_posix(),
            "run_config_sha256": sha256(config_path),
            "run_overrides": {
                "telescope.pointing_el_deg": elevation,
                "source.beam_direction": beam_direction_local,
                "obstruction.mark_only": include_blocked,
                "output.csv": hit_path.as_posix(),
                "output.pixel_csv": camera_path.as_posix(),
            },
            "hits_csv": hit_path.as_posix(),
            "hits_csv_sha256": sha256(hit_path),
            "camera_csv": camera_path.as_posix(),
            "camera_csv_sha256": sha256(camera_path),
            "log": log_path.as_posix(),
            "log_sha256": sha256(log_path),
            "full_output_rows": len(clear_rows),
            "mark_only_rows": len(hit_rows),
            "ray_sample_rows": len(rays),
            "output_sample_rows": len(output_points),
            "full_output_uv_rows": len(full_output_uv_m),
            "full_camera_hit_uv_rows": len(full_camera_hit_uv_m),
            "full_reflection_point_rows": len(full_reflection_points_global_m),
            "full_incoming_blocked_endpoint_rows": len(
                full_incoming_blocked_mirror_endpoints_global_m
            ),
            "full_reflected_blocked_endpoint_rows": len(
                full_reflected_blocked_surface_endpoints_global_m
            ),
            "ray_sample_method": (
                "uniform clear/incoming-blocked/reflected-blocked strata"
                if include_blocked else "uniform indices over physical hit_surface rows"
            ),
            "coordinate_note": (
                "mirror/surface are untouched generic-global C++ CSV values; "
                "input is the exact buildTelescopeFrame expansion of the untouched "
                "input_local columns so the incoming segment can be drawn in the same frame"
            ),
        },
    })


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-root",
        default="run_logs/coordinate_workbench/server_current_afd630ae/run_logs/coordinate_workbench",
    )
    parser.add_argument("--output", default="docs/assets/data/coordinate-parallel-cases.json")
    parser.add_argument("--source-archive-sha256", default="afd630aedc69825f55eea3961a942ae7904c2b0bb6061f016775f3f10bc49afe")
    parser.add_argument("--results-archive-sha256", default="47493ae480d4e5f00d38c656b39926f956c47da35415f7836c4802e6300cc77f")
    parser.add_argument("--sample-rays", type=int, default=48)
    parser.add_argument("--sample-output", type=int, default=400)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    root = Path(args.input_root)
    baseline = load_case(
        "parallel_baseline_el70", "理想镜片：轴上平行光", 70.0,
        root / "parallel_baseline_hits.csv",
        root / "parallel_baseline_camera.csv",
        root / "parallel_baseline.log",
        Path("configs/examples/coordinate_parallel_camera.cfg"),
        args.sample_rays, args.sample_output,
    )
    cases = []
    for elevation in range(0, 91, 10):
        tag = str(elevation)
        cases.append(load_case(
            f"parallel_deformation_el{tag}", f"支架形变：仰角 {tag}°", float(elevation),
            root / "elevation" / f"hits_el_{tag}.csv",
            root / "elevation" / f"camera_el_{tag}.csv",
            root / "elevation" / f"run_el_{tag}.log",
            Path("configs/examples/coordinate_parallel_deformation_camera.cfg"),
            args.sample_rays, args.sample_output,
        ))
    four_root = Path("run_logs/coordinate_workbench/server_sky_angle/sky-angle")
    four_direction_cases = []
    four_specs = (
        (
            "up", "天区上：光源 az=0° / el=71°", 0.0, 71.0,
            [0.0, -0.0174524064372835, -0.999847695156391], 1.0, 270.0,
        ),
        (
            "down", "天区下：光源 az=0° / el=69°", 0.0, 69.0,
            [0.0, 0.0174524064372835, -0.999847695156391], 1.0, 90.0,
        ),
        (
            "left", "天区左：光源 az=-1° / el=70°", -1.0, 70.0,
            [0.00596907455105753, -0.0000489498331834315, -0.999982183717749],
            0.3420163100786891, 359.53015229443344,
        ),
        (
            "right", "天区右：光源 az=+1° / el=70°", 1.0, 70.0,
            [-0.00596907455105753, -0.0000489498331834315, -0.999982183717749],
            0.3420163100786891, 180.4698477055665,
        ),
    )
    for name, label, source_az, source_el, direction_local, theta, phi in four_specs:
        four_direction_cases.append(load_case(
            f"parallel_sky_{name}", label, 70.0,
            four_root / f"{name}_hits.csv",
            four_root / f"{name}_camera.csv",
            four_root / f"{name}.log",
            Path(f"configs/examples/coordinate_parallel_sky_{name}.cfg"),
            args.sample_rays, args.sample_output, theta, phi,
            source_sky={"az_deg": source_az, "el_deg": source_el},
            beam_direction_local=direction_local,
            include_blocked=True,
        ))
    four_files = [
        four_root / f"{name}{suffix}"
        for name, *_ in four_specs
        for suffix in ("_hits.csv", "_camera.csv", ".log")
    ] + [
        Path(f"configs/examples/coordinate_parallel_sky_{name}.cfg")
        for name, *_ in four_specs
    ]
    output = {
        "status": "ready",
        "source": "current C++ source compiled in an isolated server /tmp directory",
        "source_archive_sha256": args.source_archive_sha256,
        "results_archive_sha256": args.results_archive_sha256,
        "four_direction_results_archive_sha256": combined_sha256(four_files),
        "baseline": baseline,
        "four_direction_cases": four_direction_cases,
        "elevation_cases": cases,
    }
    destination = Path(args.output)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(compact(output), ensure_ascii=False, separators=(",", ":")), encoding="utf-8")
    print(f'Saved {destination} ({destination.stat().st_size} bytes)')


if __name__ == "__main__":
    main()
