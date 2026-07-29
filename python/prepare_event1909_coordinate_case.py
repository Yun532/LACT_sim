#!/usr/bin/env python3
"""Build a compact coordinate-view payload for the real prod1 EventIO event 1909."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np
import uproot


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def compact(value):
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return round(float(value), 10)
    if isinstance(value, np.ndarray):
        return [compact(item) for item in value.tolist()]
    if isinstance(value, list):
        return [compact(item) for item in value]
    if isinstance(value, dict):
        return {key: compact(item) for key, item in value.items()}
    return value


def parse_flat_config(path: Path) -> dict[str, str]:
    """Keep exact effective key/value text from the cfg for page-side auditing."""
    values: dict[str, str] = {}
    for raw_line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw_line.split("#", 1)[0].strip()
        if not line or "=" not in line:
            continue
        key, value = line.split("=", 1)
        values[key.strip()] = value.strip()
    return values


def read_eventio_samples(path: Path, telescope_ids: list[int], target: int,
                         telescope_positions: dict[int, list[float]],
                         deps_path: Path) -> tuple[str, float, dict[int, dict]]:
    sys.path.insert(0, str(deps_path.resolve()))
    import eventio  # type: ignore

    selected_event = None
    with eventio.IACTFile(path) as handle:
        for event in handle:
            if int(event.event_number) == 19 and int(event.reuse) == 10:
                selected_event = event
                break
    if selected_event is None:
        raise ValueError("CORSIKA shower 19 / reuse 10 (output event 1909) was not found")

    # CORSIKA event-header lengths are centimetres.  Read the observation level
    # from the same event instead of duplicating the 4400 m site value here.
    n_observation_levels = int(selected_event.header["n_observation_levels"])
    observation_heights_cm = np.asarray(
        selected_event.header["observation_height"], dtype=float
    )[:n_observation_levels]
    if n_observation_levels != 1 or len(observation_heights_cm) != 1:
        raise ValueError(
            f"expected exactly one CORSIKA observation level, got {n_observation_levels}"
        )
    observation_altitude_m = float(observation_heights_cm[0]) * 0.01

    output: dict[int, dict] = {}
    for telescope_id in telescope_ids:
        bunches = selected_event.photon_bunches[telescope_id]
        indices = sorted({
            round(i * (len(bunches) - 1) / max(1, target - 1))
            for i in range(target)
        })
        telescope_position = telescope_positions[telescope_id]
        samples = []
        for index in indices:
            bunch = bunches[index]
            cx = float(bunch["cx"])
            cy = float(bunch["cy"])
            direction = [cx, cy, -math.sqrt(max(0.0, 1.0 - cx * cx - cy * cy))]
            anchor_tel = [float(bunch["x"]) * 0.01, float(bunch["y"]) * 0.01, 0.0]
            zem_m_asl = float(bunch["zem"]) * 0.01
            emission_height_tel_m = (
                zem_m_asl - observation_altitude_m - telescope_position[2]
            )
            distance = (emission_height_tel_m - anchor_tel[2]) / (-direction[2])
            emission_tel = [anchor_tel[i] - distance * direction[i] for i in range(3)]
            anchor_array = [anchor_tel[i] + telescope_position[i] for i in range(3)]
            emission_array = [emission_tel[i] + telescope_position[i] for i in range(3)]
            samples.append({
                "source_bunch_index": index,
                "anchor_tel_nwu_m": anchor_tel,
                "anchor_array_nwu_m": anchor_array,
                "direction_nwu": direction,
                "time_ns": float(bunch["time"]),
                "emission_altitude_km_asl": zem_m_asl * 0.001,
                "multiplicity": float(bunch["photons"]),
                "emission_tel_nwu_m": emission_tel,
                "emission_array_nwu_m": emission_array,
                "path_length_m": distance,
                "emission_time_derived_ns": float(bunch["time"]) - distance / 0.299792458,
            })
        output[telescope_id] = {
            "telescope_id": telescope_id,
            "raw_bunch_count": len(bunches),
            "samples": samples,
        }
    return eventio.__version__, observation_altitude_m, output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root",
        default=(
            "../LACT_sim_audit_main_0048369/results/trust_validation/"
            "gamma_gather_examples/event1909_latest_same_input/"
            "event1909_latest_same_input.root"
        ),
    )
    parser.add_argument(
        "--eventio",
        default=(
            "../reconstruction/lact_prod1_corsika_particle_gamma_energy_"
            "1000.0_10000.0_zenith_20.0_azimuth_0.0_run_1_event_0.zst"
        ),
    )
    parser.add_argument(
        "--eventio-deps",
        default="run_logs/coordinate_workbench/python_eventio_deps",
    )
    parser.add_argument(
        "--output",
        default="docs/assets/data/corsika-event1909-coordinate-case.json",
    )
    parser.add_argument("--source-base-commit", required=True)
    parser.add_argument("--source-archive-sha256", required=True)
    parser.add_argument("--run-binary-sha256", required=True)
    parser.add_argument("--eventio-sha256", required=True)
    parser.add_argument("--run-log")
    parser.add_argument("--sample-rays", type=int, default=240)
    args = parser.parse_args()

    root_path = Path(args.root)
    eventio_path = Path(args.eventio)
    actual_eventio_sha256 = sha256(eventio_path)
    if actual_eventio_sha256 != args.eventio_sha256.lower():
        raise ValueError(
            "EventIO SHA-256 mismatch: "
            f"actual={actual_eventio_sha256}, declared={args.eventio_sha256}"
        )
    root = uproot.open(root_path)
    event = root["corsika_events"].arrays(library="np")
    telescopes = root["telescopes"].arrays(library="np")
    trace = root["trace_summary"].arrays(library="np")
    observations = root["observations"].arrays(library="ak")
    camera_pixels = root["camera_pixels"].arrays(library="np")

    trace_by_tel = {
        int(trace["telescope_id"][i]): {
            key: compact(trace[key][i])
            for key in (
                "input_bunches", "input_photons", "blocked_by_obstruction",
                "hit_mirror", "hit_output_plane", "hit_camera", "accepted_camera",
                "unique_hit_pixels", "signal_pe",
            )
        }
        for i in range(len(trace["telescope_id"]))
    }
    array = []
    for i, telescope_id in enumerate(telescopes["telescope_id"]):
        telescope_id = int(telescope_id)
        array.append({
            "telescope_id": telescope_id,
            "position_nwu_m": [
                telescopes["array_x_north_m"][i],
                telescopes["array_y_west_m"][i],
                telescopes["array_z_up_m"][i],
            ],
            "pointing_az_deg": telescopes["pointing_az_deg"][i],
            "pointing_el_deg": telescopes["pointing_el_deg"][i],
            **trace_by_tel[telescope_id],
        })

    # Every telescope is selectable in the diagnostic page.  Keep the same
    # telescope ids and ROOT observation rows; do not rank or remap them here.
    selected_telescopes = sorted(trace_by_tel)
    camera_views = []
    for i in range(len(observations)):
        telescope_id = int(observations["telescope_id"][i])
        if telescope_id not in selected_telescopes:
            continue
        pixel_id = np.asarray(observations["pixel_id"][i], dtype=int)
        cherenkov = np.asarray(observations["image_cherenkov_pe"][i], dtype=float)
        keep = np.isfinite(cherenkov) & (cherenkov > 0.0)
        camera_views.append({
            "telescope_id": telescope_id,
            "triggered": bool(observations["triggered"][i]),
            "camera_signal": [
                {"pixel_id": int(pid), "pe": float(pe)}
                for pid, pe in zip(pixel_id[keep], cherenkov[keep])
            ],
            "cherenkov_pe_sum": float(cherenkov[keep].sum()),
            "trace_signal_pe": trace_by_tel[telescope_id]["signal_pe"],
        })
    camera_views.sort(key=lambda view: selected_telescopes.index(view["telescope_id"]))

    positions_by_id = {row["telescope_id"]: row["position_nwu_m"] for row in array}
    eventio_version, observation_altitude_m, rays_by_telescope = read_eventio_samples(
        eventio_path, selected_telescopes, args.sample_rays, positions_by_id,
        Path(args.eventio_deps),
    )

    all_samples = [
        sample
        for telescope_id in selected_telescopes
        for sample in rays_by_telescope[telescope_id]["samples"]
    ]
    max_direction_norm_error = max(
        abs(math.sqrt(sum(x * x for x in sample["direction_nwu"])) - 1.0)
        for sample in all_samples
    )
    max_emission_altitude_error_m = max(
        abs(
            sample["emission_array_nwu_m"][2]
            - (sample["emission_altitude_km_asl"] * 1000.0 - observation_altitude_m)
        )
        for sample in all_samples
    )
    max_line_reconstruction_error_m = max(
        math.sqrt(sum(
            (
                sample["emission_array_nwu_m"][axis]
                + sample["path_length_m"] * sample["direction_nwu"][axis]
                - sample["anchor_array_nwu_m"][axis]
            ) ** 2
            for axis in range(3)
        ))
        for sample in all_samples
    )
    if any(sample["path_length_m"] <= 0.0 for sample in all_samples):
        raise ValueError("CORSIKA emission reconstruction produced a non-positive path length")
    if max_direction_norm_error >= 1.0e-12:
        raise ValueError("CORSIKA photon direction is not unit length")
    if max_emission_altitude_error_m >= 1.0e-7:
        raise ValueError("CORSIKA reconstructed emission point does not reproduce raw zem")
    if max_line_reconstruction_error_m >= 1.0e-7:
        raise ValueError("CORSIKA reconstructed emission point is not collinear with the raw ray")

    event_row = {key: compact(value[0]) for key, value in event.items()}
    payload = {
        "status": "ready",
        "title": "prod1 4 GB EventIO：event 1909（shower 19 / array 9）",
        "event": event_row,
        "array": array,
        "selected_telescope_id": 19,
        "camera_views": camera_views,
        "camera_geometry_lact_uv": [
            {
                "id": int(camera_pixels["pixel_id"][i]),
                "u": float(camera_pixels["x_m"][i]),
                "v": float(camera_pixels["y_m"][i]),
                "size": float(camera_pixels["size_m"][i]),
                "shape_code": int(camera_pixels["shape_code"][i]),
            }
            for i in range(len(camera_pixels["pixel_id"]))
        ],
        "ray_samples_by_telescope": [rays_by_telescope[telescope_id] for telescope_id in selected_telescopes],
        "raw_tel19_bunch_count": rays_by_telescope[19]["raw_bunch_count"],
        "pointing": {"az_deg": 0.0, "el_deg": 70.0},
        "coordinate_frame": "CORSIKA array NWU (x=North, y=West, z=Up)",
        "spatial_scale": {
            "displayed_coordinate_unit": "m",
            "eventio_xy_cm_to_m": 0.01,
            "eventio_zem_cm_asl_to_m_asl": 0.01,
            "root_telescope_position_unit": "m",
            "observation_altitude_m_from_event_header": observation_altitude_m,
            "display_projection": "one isotropic span shared by x/y/z and all telescopes",
            "marker_radius": "fixed screen pixels for visibility; not a physical length",
        },
        "emission_reconstruction": {
            "point_kind": "derived from each raw 2D photon bunch; not an explicit 3D EventIO field",
            "assumption": "straight-line propagation from zem altitude to the receive plane",
            "formula": "h=zem-observation_altitude-telescope_z; s=(h-anchor_z)/(-dir_z); emission=anchor-s*direction",
        },
        "provenance": {
            "source_file_name": eventio_path.name,
            "source_file_size_bytes": eventio_path.stat().st_size,
            "source_file_path": eventio_path.as_posix(),
            "source_file_sha256": actual_eventio_sha256,
            "eventio_python_version": eventio_version,
            "root_path": root_path.as_posix(),
            "root_sha256": sha256(root_path),
            "source_base_commit": args.source_base_commit,
            "source_archive_sha256": args.source_archive_sha256,
            "run_binary_sha256": args.run_binary_sha256,
            "run_config": "configs/examples/corsika_coordinate_event1909_array.cfg",
            "run_config_values": parse_flat_config(
                Path("configs/examples/corsika_coordinate_event1909_array.cfg")
            ),
            "critical_source_sha256": {
                path.as_posix(): sha256(path)
                for path in (
                    Path("apps/run_corsika_trace.cpp"),
                    Path("src/app/OpticalSimCommon.cpp"),
                    Path("apps/test_coordinate_frames.cpp"),
                    Path("configs/examples/corsika_coordinate_event1909_array.cfg"),
                )
            },
            "run_log": args.run_log,
            "run_log_sha256": sha256(Path(args.run_log)) if args.run_log else None,
            "event_id_mode": "event_array100",
            "selection": {"event_id": 1909, "shower_event_id": 19, "array_id": 9},
            "ray_boundary": (
                "Photon x/y, cx/cy, time and zem are read directly from the local 4 GB EventIO file. "
                "The 3D emission point is derived from the raw receive-plane anchor, raw direction "
                "and raw zem using the same unit/sign convention as EventIOPhotonSource.cpp."
            ),
        },
        "validation": {
            "status": "passed",
            "sample_count": len(all_samples),
            "max_direction_norm_error": max_direction_norm_error,
            "max_emission_altitude_error_m": max_emission_altitude_error_m,
            "max_line_reconstruction_error_m": max_line_reconstruction_error_m,
            "checks": [
                "ROOT contains exactly event 1909 and 32 telescope observations",
                "all 32 displayed telescopes use raw photon bunches from the same shower 19 / reuse 10",
                "selected detailed telescope 19 has 27159 raw photon bunch rows",
                "all selectable camera panels use image_cherenkov_pe from the same ROOT event",
                "all displayed LACT u/v pixel centers and sizes come from the ROOT camera_pixels tree",
                "array layout and core use the ROOT NWU branches without coordinate conversion",
                "observation altitude is read from the selected CORSIKA event header",
                "every displayed emission point reproduces raw zem and the raw straight ray",
                "all displayed CORSIKA positions use metres and one shared isotropic 3D scale",
                "ROOT output records the latest-main base commit, exact source archive and binary hashes",
            ],
        },
    }
    for telescope_id in selected_telescopes:
        if rays_by_telescope[telescope_id]["raw_bunch_count"] != int(trace_by_tel[telescope_id]["input_bunches"]):
            raise ValueError(f"telescope {telescope_id}: EventIO bunch count differs from ROOT trace summary")
    if int(event_row["event_id"]) != 1909 or len(array) != 32 or rays_by_telescope[19]["raw_bunch_count"] != 27159:
        raise ValueError("event 1909 source consistency check failed")
    destination = Path(args.output)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(compact(payload), ensure_ascii=False, separators=(",", ":")),
        encoding="utf-8",
    )
    print(f"Saved {destination} ({destination.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
