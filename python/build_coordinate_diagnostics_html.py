#!/usr/bin/env python3
"""Build the self-contained Chinese LACT coordinate diagnostic workbench."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(line for line in handle if not line.lstrip().startswith("#")))


def number(row: dict[str, str], key: str, default: float = 0.0) -> float:
    value = row.get(key, "")
    return default if value in (None, "") else float(value)


def vec(row: dict[str, str], prefix: str) -> list[float]:
    return [number(row, prefix + axis) for axis in ("x", "y", "z")]


def vec_m(row: dict[str, str], prefix: str) -> list[float]:
    return [number(row, prefix + axis + "_m") for axis in ("x", "y", "z")]


def unit(v: list[float]) -> list[float]:
    length = math.sqrt(sum(x * x for x in v))
    return [x / length for x in v]


def cross(a: list[float], b: list[float]) -> list[float]:
    return [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]


def local_basis(normal: list[float]) -> tuple[list[float], list[float]]:
    n = unit(normal)
    ref = [0.0, 0.0, 1.0] if abs(n[2]) < 0.9 else [0.0, 1.0, 0.0]
    u = unit(cross(ref, n))
    return u, cross(n, u)


def add(a: list[float], b: list[float], scale: float = 1.0) -> list[float]:
    return [a[i] + scale * b[i] for i in range(3)]


def facet_polygon(center: list[float], normal: list[float], size: float) -> list[list[float]]:
    u, v = local_basis(normal)
    radius = size / 2.0
    xy = [
        (radius, radius / math.sqrt(3.0)),
        (0.0, 2.0 * radius / math.sqrt(3.0)),
        (-radius, radius / math.sqrt(3.0)),
        (-radius, -radius / math.sqrt(3.0)),
        (0.0, -2.0 * radius / math.sqrt(3.0)),
        (radius, -radius / math.sqrt(3.0)),
    ]
    return [add(add(center, u, x), v, y) for x, y in xy]


def box_edges(center: list[float], half: list[float]) -> list[list[list[float]]]:
    x, y, z = center
    hx, hy, hz = half
    c = [
        [x-hx, y-hy, z-hz], [x+hx, y-hy, z-hz],
        [x+hx, y+hy, z-hz], [x-hx, y+hy, z-hz],
        [x-hx, y-hy, z+hz], [x+hx, y-hy, z+hz],
        [x+hx, y+hy, z+hz], [x-hx, y+hy, z+hz],
    ]
    return [
        [c[0], c[1], c[2], c[3], c[0]],
        [c[4], c[5], c[6], c[7], c[4]],
        [c[0], c[4]], [c[1], c[5]], [c[2], c[6]], [c[3], c[7]],
    ]


def prism_edges(center: list[float], radius: float, height: float,
                rotation: float, sides: int) -> list[list[list[float]]]:
    angles = [rotation + i * 2.0 * math.pi / sides for i in range(sides)]
    bottom = [[center[0] + radius * math.cos(a), center[1] + radius * math.sin(a),
               center[2] - height / 2.0] for a in angles]
    top = [[p[0], p[1], center[2] + height / 2.0] for p in bottom]
    return [bottom + [bottom[0]], top + [top[0]]] + [[a, b] for a, b in zip(bottom, top)]


def cylinder_edges(p0: list[float], p1: list[float], radius: float,
                   sides: int = 8) -> list[list[list[float]]]:
    if radius <= 0.0:
        return [[p0, p1]]
    axis = unit([p1[i] - p0[i] for i in range(3)])
    ref = [0.0, 0.0, 1.0] if abs(axis[2]) < 0.9 else [0.0, 1.0, 0.0]
    u = unit(cross(axis, ref))
    v = cross(axis, u)
    ring0, ring1 = [], []
    for index in range(sides):
        angle = 2.0 * math.pi * index / sides
        offset = [radius * (math.cos(angle) * u[i] + math.sin(angle) * v[i])
                  for i in range(3)]
        ring0.append(add(p0, offset))
        ring1.append(add(p1, offset))
    return [ring0 + [ring0[0]], ring1 + [ring1[0]]] + [
        [ring0[i], ring1[i]] for i in range(sides)
    ]


def read_geometry(ideal_path: Path, deformation_path: Path,
                  primitive_path: Path) -> tuple[list[dict], list[dict], list[dict]]:
    ideal = []
    by_id: dict[int, dict] = {}
    for row in read_csv(ideal_path):
        item = {
            "id": int(row["id"]),
            "center": vec(row, "center_"),
            "normal": vec(row, "normal_"),
            "size": number(row, "size1", 0.8),
        }
        item["polygon"] = facet_polygon(item["center"], item["normal"], item["size"])
        ideal.append(item)
        by_id[item["id"]] = item

    series: dict[float, list[dict]] = defaultdict(list)
    for row in read_csv(deformation_path):
        elevation = number(row, "elevation_deg")
        base = by_id[int(row["id"])]
        center = vec(row, "center_")
        normal = vec(row, "normal_")
        delta = [center[i] - base["center"][i] for i in range(3)]
        dot_value = max(-1.0, min(1.0, sum(
            unit(base["normal"])[i] * unit(normal)[i] for i in range(3)
        )))
        series[elevation].append({
            "id": base["id"], "center": center, "normal": normal,
            "delta": delta,
            "delta_mm": math.sqrt(sum(x * x for x in delta)) * 1000.0,
            "normal_delta_mdeg": math.degrees(math.acos(dot_value)) * 1000.0,
            "polygon": facet_polygon(center, normal, base["size"]),
        })

    deformation = []
    for elevation in sorted(series):
        facets = sorted(series[elevation], key=lambda item: item["id"])
        components = [[f["delta"][i] * 1000.0 for f in facets] for i in range(3)]
        values = [f["delta_mm"] for f in facets]
        deformation.append({
            "elevation_deg": elevation,
            "zenith_deg": 90.0 - elevation,
            "facets": facets,
            "stats": {
                "max_mm": max(values),
                "rms_mm": math.sqrt(sum(x * x for x in values) / len(values)),
                "mean_dx_mm": sum(components[0]) / len(facets),
                "mean_dy_mm": sum(components[1]) / len(facets),
                "mean_dz_mm": sum(components[2]) / len(facets),
                "max_normal_mdeg": max(f["normal_delta_mdeg"] for f in facets),
            },
        })

    primitive_lines = []
    for row in read_csv(primitive_path):
        role = row.get("role", "")
        name = row.get("name", "")
        kind = row.get("type", "").lower()
        if kind == "cylinder":
            p0 = [number(row, f"{axis}0_m") for axis in "xyz"]
            p1 = [number(row, f"{axis}1_m") for axis in "xyz"]
            for edge in cylinder_edges(p0, p1, number(row, "radius_m")):
                primitive_lines.append({"points": edge, "role": role, "name": name})
        elif kind in {"box", "aabb"}:
            center = [number(row, f"center_{axis}_m", number(row, f"{axis}0_m")) for axis in "xyz"]
            half = [number(row, f"half_{axis}_m") for axis in "xyz"]
            for edge in box_edges(center, half):
                primitive_lines.append({"points": edge, "role": role, "name": name})
            if number(row, "hole_radius_m") > 0.0:
                for edge in prism_edges(
                    center, number(row, "hole_radius_m"), half[2] * 2.0,
                    number(row, "hole_rotation_rad"), int(number(row, "hole_sides", 6)),
                ):
                    primitive_lines.append({
                        "points": edge, "role": "camera_adapter_hole",
                        "name": name + " opening",
                    })
        elif kind == "polygon_prism":
            center = [number(row, f"center_{axis}_m", number(row, f"{axis}0_m")) for axis in "xyz"]
            for edge in prism_edges(center, number(row, "radius_m"), number(row, "height_m"),
                                    number(row, "rotation_rad"), int(number(row, "sides", 6))):
                primitive_lines.append({"points": edge, "role": role, "name": name})
    return ideal, deformation, primitive_lines


def read_traces(path: Path) -> list[dict]:
    traces = []
    for row in read_csv(path):
        traces.append({
            "photon_index": int(row["photon_index"]),
            "mirror_id": int(row["mirror_id"]),
            "source_bunch_index": int(row["source_bunch_index"]),
            "input": vec_m(row, "input_"),
            "input_dir": vec(row, "input_dir_"),
            "mirror": vec_m(row, "mirror_"),
            "surface": vec_m(row, "surface_"),
            "u_m": number(row, "u_m"), "v_m": number(row, "v_m"),
            "raw_input": vec_m(row, "raw_input_"),
            "raw_input_dir": vec(row, "raw_input_dir_"),
            "output_dir": vec(row, "dir_"),
            "time_ns": number(row, "time_ns"),
            "wavelength_nm": number(row, "wavelength_nm"),
        })
    return traces


def read_output_plane(path: Path) -> list[dict]:
    return [{
        "photon_index": int(row["photon_index"]),
        "source_bunch_index": int(row["source_bunch_index"]),
        "u_m": number(row, "u_m"), "v_m": number(row, "v_m"),
        "signal_weight": number(row, "signal_weight"),
        "time_ns": number(row, "time_ns"),
    } for row in read_csv(path)]


def read_shower(path: Path, observation_altitude_m: float,
                telescope_z_m: float) -> list[dict]:
    rows = read_csv(path)
    target = min(160, len(rows))
    indices = {round(i * (len(rows) - 1) / max(1, target - 1)) for i in range(target)}
    shower = []
    for index, row in enumerate(rows):
        if index not in indices:
            continue
        anchor = vec_m(row, "input_")
        direction = vec(row, "input_dir_")
        height = number(row, "emission_altitude_km") * 1000.0 - observation_altitude_m - telescope_z_m
        distance = (height - anchor[2]) / (-direction[2])
        emission = [anchor[i] - distance * direction[i] for i in range(3)]
        path_length_m = math.sqrt(sum((emission[i] - anchor[i]) ** 2 for i in range(3)))
        record_time_ns = number(row, "time_ns")
        shower.append({
            "source_bunch_index": int(row["source_bunch_index"]),
            "anchor": anchor, "direction": direction, "emission": emission,
            "emission_altitude_km": number(row, "emission_altitude_km"),
            "multiplicity": number(row, "multiplicity"),
            "time_ns": record_time_ns,
            "path_length_m": path_length_m,
            "emission_time_ns": record_time_ns - path_length_m / 0.299792458,
        })
    return shower


def percentile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    return ordered[round((len(ordered) - 1) * fraction)]


def annotate_shower_consistency(example: dict, raw_rows: list[dict[str, str]]) -> None:
    """Compare raw bunch directions/derived points with the selected shower axis."""
    provenance = example["provenance"]
    event = provenance["corsika_event"]
    telescope = provenance["telescope_eventio"]["position_nwu_m"]
    reference = provenance["corsika_input_reference"]
    azimuth = math.radians(event["arrival_azimuth_north_to_east_deg"])
    altitude = math.radians(event["altitude_deg"])
    sky = [
        math.cos(altitude) * math.cos(azimuth),
        -math.cos(altitude) * math.sin(azimuth),
        math.sin(altitude),
    ]
    core = [event["core_x_north_m"], event["core_y_west_m"], 0.0]
    angle_offsets: list[float] = []
    emission_residuals: list[float] = []
    anchor_residuals: list[float] = []

    def annotate(anchor: list[float], direction: list[float], emission: list[float]) -> dict:
        upward = [-value for value in direction]
        dot_value = max(-1.0, min(1.0, sum(upward[i] * sky[i] for i in range(3))))
        angle_deg = math.degrees(math.acos(dot_value))
        anchor_array = add(telescope, anchor)
        emission_array = add(telescope, emission)
        axis_t = (emission_array[2] - core[2]) / sky[2]
        axis_point = add(core, sky, axis_t)
        anchor_axis_t = (anchor_array[2] - core[2]) / sky[2]
        anchor_axis_point = add(core, sky, anchor_axis_t)
        emission_residual = math.hypot(
            emission_array[0] - axis_point[0], emission_array[1] - axis_point[1]
        )
        anchor_residual = math.hypot(
            anchor_array[0] - anchor_axis_point[0],
            anchor_array[1] - anchor_axis_point[1],
        )
        return {
            "arrival_axis_up": sky,
            "anchor_array_nwu": anchor_array,
            "emission_array_nwu": emission_array,
            "axis_point_same_height_nwu": axis_point,
            "direction_offset_deg": angle_deg,
            "emission_axis_residual_m": emission_residual,
            "anchor_axis_residual_m": anchor_residual,
        }

    for row in raw_rows:
        anchor = vec_m(row, "input_")
        direction = vec(row, "input_dir_")
        height = (
            number(row, "emission_altitude_km") * 1000.0
            - reference["observation_altitude_m"]
            - reference["telescope_z_m"]
        )
        distance = (height - anchor[2]) / (-direction[2])
        emission = [anchor[i] - distance * direction[i] for i in range(3)]
        result = annotate(anchor, direction, emission)
        angle_offsets.append(result["direction_offset_deg"])
        emission_residuals.append(result["emission_axis_residual_m"])
        anchor_residuals.append(result["anchor_axis_residual_m"])

    for bunch in example["shower"]:
        bunch.update(annotate(bunch["anchor"], bunch["direction"], bunch["emission"]))

    example["shower_consistency"] = {
        "checked_bunches": len(raw_rows),
        "direction_offset_deg_median": percentile(angle_offsets, 0.5),
        "direction_offset_deg_p95": percentile(angle_offsets, 0.95),
        "direction_offset_deg_max": max(angle_offsets),
        "emission_axis_residual_m_median": percentile(emission_residuals, 0.5),
        "emission_axis_residual_m_p95": percentile(emission_residuals, 0.95),
        "anchor_axis_residual_m_median": percentile(anchor_residuals, 0.5),
        "anchor_axis_residual_m_p95": percentile(anchor_residuals, 0.95),
        "interpretation": (
            "core is the shower-axis crossing at array z=0; photon anchors are individual "
            "Cherenkov arrivals relative to the telescope and are not expected to equal core"
        ),
    }


def validate_sources(ideal: list[dict], deformation: list[dict], examples: dict,
                     camera_geometry: list[dict]) -> dict:
    if [row["elevation_deg"] for row in deformation] != [float(x) for x in range(0, 91, 10)]:
        raise ValueError("deformation series must contain the exact 0..90 deg anchors")
    ideal_by_id = {row["id"]: row for row in ideal}
    for state in deformation:
        if len(state["facets"]) != len(ideal):
            raise ValueError("deformation facet count does not match ideal layout")
        for facet in state["facets"]:
            base = ideal_by_id[facet["id"]]
            for i in range(3):
                expected = facet["center"][i] - base["center"][i]
                if abs(expected - facet["delta"][i]) > 1e-12:
                    raise ValueError("deformation delta is not a direct center subtraction")
    pixel_ids = {row["id"] for row in camera_geometry}
    for name, example in examples.items():
        provenance = example["provenance"]
        event = provenance["corsika_event"]
        telescope = provenance["telescope_eventio"]
        if event["output_event_id"] != int(provenance["selection"]["event_id"]):
            raise ValueError(f"{name}: core metadata event does not match selected output event")
        if event["array_id"] != event["output_event_id"] % 100:
            raise ValueError(f"{name}: array id does not match event_array100 output id")
        if event["core_source"] != "negative_MC_TELOFF_array_offset":
            raise ValueError(f"{name}: displayed shower core is not sourced from MC_TELOFF")
        if telescope["position_nwu_m"] != [0.0, 0.0, 4.0]:
            raise ValueError(f"{name}: unexpected EventIO telescope position")
        consistency = example["shower_consistency"]
        if consistency["checked_bunches"] != len(example["raw_rows"]):
            raise ValueError(f"{name}: shower consistency did not use every raw bunch")
        if consistency["direction_offset_deg_max"] >= 2.0:
            raise ValueError(f"{name}: photon directions disagree with the shower header axis")
        raw_by_index = {int(row["source_bunch_index"]): row for row in example["raw_rows"]}
        for trace in example["traces"]:
            if abs(trace["surface"][0] - trace["u_m"]) > 1e-9 or abs(trace["surface"][1] - trace["v_m"]) > 1e-9:
                raise ValueError(f"{name}: trace surface x/y do not equal raw u/v")
            raw = raw_by_index[trace["source_bunch_index"]]
            expected_position = vec_m(raw, "input_")
            expected_direction = vec(raw, "input_dir_")
            if any(abs(a - b) > 1e-9 for a, b in zip(trace["raw_input"], expected_position)):
                raise ValueError(f"{name}: trace raw position does not match raw-input CSV")
            if any(abs(a - b) > 1e-9 for a, b in zip(trace["raw_input_dir"], expected_direction)):
                raise ValueError(f"{name}: trace raw direction does not match raw-input CSV")
        if len(example["output_points"]) != int(example["summary"]["hit_output_plane"]):
            raise ValueError(f"{name}: output-plane point count does not match summary")
        if any(row["pixel_id"] not in pixel_ids for row in example["camera_signal"]):
            raise ValueError(f"{name}: camera output references a missing pixel id")
    return {
        "status": "passed",
        "checks": [
            "10 exact deformation anchors and 54 facet ids",
            "delta center is direct deformed minus ideal subtraction",
            "both real examples keep whiteboard surface x/y equal to u/v",
            "both examples join raw input exactly by source_bunch_index",
            "full output-plane row counts equal each C++ summary",
            "camera output pixel ids exist in the 1616-pixel geometry",
            "global shower core and EventIO telescope position match the selected event metadata",
            "every raw bunch direction is within 2 deg of the selected shower arrival axis",
        ],
    }


def validate_parallel_cases(parallel_cases: dict) -> dict:
    if parallel_cases.get("status") != "ready":
        return {"status": "pending", "checks": []}
    expected_angles = [float(x) for x in range(0, 91, 10)]
    cases = parallel_cases.get("elevation_cases", [])
    actual_angles = [float(case["pointing"]["el_deg"]) for case in cases]
    if actual_angles != expected_angles:
        raise ValueError("parallel elevation cases must contain the exact 0..90 deg scan")
    four_cases = parallel_cases.get("four_direction_cases", [])
    expected_sky = [(0.0, 71.0), (0.0, 69.0), (-1.0, 70.0), (1.0, 70.0)]
    actual_sky = [
        (float(case["source_sky"]["az_deg"]), float(case["source_sky"]["el_deg"]))
        for case in four_cases
    ]
    if actual_sky != expected_sky:
        raise ValueError("parallel sky cases must be el=71/69 and az=-1/+1 around az=0/el=70")
    pointing_el = math.radians(70.0)
    basis_x = [0.0, 1.0, 0.0]
    basis_y = [-math.sin(pointing_el), 0.0, math.cos(pointing_el)]
    basis_z = [math.cos(pointing_el), 0.0, math.sin(pointing_el)]
    for case, (source_az_deg, source_el_deg) in zip(four_cases, expected_sky):
        source_az = math.radians(source_az_deg)
        source_el = math.radians(source_el_deg)
        propagation_global = [
            -math.cos(source_el) * math.cos(source_az),
            -math.cos(source_el) * math.sin(source_az),
            -math.sin(source_el),
        ]
        expected_local = [
            sum(propagation_global[i] * axis[i] for i in range(3))
            for axis in (basis_x, basis_y, basis_z)
        ]
        actual_local = case["beam"]["direction_local"]
        if any(abs(a - b) > 1e-9 for a, b in zip(actual_local, expected_local)):
            raise ValueError("parallel source sky angle does not reproduce source.beam_direction")
        summary = case["summary"]
        provenance = case["provenance"]
        if int(provenance["mark_only_rows"]) != int(summary["hit_output_before_obstruction"]):
            raise ValueError("parallel mark-only row count does not match pre-obstruction output")
        if int(provenance["full_output_rows"]) != int(summary["hit_output_plane"]):
            raise ValueError("parallel clear row count does not match physical output")
        if len(case.get("full_output_uv_m", [])) != int(summary["hit_output_plane"]):
            raise ValueError("parallel full output u/v payload is incomplete")
        if len(case.get("full_camera_hit_uv_m", [])) != int(summary["hit_camera"]):
            raise ValueError("parallel full camera-hit u/v payload is incomplete")
        if len(case.get("mirror_stats", [])) != 54:
            raise ValueError("parallel mirror statistics must cover all 54 facets")
    centroid_signs = [
        four_cases[0]["camera_summary"]["output_uv_centroid_m"][1] < 0,
        four_cases[1]["camera_summary"]["output_uv_centroid_m"][1] > 0,
        four_cases[2]["camera_summary"]["output_uv_centroid_m"][0] > 0,
        four_cases[3]["camera_summary"]["output_uv_centroid_m"][0] < 0,
    ]
    if not all(centroid_signs):
        raise ValueError("parallel camera centroids do not show the expected optical inversion")
    all_cases = [parallel_cases.get("baseline"), *four_cases, *cases]
    if any(case is None for case in all_cases):
        raise ValueError("parallel cases are marked ready but a case is missing")
    if any(case.get("validation", {}).get("status") != "passed" for case in all_cases):
        raise ValueError("a parallel case failed the raw surface/camera validation")
    for case in all_cases:
        summary = case["summary"]
        validation = case["validation"]
        if len(case.get("full_output_uv_m", [])) != int(summary["hit_output_plane"]):
            raise ValueError("parallel case does not embed every physical output u/v")
        if len(case.get("full_camera_hit_uv_m", [])) != int(summary["hit_camera"]):
            raise ValueError("parallel case does not embed every accepted camera-hit u/v")
        if validation["camera_pixel_photon_sum"] != int(summary["hit_camera"]):
            raise ValueError("parallel camera pixel sum does not equal the C++ hit count")
        if validation["max_surface_to_uv_error_m"] >= 1e-7:
            raise ValueError("parallel global-surface projection does not reproduce raw u/v")
        if validation["max_camera_xy_to_uv_error_m"] >= 1e-10:
            raise ValueError("parallel camera x/y do not reproduce raw u/v")
    if not parallel_cases.get("source_archive_sha256") or not parallel_cases.get("results_archive_sha256"):
        raise ValueError("parallel output is missing source/results archive provenance hashes")
    return {
        "status": "passed",
        "checks": [
            "baseline and exact 0..90 deg elevation cases come from captured C++ runs",
            "four sky sources are exact az/el inputs transformed to the logged local beam vectors",
            "mark-only and clear rows reproduce the C++ pre/post-obstruction counts",
            "camera centroids show upper/lower and left/right optical inversion",
            "raw global output points project back to the logged u/v coordinates",
            "camera x/y equal raw u/v and pixel sums equal the C++ hit_camera count",
            "all baseline, sky-direction, and elevation cases embed every physical output u/v and every accepted camera-hit u/v",
            "source and result archives retain SHA-256 provenance",
        ],
    }


def compact(value):
    if isinstance(value, float):
        return round(value, 10)
    if isinstance(value, list):
        return [compact(x) for x in value]
    if isinstance(value, dict):
        return {k: compact(v) for k, v in value.items()}
    return value


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True)
    parser.add_argument("--ideal-mirror", default="configs/mirror_1229_facets.csv")
    parser.add_argument("--deformation-series", default="configs/mirror_1229_elevation_series.csv")
    parser.add_argument("--obstruction-primitives", default="configs/obstructions/raytrace_final_structure_primitives.csv")
    parser.add_argument("--camera-geometry", default="configs/cameras/new_camera_pixels.csv")
    parser.add_argument("--north-trace-csv", default="docs/assets/data/corsika-north-example-rays.csv")
    parser.add_argument("--north-raw-input-csv", default="docs/assets/data/corsika-north-example-raw-input.csv")
    parser.add_argument("--north-output-plane-csv", default="docs/assets/data/corsika-north-output-plane.csv")
    parser.add_argument("--north-summary-csv", default="docs/assets/data/corsika-north-example-summary.csv")
    parser.add_argument("--north-provenance", default="docs/assets/data/corsika-north-example-provenance.json")
    parser.add_argument("--aligned-trace-csv", default="docs/assets/data/corsika-aligned-example-rays.csv")
    parser.add_argument("--aligned-raw-input-csv", default="docs/assets/data/corsika-aligned-example-raw-input.csv")
    parser.add_argument("--aligned-output-plane-csv", default="docs/assets/data/corsika-aligned-output-plane.csv")
    parser.add_argument("--aligned-summary-csv", default="docs/assets/data/corsika-aligned-example-summary.csv")
    parser.add_argument("--aligned-provenance", default="docs/assets/data/corsika-aligned-example-provenance.json")
    parser.add_argument("--camera-pixel-csv", default="docs/assets/data/corsika-aligned-camera-pixels.csv")
    parser.add_argument("--camera-summary-csv", default="docs/assets/data/corsika-aligned-camera-summary.csv")
    parser.add_argument("--north-camera-summary-csv", default="docs/assets/data/corsika-north-camera-summary.csv")
    parser.add_argument(
        "--parallel-cases",
        default="docs/assets/data/coordinate-parallel-cases.json",
        help="compact real run_optical_sim baseline/elevation cases",
    )
    parser.add_argument(
        "--event1909-case",
        default="docs/assets/data/corsika-event1909-coordinate-case.json",
        help="compact real prod1 EventIO event 1909 array/camera/ray payload",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    root = Path.cwd()
    ideal, deformation, primitives = read_geometry(
        root / args.ideal_mirror, root / args.deformation_series,
        root / args.obstruction_primitives,
    )
    north_provenance = json.loads((root / args.north_provenance).read_text(encoding="utf-8"))
    aligned_provenance = json.loads((root / args.aligned_provenance).read_text(encoding="utf-8"))
    camera_geometry = [{
        "id": int(row["id"]), "u": number(row, "x_m"), "v": number(row, "y_m"),
        "size": number(row, "size_m"), "shape": row.get("shape", "square"),
    } for row in read_csv(root / args.camera_geometry)]
    camera_signal = [{
        "pixel_id": int(row["pixel_id"]), "photon_count": int(row["photon_count"]),
        "pe": number(row, "pe"), "time_mean_ns": number(row, "time_mean_ns"),
    } for row in read_csv(root / args.camera_pixel_csv)]
    camera_summary = read_csv(root / args.camera_summary_csv)[0]
    north_camera_summary = read_csv(root / args.north_camera_summary_csv)[0]
    parallel_path = root / args.parallel_cases
    parallel_cases = (
        json.loads(parallel_path.read_text(encoding="utf-8"))
        if parallel_path.exists()
        else {
            "status": "pending_current_cpp_run",
            "reason": (
                "必须先编译并运行当前 C++ 源码，才能把平行光光路和相机图标记为程序真实输出；"
                "目前没有使用旧二进制或其他事例替代。"
            ),
            "baseline": None,
            "elevation_cases": [],
        }
    )
    event1909 = json.loads((root / args.event1909_case).read_text(encoding="utf-8"))
    if event1909.get("validation", {}).get("status") != "passed":
        raise ValueError("event 1909 coordinate payload did not pass its source checks")
    examples = {}
    example_specs = {
        "north": {
            "title": "望远镜正北",
            "pointing": {"az_deg": 0.0, "el_deg": 70.0},
            "trace": args.north_trace_csv, "raw": args.north_raw_input_csv,
            "output": args.north_output_plane_csv, "summary": args.north_summary_csv,
            "provenance": north_provenance, "camera_summary": north_camera_summary,
            "camera_signal": [],
        },
        "aligned": {
            "title": "望远镜对准簇射",
            "pointing": {"az_deg": 300.027133, "el_deg": 88.282787},
            "trace": args.aligned_trace_csv, "raw": args.aligned_raw_input_csv,
            "output": args.aligned_output_plane_csv, "summary": args.aligned_summary_csv,
            "provenance": aligned_provenance, "camera_summary": camera_summary,
            "camera_signal": camera_signal,
        },
    }
    for name, spec in example_specs.items():
        raw_rows = read_csv(root / spec["raw"])
        reference = spec["provenance"]["corsika_input_reference"]
        examples[name] = {
            "title": spec["title"], "pointing": spec["pointing"],
            "traces": read_traces(root / spec["trace"]),
            "output_points": read_output_plane(root / spec["output"]),
            "summary": read_csv(root / spec["summary"])[0],
            "camera_summary": spec["camera_summary"],
            "camera_signal": spec["camera_signal"],
            "provenance": spec["provenance"],
            "shower": read_shower(root / spec["raw"],
                                   reference["observation_altitude_m"],
                                   reference["telescope_z_m"]),
            "raw_rows": raw_rows,
        }
        annotate_shower_consistency(examples[name], raw_rows)
    validation = validate_sources(ideal, deformation, examples, camera_geometry)
    validation["parallel"] = validate_parallel_cases(parallel_cases)
    for example in examples.values():
        del example["raw_rows"]
    data = compact({
        "ideal": ideal,
        "deformation": deformation,
        "primitives": primitives,
        "examples": examples,
        "parallel": parallel_cases,
        "event1909": event1909,
        "camera_geometry": camera_geometry,
        "camera_signal": camera_signal,
        "validation": validation,
        "provenance": north_provenance,
        "source_files": {
            "ideal": args.ideal_mirror, "deformation": args.deformation_series,
            "obstruction": args.obstruction_primitives, "camera": args.camera_geometry,
            "north_trace": args.north_trace_csv,
            "north_raw_input": args.north_raw_input_csv,
            "north_output_plane": args.north_output_plane_csv,
            "aligned_trace": args.aligned_trace_csv,
            "aligned_raw_input": args.aligned_raw_input_csv,
            "aligned_output_plane": args.aligned_output_plane_csv,
            "camera_signal": args.camera_pixel_csv,
            "event1909": args.event1909_case,
        },
    })

    template_path = Path(__file__).with_name("coordinate_diagnostics_template.html")
    html = template_path.read_text(encoding="utf-8").replace(
        "__DATA__", json.dumps(data, ensure_ascii=False, separators=(",", ":"))
    )
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(html, encoding="utf-8")
    print(f"Saved coordinate diagnostic HTML = {output}")


if __name__ == "__main__":
    main()
