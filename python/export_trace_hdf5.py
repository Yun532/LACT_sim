#!/usr/bin/env python3
"""Pack LACT_sim CSV outputs and static geometry into one HDF5 file.

The exporter keeps CSV as the debugging format, but creates a self-contained
file for large-scale analysis:

  - static geometry is stored once: camera pixels, mirror facets, telescope
  - image data can be sparse, dense, or both
  - each image is one event/telescope pair
"""

import argparse
import csv
from pathlib import Path

import h5py
import numpy as np

from config_io import (
    expand_component_config,
    load_facets_from_config,
    parse_vec3,
    resolve_workspace_path,
)


def read_camera_pixels(path):
    rows = []
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            rows.append((
                int(row["id"]),
                float(row["x_m"]),
                float(row["y_m"]),
                row.get("shape", "Square").strip(),
                float(row["size_m"]),
            ))
    return rows


def shape_code(shape):
    text = shape.strip().lower()
    if text in {"square", "sq"}:
        return 1
    if text in {"hex", "hexagon", "hexagonal"}:
        return 2
    if text in {"circle", "circular"}:
        return 3
    return 0


def read_pixel_csv(path):
    rows = []
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            rows.append({
                "event_id": int(row["event_id"]),
                "telescope_id": int(row["telescope_id"]),
                "pixel_id": int(row["pixel_id"]),
                "photon_count": int(float(row["photon_count"])),
                "pe": float(row.get("pe", row.get("signal", 0.0))),
                "signal": float(row.get("signal", row.get("pe", 0.0))),
                "time_mean_ns": float(row.get("time_mean_ns", 0.0)),
                "time_rms_ns": float(row.get("time_rms_ns", 0.0)),
            })
    return rows


def read_summary_csv(path):
    if not path:
        return {}
    p = Path(path)
    if not p.exists():
        return {}
    summary = {}
    with p.open(newline="") as f:
        for row in csv.DictReader(f):
            key = (int(row["event_id"]), int(row["telescope_id"]))
            summary[key] = row
    return summary


def image_keys(pixel_rows, summary):
    keys = set(summary)
    keys.update((row["event_id"], row["telescope_id"]) for row in pixel_rows)
    return sorted(keys)


def make_image_index(pixel_rows, summary):
    by_image = {}
    for row in pixel_rows:
        by_image.setdefault((row["event_id"], row["telescope_id"]), []).append(row)

    keys = image_keys(pixel_rows, summary)
    sparse = []
    index = []
    start = 0
    for image_index, key in enumerate(keys):
        rows = sorted(by_image.get(key, []), key=lambda r: r["pixel_id"])
        count = len(rows)
        total_signal = sum(r["signal"] for r in rows)
        total_pe = sum(r["pe"] for r in rows)
        total_photons = sum(r["photon_count"] for r in rows)

        srow = summary.get(key, {})
        if srow:
            total_signal = float(srow.get("signal", total_signal))
            total_pe = float(srow.get("pe", total_pe))
            total_photons = float(srow.get("hit_camera", total_photons))

        index.append((
            image_index,
            key[0],
            key[1],
            start,
            count,
            total_photons,
            total_pe,
            total_signal,
            float(srow.get("time_mean_ns", 0.0)) if srow else 0.0,
            float(srow.get("time_rms_ns", 0.0)) if srow else 0.0,
        ))
        for r in rows:
            sparse.append((
                r["pixel_id"],
                r["photon_count"],
                r["pe"],
                r["signal"],
                r["time_mean_ns"],
                r["time_rms_ns"],
            ))
        start += count
    return index, sparse


def write_string_dataset(group, name, text):
    group.create_dataset(name, data=text.encode("utf-8"))


def write_config_group(h5, config_path, cfg, component_paths):
    g = h5.create_group("config")
    if config_path:
        config_path = Path(config_path)
        g.attrs["main_config_path"] = str(config_path)
        if config_path.exists():
            write_string_dataset(g, "main_config_text", config_path.read_text())
    for key, value in sorted(cfg.items()):
        g.attrs[key] = str(value)
    components = g.create_group("components")
    for label, path in sorted(component_paths.items()):
        components.attrs[label] = str(path)
        p = Path(path)
        if p.exists():
            write_string_dataset(components, f"{label}_text", p.read_text())


def write_telescope_group(h5, cfg, image_index):
    telescope_ids = sorted(set(int(row[2]) for row in image_index))
    dtype = np.dtype([
        ("telescope_id", "i4"),
        ("x_m", "f8"),
        ("y_m", "f8"),
        ("z_m", "f8"),
        ("pointing_az_deg", "f8"),
        ("pointing_el_deg", "f8"),
        ("focal_length_m", "f8"),
    ])
    pos = parse_vec3(cfg.get("telescope.position_m"), [0, 0, 0])
    table = np.zeros(len(telescope_ids), dtype=dtype)
    for i, tid in enumerate(telescope_ids):
        table[i] = (
            tid,
            pos[0],
            pos[1],
            pos[2],
            float(cfg.get("telescope.pointing_az_deg", 0.0)),
            float(cfg.get("telescope.pointing_el_deg", 90.0)),
            float(cfg.get("telescope.focal_length_m", 8.0)),
        )
    h5.create_group("telescopes").create_dataset("table", data=table)


def write_camera_group(h5, camera_pixels):
    dtype = np.dtype([
        ("pixel_id", "i4"),
        ("x_m", "f4"),
        ("y_m", "f4"),
        ("size_m", "f4"),
        ("shape_code", "i2"),
    ])
    data = np.zeros(len(camera_pixels), dtype=dtype)
    for i, (pid, x, y, shape, size) in enumerate(camera_pixels):
        data[i] = (pid, x, y, size, shape_code(shape))
    g = h5.create_group("camera")
    g.attrs["shape_code_map"] = "0=unknown,1=square,2=hexagon,3=circular"
    g.create_dataset("pixels", data=data)


def write_mirror_group(h5, facets):
    dtype = np.dtype([
        ("mirror_id", "i4"),
        ("center_x_m", "f4"),
        ("center_y_m", "f4"),
        ("center_z_m", "f4"),
        ("normal_x", "f4"),
        ("normal_y", "f4"),
        ("normal_z", "f4"),
        ("radius_of_curvature_m", "f4"),
        ("size1_m", "f4"),
        ("size2_m", "f4"),
        ("aperture_rotation_rad", "f4"),
        ("shape_code", "i2"),
    ])
    data = np.zeros(len(facets), dtype=dtype)
    for i, f in enumerate(facets):
        data[i] = (
            int(f["id"]),
            *np.asarray(f["center"], dtype=float),
            *np.asarray(f["normal"], dtype=float),
            float(f.get("radius_of_curvature", 0.0)),
            float(f.get("size1", 0.0)),
            float(f.get("size2", 0.0)),
            float(f.get("rotation", 0.0)),
            shape_code(f.get("shape", "")),
        )
    g = h5.create_group("mirrors")
    g.attrs["shape_code_map"] = "0=unknown,1=square,2=hexagon,3=circular"
    g.create_dataset("facets", data=data)


def write_events_group(h5, image_index, event_id_mode):
    event_ids = sorted(set(int(row[1]) for row in image_index))
    dtype = np.dtype([
        ("event_index", "i4"),
        ("event_id", "i8"),
        ("shower_event_id", "i8"),
        ("array_id", "i4"),
    ])
    data = np.zeros(len(event_ids), dtype=dtype)
    is_array100 = event_id_mode.lower() in {"event_array100", "runid"}
    for i, event_id in enumerate(event_ids):
        shower = event_id // 100 if is_array100 else event_id
        array = event_id % 100 if is_array100 else 0
        data[i] = (i, event_id, shower, array)
    h5.create_group("events").create_dataset("table", data=data)


def write_images_group(h5, image_index, sparse_rows, camera_pixels, storage, compression):
    g = h5.create_group("images")
    index_dtype = np.dtype([
        ("image_index", "i4"),
        ("event_id", "i8"),
        ("telescope_id", "i4"),
        ("start", "i8"),
        ("count", "i4"),
        ("total_photons", "f8"),
        ("total_pe", "f8"),
        ("total_signal", "f8"),
        ("time_mean_ns", "f4"),
        ("time_rms_ns", "f4"),
    ])
    sparse_dtype = np.dtype([
        ("pixel_id", "i4"),
        ("photon_count", "i4"),
        ("pe", "f4"),
        ("signal", "f4"),
        ("time_mean_ns", "f4"),
        ("time_rms_ns", "f4"),
    ])
    index_arr = np.array(image_index, dtype=index_dtype)
    sparse_arr = np.array(sparse_rows, dtype=sparse_dtype)
    g.create_dataset("index", data=index_arr, compression=compression)

    if storage in {"sparse", "both"}:
        sg = g.create_group("sparse")
        sg.create_dataset("pixels", data=sparse_arr, compression=compression)

    if storage in {"dense", "both"}:
        pixel_ids = np.array([p[0] for p in camera_pixels], dtype=np.int32)
        pixel_to_col = {int(pid): i for i, pid in enumerate(pixel_ids)}
        n_images = len(index_arr)
        n_pixels = len(pixel_ids)
        signal = np.zeros((n_images, n_pixels), dtype=np.float32)
        pe = np.zeros((n_images, n_pixels), dtype=np.float32)
        photon_count = np.zeros((n_images, n_pixels), dtype=np.int32)
        for image in index_arr:
            img = int(image["image_index"])
            start = int(image["start"])
            count = int(image["count"])
            for row in sparse_arr[start:start + count]:
                col = pixel_to_col.get(int(row["pixel_id"]))
                if col is None:
                    continue
                signal[img, col] = row["signal"]
                pe[img, col] = row["pe"]
                photon_count[img, col] = row["photon_count"]
        dg = g.create_group("dense")
        dg.create_dataset("pixel_id_axis", data=pixel_ids, compression=compression)
        dg.create_dataset("signal", data=signal, compression=compression, chunks=True)
        dg.create_dataset("pe", data=pe, compression=compression, chunks=True)
        dg.create_dataset("photon_count", data=photon_count, compression=compression, chunks=True)


def main():
    parser = argparse.ArgumentParser(description="Export LACT_sim camera CSV output to self-contained HDF5.")
    parser.add_argument("--pixel-csv", required=True, help="camera pixel aggregate CSV")
    parser.add_argument("--output", required=True, help="output HDF5 path")
    parser.add_argument("--summary-csv", default=None, help="optional run summary CSV")
    parser.add_argument("--config", default=None, help="run cfg; used to embed geometry/config")
    parser.add_argument("--camera-csv", default=None, help="camera pixel CSV if not found in config")
    parser.add_argument("--storage", choices=("sparse", "dense", "both"), default="sparse")
    parser.add_argument("--compression", default="gzip", help="HDF5 compression filter, or none")
    args = parser.parse_args()

    cfg = {}
    component_paths = {}
    config_path = Path(args.config).resolve() if args.config else None
    if config_path:
        cfg, component_paths = expand_component_config(config_path)

    pixel_csv = Path(args.pixel_csv)
    pixel_rows = read_pixel_csv(pixel_csv)
    summary = read_summary_csv(args.summary_csv)
    image_index, sparse_rows = make_image_index(pixel_rows, summary)

    camera_csv = args.camera_csv or cfg.get("camera.csv_path")
    if not camera_csv:
        raise SystemExit("--camera-csv is required when config does not define camera.csv_path")
    camera_path = resolve_workspace_path(config_path or pixel_csv, camera_csv)
    camera_pixels = read_camera_pixels(camera_path)

    facets = []
    if config_path and cfg.get("mirror.mode"):
        facets = load_facets_from_config(config_path, cfg)

    compression = None if args.compression.lower() in {"none", "off", "false"} else args.compression
    out = Path(args.output)
    if out.parent:
        out.parent.mkdir(parents=True, exist_ok=True)

    with h5py.File(out, "w") as h5:
        h5.attrs["format"] = "LACT_sim trace HDF5"
        h5.attrs["format_version"] = "0.1"
        h5.attrs["image_storage"] = args.storage
        h5.attrs["source_pixel_csv"] = str(pixel_csv)
        if args.summary_csv:
            h5.attrs["source_summary_csv"] = str(args.summary_csv)
        if config_path:
            write_config_group(h5, config_path, cfg, component_paths)
        write_camera_group(h5, camera_pixels)
        if facets:
            write_mirror_group(h5, facets)
        write_events_group(h5, image_index, cfg.get("source.event_id_mode", "event"))
        write_telescope_group(h5, cfg, image_index)
        write_images_group(h5, image_index, sparse_rows, camera_pixels, args.storage, compression)

    print(f"Saved HDF5 = {out}")
    print(f"Images = {len(image_index)}")
    print(f"Sparse pixel rows = {len(sparse_rows)}")
    print(f"Camera pixels = {len(camera_pixels)}")
    print(f"Storage = {args.storage}")


if __name__ == "__main__":
    main()
