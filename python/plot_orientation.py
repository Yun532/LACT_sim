"""Shared 2D orientation helpers for LACT plots.

The display convention implemented here is:

  display +y = projection of global +z/up onto the plotted camera/mirror plane.

If the projection is degenerate, the original plot y axis is kept.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from config_io import (
    expand_component_config,
    parse_vec3,
    rotate_local_vector,
    source_coordinate_frame_name_from_config,
    source_telescope_frame_from_config,
    telescope_frame_from_config,
)


def display_basis_from_up_vector(up_in_plot_xy: np.ndarray) -> tuple[np.ndarray, np.ndarray, bool]:
    y_axis = np.asarray(up_in_plot_xy, dtype=float)
    norm = float(np.linalg.norm(y_axis))
    if norm < 1e-12:
        return np.array([1.0, 0.0]), np.array([0.0, 1.0]), False
    y_axis /= norm
    x_axis = np.array([y_axis[1], -y_axis[0]])
    return x_axis, y_axis, True


def basis_from_optical_config(config_path: str | Path) -> tuple[np.ndarray, np.ndarray, bool]:
    cfg, _ = expand_component_config(Path(config_path).resolve())
    frame = telescope_frame_from_config(cfg)
    local_u = parse_vec3(cfg.get("output.plane_u_axis"), [1.0, 0.0, 0.0])
    local_v = parse_vec3(cfg.get("output.plane_v_axis"), [0.0, 1.0, 0.0])
    global_u = rotate_local_vector(local_u, frame)
    global_v = rotate_local_vector(local_v, frame)
    global_up = np.array([0.0, 0.0, 1.0])
    up_in_uv = np.array([np.dot(global_up, global_u), np.dot(global_up, global_v)])
    return display_basis_from_up_vector(up_in_uv)


def basis_from_config(config_path: str | Path) -> tuple[np.ndarray, np.ndarray, bool]:
    """Use source adapters for external frames and layout pointing for local input."""
    cfg, _ = expand_component_config(Path(config_path).resolve())
    source_frame = source_coordinate_frame_name_from_config(cfg)
    frame = (
        telescope_frame_from_config(cfg)
        if source_frame == "telescope_local"
        else source_telescope_frame_from_config(cfg)
    )
    local_u = parse_vec3(cfg.get("output.plane_u_axis"), [1.0, 0.0, 0.0])
    local_v = parse_vec3(cfg.get("output.plane_v_axis"), [0.0, 1.0, 0.0])
    global_u = rotate_local_vector(local_u, frame)
    global_v = rotate_local_vector(local_v, frame)
    global_up = np.array([0.0, 0.0, 1.0])
    up_in_uv = np.array([np.dot(global_up, global_u), np.dot(global_up, global_v)])
    return display_basis_from_up_vector(up_in_uv)


def basis_from_corsika_pointing(az_deg: float, el_deg: float) -> tuple[np.ndarray, np.ndarray, bool]:
    """Return display basis for CORSIKA/EventIO HDF5 camera images.

    CORSIKA uses x=magnetic North, y=West, z=Up. The local x/y basis here follows
    the `corsikaIactFrame()` transform used by run_corsika_trace.
    """

    az = np.radians(float(az_deg))
    el = np.radians(float(el_deg))
    x_axis = np.array([
        -np.sin(el) * np.cos(az),
        np.sin(el) * np.sin(az),
        np.cos(el),
    ])
    y_axis = np.array([
        -np.sin(az),
        -np.cos(az),
        0.0,
    ])
    global_up = np.array([0.0, 0.0, 1.0])
    up_in_xy = np.array([np.dot(global_up, x_axis), np.dot(global_up, y_axis)])
    return display_basis_from_up_vector(up_in_xy)


def orient_xy(x, y, x_axis: np.ndarray, y_axis: np.ndarray):
    x_arr = np.asarray(x, dtype=float)
    y_arr = np.asarray(y, dtype=float)
    return x_axis[0] * x_arr + x_axis[1] * y_arr, y_axis[0] * x_arr + y_axis[1] * y_arr


def orient_points(points: np.ndarray, x_axis: np.ndarray, y_axis: np.ndarray) -> np.ndarray:
    pts = np.asarray(points, dtype=float)
    out_x, out_y = orient_xy(pts[:, 0], pts[:, 1], x_axis, y_axis)
    return np.column_stack([out_x, out_y])
