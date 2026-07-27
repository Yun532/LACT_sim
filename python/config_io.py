import csv
import math
from pathlib import Path

import numpy as np


def read_key_value_config(path):
    values = {}
    if path is None:
        return values
    with open(path, encoding="utf-8-sig") as f:
        for line in f:
            line = line.split("#", 1)[0].strip()
            if not line or "=" not in line:
                continue
            key, value = line.split("=", 1)
            values[key.strip().lower()] = value.strip()
    return values


def _scoped_key(key, prefix):
    scoped_prefixes = (
        "telescope.",
        "mirror.",
        "source.",
        "output.",
        "camera.",
        "sipm.",
        "electronics.",
        "efficiency.",
        "atmosphere.",
        "nsb.",
        "trigger.",
        "error.",
        "obstruction.",
        "dish.",
        "facet.",
    )
    if key.startswith(scoped_prefixes):
        return key
    return prefix + key


def _resolve_include_path(main_config_path, include_path):
    path = Path(include_path)
    if path.is_absolute():
        return path
    return main_config_path.parent / path


def _is_include_config_key(key):
    return key in {
        "mirror.config",
        "source.config",
        "output.config",
        "camera.config",
        "sipm.config",
        "electronics.config",
        "efficiency.config",
        "atmosphere.config",
        "nsb.config",
        "trigger.config",
        "error.config",
        "obstruction.config",
    }


def expand_component_config(main_config_path):
    main_config_path = Path(main_config_path).resolve()
    main_cfg = read_key_value_config(main_config_path)
    expanded = {}
    component_paths = {}
    assembly_cfg = {}

    telescope_include = main_cfg.get("telescope.config")
    if telescope_include:
        telescope_path = _resolve_include_path(main_config_path, telescope_include)
        component_paths["telescope"] = telescope_path
        telescope_cfg = read_key_value_config(telescope_path)
        for key, value in telescope_cfg.items():
            scoped = _scoped_key(key, "telescope.")
            if _is_include_config_key(scoped):
                value = str(_resolve_include_path(telescope_path, value))
            assembly_cfg[scoped] = value

    assembly_cfg.update(main_cfg)

    for include_key, prefix, label in (
        ("mirror.config", "mirror.", "mirror"),
        ("source.config", "source.", "source"),
        ("output.config", "output.", "output"),
        ("camera.config", "camera.", "camera"),
        ("sipm.config", "sipm.", "sipm"),
        ("electronics.config", "electronics.", "electronics"),
        ("efficiency.config", "efficiency.", "efficiency"),
        ("atmosphere.config", "atmosphere.", "atmosphere"),
        ("nsb.config", "nsb.", "nsb"),
        ("trigger.config", "trigger.", "trigger"),
        ("error.config", "error.", "error"),
        ("obstruction.config", "obstruction.", "obstruction"),
    ):
        include = assembly_cfg.get(include_key)
        if not include:
            continue
        component_path = _resolve_include_path(main_config_path, include)
        component_paths[label] = component_path
        component_cfg = read_key_value_config(component_path)
        for key, value in component_cfg.items():
            scoped = _scoped_key(key, prefix)
            if scoped == "error.structural_deformation_config":
                value = str(_resolve_include_path(component_path, value))
            if scoped in {"obstruction.mask_csv", "obstruction.primitives_csv"}:
                value = str(_resolve_include_path(component_path, value))
            expanded[scoped] = value

    expanded.update(assembly_cfg)
    return expanded, component_paths


def workspace_root_from_config(config_path):
    config_path = Path(config_path).resolve()
    for parent in (config_path.parent, *config_path.parents):
        if parent.name == "configs":
            return parent.parent
    return Path.cwd()


def resolve_workspace_path(config_path, value):
    path = Path(value)
    if path.is_absolute():
        return path
    return workspace_root_from_config(config_path) / path


def parse_vec3(text, fallback=None):
    if text is None:
        if fallback is None:
            raise ValueError("missing vector value")
        return np.array(fallback, dtype=float)
    parts = [float(x.strip()) for x in text.split(",")]
    if len(parts) != 3:
        raise ValueError(f"expected 3 comma-separated values, got: {text}")
    return np.array(parts, dtype=float)


def parse_float_list(text):
    return [float(x.strip()) for x in text.split(",") if x.strip()]


def format_angle_token(angle_deg):
    if float(angle_deg).is_integer():
        return str(int(angle_deg))
    return f"{float(angle_deg):g}".replace(".", "p")


def expand_angle_pattern(pattern, angle_deg):
    return pattern.replace("{angle}", format_angle_token(angle_deg))


def _normalize(v):
    v = np.asarray(v, dtype=float)
    n = np.linalg.norm(v)
    if n <= 0.0:
        raise ValueError("cannot normalize zero vector")
    return v / n


def _slerp(a, b, t):
    a = _normalize(a)
    b = _normalize(b)
    dot = float(np.clip(np.dot(a, b), -1.0, 1.0))
    if dot > 0.999999:
        return _normalize((1.0 - t) * a + t * b)
    if dot < -0.999999:
        return _normalize((1.0 - t) * a - t * b)
    theta = math.acos(dot)
    s = math.sin(theta)
    return _normalize(
        a * (math.sin((1.0 - t) * theta) / s) +
        b * (math.sin(t * theta) / s)
    )


def _read_facets_csv(path):
    facets = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            facets.append({
                "id": int(row["id"]),
                "center": np.array([
                    float(row["center_x"]),
                    float(row["center_y"]),
                    float(row["center_z"]),
                ], dtype=float),
                "normal": _normalize([
                    float(row["normal_x"]),
                    float(row["normal_y"]),
                    float(row["normal_z"]),
                ]),
                "shape": row["aperture_shape"].strip().lower(),
                "size1": float(row["size1"]),
                "size2": float(row.get("size2", 0.0) or 0.0),
                "rotation": float(row.get("aperture_rotation_rad", 0.0) or 0.0),
                "surface_type": row.get("surface_type", "Spherical").strip(),
                "radius_of_curvature": float(row.get("radius_of_curvature", 16.0) or 16.0),
            })
    return facets


def _parse_vec3_text(text, fallback):
    if text is None:
        return np.array(fallback, dtype=float)
    parts = [part.strip() for part in text.split(",") if part.strip()]
    if len(parts) != 3:
        raise ValueError(f"expected 3 vector components, got: {text}")
    return np.array([float(part) for part in parts], dtype=float)


def _design_normal_to_focus(center, focus):
    to_source = np.array([0.0, 0.0, 1.0], dtype=float)
    to_focus = _normalize(focus - center)
    return _normalize(to_source + to_focus)


def _build_generated_facets(cfg):
    dish_type = cfg.get("dish.type", "DaviesCotton").strip().lower()
    telescope_focal_length = float(cfg.get("dish.telescope_focal_length", 5.0))
    dish_shape_length = float(cfg.get("dish.dish_shape_length", telescope_focal_length))
    dish_radius = float(cfg.get("dish.dish_radius", 2.0))
    vertex = _parse_vec3_text(cfg.get("dish.vertex"), [0.0, 0.0, 0.0])

    facet_spacing = float(cfg.get("facet.spacing", 0.55))
    facet_radius = float(cfg.get("facet.radius", 0.22))
    shape = cfg.get("facet.aperture_shape", "Circular").strip().lower()
    surface_type = cfg.get("facet.surface_type", "Spherical").strip()
    focus = vertex + np.array([0.0, 0.0, telescope_focal_length], dtype=float)

    if facet_spacing <= 0.0:
        raise ValueError("facet.spacing must be > 0 for generated mirrors")

    facets = []
    nmax = int(math.ceil(dish_radius / facet_spacing))
    for iy in range(-nmax, nmax + 1):
        for ix in range(-nmax, nmax + 1):
            x = ix * facet_spacing
            y = iy * facet_spacing
            if math.hypot(x, y) > dish_radius:
                continue

            if dish_type in {"daviescotton", "davies-cotton", "dc"}:
                z = vertex[2]
                if x * x + y * y < dish_shape_length * dish_shape_length:
                    z = focus[2] - math.sqrt(dish_shape_length * dish_shape_length - x * x - y * y)
                center = vertex + np.array([x, y, z - vertex[2]], dtype=float)
                radius_of_curvature = 2.0 * telescope_focal_length
            elif dish_type in {"parabolic", "paraboloid"}:
                z = vertex[2] + (x * x + y * y) / (4.0 * dish_shape_length)
                center = np.array([vertex[0] + x, vertex[1] + y, z], dtype=float)
                radius_of_curvature = 2.0 * float(np.linalg.norm(focus - center))
            else:
                raise ValueError(f"unsupported dish.type for python plotting: {cfg.get('dish.type')}")

            facets.append({
                "id": len(facets),
                "center": center,
                "normal": _design_normal_to_focus(center, focus),
                "shape": shape,
                "size1": facet_radius,
                "size2": 0.0,
                "rotation": 0.0,
                "surface_type": surface_type,
                "radius_of_curvature": radius_of_curvature,
            })
    return facets


def load_facets_csv(path):
    return _read_facets_csv(Path(path))


def _read_raw1229_states(path, swap_xy=True):
    states = []
    with open(path, newline="") as f:
        reader = csv.reader(f)
        for row in reader:
            if not row:
                continue
            values = [float(cell) for cell in row]
            if len(values) != 8:
                raise ValueError(f"raw 1229 file {path} expected 8 columns, got {len(values)}")
            center = np.array(values[2:5], dtype=float) * 0.001
            curvature_center = np.array(values[5:8], dtype=float) * 0.001
            if swap_xy:
                center[[0, 1]] = center[[1, 0]]
                curvature_center[[0, 1]] = curvature_center[[1, 0]]
            states.append({
                "cell_index": int(round(values[0])),
                "ring_index": int(round(values[1])),
                "center": center,
                "normal": _normalize(curvature_center - center),
            })
    if not states:
        raise ValueError(f"raw 1229 mirror file is empty: {path}")
    return states


def _read_elevation_series_csv(path):
    by_angle = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            angle = float(row["elevation_deg"])
            state = {
                "id": int(row["id"]),
                "center": np.array([
                    float(row["center_x"]),
                    float(row["center_y"]),
                    float(row["center_z"]),
                ], dtype=float),
                "normal": _normalize([
                    float(row["normal_x"]),
                    float(row["normal_y"]),
                    float(row["normal_z"]),
                ]),
            }
            if row.get("surface_type"):
                state["surface_type"] = row["surface_type"].strip()
            if row.get("radius_of_curvature"):
                state["radius_of_curvature"] = float(row["radius_of_curvature"])
            if row.get("aperture_shape"):
                state["shape"] = row["aperture_shape"].strip().lower()
            if row.get("size1"):
                state["size1"] = float(row["size1"])
            if row.get("size2"):
                state["size2"] = float(row["size2"])
            if row.get("aperture_rotation_rad"):
                state["rotation"] = float(row["aperture_rotation_rad"])
            by_angle.setdefault(angle, []).append(state)
    if not by_angle:
        raise ValueError(f"elevation-series CSV is empty: {path}")
    for angle in by_angle:
        by_angle[angle].sort(key=lambda row: row["id"])
    return by_angle


def load_facets_from_scoped_mirror_cfg(cfg_path, cfg):
    mode = cfg.get("mirror.mode", "csv").strip().lower()
    if mode == "csv":
        mirror_path = resolve_workspace_path(cfg_path, cfg["mirror.csv_path"])
        return _read_facets_csv(mirror_path)

    if mode in {"generated", "auto"}:
        return _build_generated_facets(cfg)

    if mode not in {"elevation_series", "series"}:
        raise ValueError(f"unsupported mirror.mode for python plotting: {mode}")

    fmt = cfg.get("mirror.series_format", "raw1229").strip().lower()
    if fmt not in {"raw1229", "facet_csv", "csv"}:
        raise ValueError(f"unsupported mirror.series_format for python plotting: {fmt}")

    elevation = float(cfg.get("mirror.series_elevation_deg", cfg.get("telescope.pointing_el_deg", 90.0)))

    shape = cfg.get("mirror.aperture_shape", "Hexagon").strip().lower()
    size1 = float(cfg.get("mirror.size1", 0.8))
    size2 = float(cfg.get("mirror.size2", 0.0))
    rotation = float(cfg.get("mirror.aperture_rotation_rad", 0.0))
    surface_type = cfg.get("mirror.surface_type", "Spherical").strip()
    radius = float(cfg.get("mirror.radius_of_curvature", 16.0))

    if fmt == "raw1229":
        pattern = cfg.get("mirror.series_csv_pattern")
        if not pattern:
            raise ValueError("mirror.series_csv_pattern is required when mirror.series_format=raw1229")
        angles = sorted(parse_float_list(cfg["mirror.series_angles_deg"]))
        swap_xy = cfg.get("mirror.series_swap_xy", "true").strip().lower() in {"1", "true", "yes", "on"}

        lower_angle = angles[0]
        upper_angle = angles[-1]
        for angle in angles:
            if angle <= elevation:
                lower_angle = angle
            if angle >= elevation:
                upper_angle = angle
                break
        if elevation <= angles[0]:
            lower_angle = upper_angle = angles[0]
        elif elevation >= angles[-1]:
            lower_angle = upper_angle = angles[-1]

        lower_path = Path(expand_angle_pattern(pattern, lower_angle))
        upper_path = Path(expand_angle_pattern(pattern, upper_angle))
        if not lower_path.is_absolute():
            lower_path = resolve_workspace_path(cfg_path, str(lower_path))
        if not upper_path.is_absolute():
            upper_path = resolve_workspace_path(cfg_path, str(upper_path))

        lower_states = _read_raw1229_states(lower_path, swap_xy=swap_xy)
        upper_states = lower_states if upper_angle == lower_angle else _read_raw1229_states(upper_path, swap_xy=swap_xy)
        if len(lower_states) != len(upper_states):
            raise ValueError("mirror elevation series files have mismatched facet counts")

        t = 0.0 if upper_angle == lower_angle else (elevation - lower_angle) / (upper_angle - lower_angle)

        facets = []
        for idx, (lo, hi) in enumerate(zip(lower_states, upper_states)):
            if (lo["cell_index"], lo["ring_index"]) != (hi["cell_index"], hi["ring_index"]):
                raise ValueError("mirror elevation series files have inconsistent facet ordering")
            facets.append({
                "id": idx,
                "center": (1.0 - t) * lo["center"] + t * hi["center"],
                "normal": _slerp(lo["normal"], hi["normal"], t),
                "shape": shape,
                "size1": size1,
                "size2": size2,
                "rotation": rotation,
                "surface_type": surface_type,
                "radius_of_curvature": radius,
            })
        return facets

    csv_path = cfg.get("mirror.series_csv_path")
    if not csv_path:
        raise ValueError("mirror.series_csv_path is required when mirror.series_format=facet_csv")
    csv_path = resolve_workspace_path(cfg_path, csv_path)
    by_angle = _read_elevation_series_csv(csv_path)
    angles = sorted(by_angle)

    lower_angle = angles[0]
    upper_angle = angles[-1]
    for angle in angles:
        if angle <= elevation:
            lower_angle = angle
        if angle >= elevation:
            upper_angle = angle
            break
    if elevation <= angles[0]:
        lower_angle = upper_angle = angles[0]
    elif elevation >= angles[-1]:
        lower_angle = upper_angle = angles[-1]

    lower_states = by_angle[lower_angle]
    upper_states = by_angle[upper_angle]
    if len(lower_states) != len(upper_states):
        raise ValueError("mirror series CSV has mismatched facet counts between elevations")
    t = 0.0 if upper_angle == lower_angle else (elevation - lower_angle) / (upper_angle - lower_angle)

    facets = []
    for lo, hi in zip(lower_states, upper_states):
        if lo["id"] != hi["id"]:
            raise ValueError("mirror series CSV has inconsistent facet ordering")
        facets.append({
            "id": lo["id"],
            "center": (1.0 - t) * lo["center"] + t * hi["center"],
            "normal": _slerp(lo["normal"], hi["normal"], t),
            "shape": lo.get("shape", shape),
            "size1": lo.get("size1", size1),
            "size2": lo.get("size2", size2),
            "rotation": lo.get("rotation", rotation),
            "surface_type": lo.get("surface_type", surface_type),
            "radius_of_curvature": lo.get("radius_of_curvature", radius),
        })
    return facets


def load_facets_from_config(cfg_path, cfg):
    facets = load_facets_from_scoped_mirror_cfg(cfg_path, cfg)

    deformation_path = cfg.get("error.structural_deformation_config", "").strip()
    if deformation_path.lower() in {"", "none", "off", "false", "no"}:
        return facets

    deformation_cfg_path = Path(deformation_path)
    if not deformation_cfg_path.is_absolute():
        deformation_cfg_path = _resolve_include_path(cfg_path, deformation_path)
    deformation_raw = read_key_value_config(deformation_cfg_path)
    deformation_cfg = {
        _scoped_key(key, "mirror."): value
        for key, value in deformation_raw.items()
    }
    deformation_cfg["telescope.pointing_el_deg"] = cfg.get(
        "telescope.pointing_el_deg",
        deformation_cfg.get("telescope.pointing_el_deg", "90"),
    )
    deformation_cfg["mirror.series_elevation_deg"] = cfg.get(
        "mirror.series_elevation_deg",
        cfg.get("telescope.pointing_el_deg", "90"),
    )

    deformation_facets = load_facets_from_scoped_mirror_cfg(deformation_cfg_path, deformation_cfg)
    if len(facets) != len(deformation_facets):
        raise ValueError("structural deformation facet count does not match base mirror layout")

    merged = []
    for base, deform in zip(facets, deformation_facets):
        item = dict(base)
        item["center"] = deform["center"]
        item["normal"] = deform["normal"]
        merged.append(item)
    return merged


def source_direction_from_config(cfg):
    theta_text = cfg.get("source.beam_theta_deg")
    phi_text = cfg.get("source.beam_phi_deg")
    if theta_text is not None or phi_text is not None:
        theta = math.radians(float(theta_text or 0.0))
        phi = math.radians(float(phi_text or 0.0))
        return np.array([
            math.sin(theta) * math.cos(phi),
            math.sin(theta) * math.sin(phi),
            -math.cos(theta),
        ], dtype=float)
    return _normalize(parse_vec3(cfg.get("source.beam_direction"), [0.0, 0.0, -1.0]))


def telescope_frame_from_config(cfg):
    az = math.radians(float(cfg.get("telescope.pointing_az_deg", 0.0)))
    el = math.radians(float(cfg.get("telescope.pointing_el_deg", 90.0)))
    origin = parse_vec3(cfg.get("telescope.position_m"), [0.0, 0.0, 0.0])
    z_axis = _normalize([
        math.cos(el) * math.cos(az),
        math.cos(el) * math.sin(az),
        math.sin(el),
    ])
    x_axis = _normalize([-math.sin(az), math.cos(az), 0.0])
    y_axis = _normalize(np.cross(z_axis, x_axis))
    return {
        "origin": origin,
        "x_axis": x_axis,
        "y_axis": y_axis,
        "z_axis": z_axis,
    }


def normalize_source_coordinate_frame(frame_name):
    """Mirror ``normalizeSourceCoordinateFrame`` in OpticalSimCommon.cpp."""
    name = str(frame_name or "").strip().lower()
    aliases = {
        "": "telescope_local",
        "local": "telescope_local",
        "optical_local": "telescope_local",
        "corsika_iact": "corsika_nwu_relative",
        "corsika": "corsika_nwu_relative",
        "simtelarray": "corsika_nwu_relative",
        "corsika_global": "corsika_nwu_global",
        "enu_relative": "enu_east_relative",
        "east_north_up_relative": "enu_east_relative",
        "east_start_relative": "enu_east_relative",
        "enu_global": "enu_east_global",
        "east_north_up_global": "enu_east_global",
        "east_start_global": "enu_east_global",
        "generic_global": "lact_generic_global",
        "array_global": "lact_generic_global",
        "global": "lact_generic_global",
    }
    name = aliases.get(name, name)
    allowed = {
        "telescope_local",
        "corsika_nwu_relative",
        "corsika_nwu_global",
        "enu_east_relative",
        "enu_east_global",
        "lact_generic_global",
    }
    if name not in allowed:
        raise ValueError(f"unsupported source.coordinate_frame: {frame_name}")
    return name


def source_coordinate_frame_name_from_config(cfg):
    """Apply the same coordinate-frame default rules as the C++ runtime."""
    if "source.coordinate_frame" in cfg:
        return normalize_source_coordinate_frame(cfg["source.coordinate_frame"])

    mode = str(cfg.get("source.mode", "ParallelBeam")).strip().lower()
    if mode in {"eventio", "corsika", "corsikaeventio", "corsika_eventio", "iact"}:
        return normalize_source_coordinate_frame(
            cfg.get("source.eventio_coordinate_frame", "corsika_iact")
        )
    if mode in {"photoncsv", "photon_csv", "csv", "file"} and "source.local_telescope_frame" in cfg:
        local = str(cfg["source.local_telescope_frame"]).strip().lower() in {
            "1", "true", "yes", "on",
        }
        return "telescope_local" if local else "lact_generic_global"
    return "telescope_local"


def source_telescope_frame_from_config(cfg):
    """Return the exact physical source-adapter basis used by the C++ runtime."""
    name = source_coordinate_frame_name_from_config(cfg)
    origin = parse_vec3(cfg.get("telescope.position_m"), [0.0, 0.0, 0.0])
    az = math.radians(float(cfg.get("telescope.pointing_az_deg", 0.0)))
    el = math.radians(float(cfg.get("telescope.pointing_el_deg", 90.0)))
    sin_az, cos_az = math.sin(az), math.cos(az)
    sin_el, cos_el = math.sin(el), math.cos(el)

    if name in {"corsika_nwu_relative", "corsika_nwu_global"}:
        axes = (
            [-sin_el * cos_az, sin_el * sin_az, cos_el],
            [-sin_az, -cos_az, 0.0],
            [cos_el * cos_az, -cos_el * sin_az, sin_el],
        )
    elif name in {"enu_east_relative", "enu_east_global"}:
        axes = (
            [-sin_el * cos_az, -sin_el * sin_az, cos_el],
            [sin_az, -cos_az, 0.0],
            [cos_el * cos_az, cos_el * sin_az, sin_el],
        )
    elif name == "telescope_local":
        axes = ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0])
    else:
        return telescope_frame_from_config(cfg)

    return {
        "origin": origin,
        "x_axis": _normalize(axes[0]),
        "y_axis": _normalize(axes[1]),
        "z_axis": _normalize(axes[2]),
    }


def rotate_local_vector(vec, frame):
    vec = np.asarray(vec, dtype=float)
    return (
        frame["x_axis"] * vec[0] +
        frame["y_axis"] * vec[1] +
        frame["z_axis"] * vec[2]
    )


def point_to_global(point, frame):
    return frame["origin"] + rotate_local_vector(point, frame)


def apply_telescope_frame_to_facets(facets, frame):
    transformed = []
    for facet in facets:
        item = dict(facet)
        local_u, local_v, _ = local_frame_for_normal(facet["normal"])
        item["center"] = point_to_global(facet["center"], frame)
        item["normal"] = _normalize(rotate_local_vector(facet["normal"], frame))
        item["u_axis"] = _normalize(rotate_local_vector(local_u, frame))
        item["v_axis"] = _normalize(rotate_local_vector(local_v, frame))
        transformed.append(item)
    return transformed


def local_frame_for_normal(normal):
    n = _normalize(normal)
    ref = np.array([0.0, 0.0, 1.0]) if abs(n[2]) < 0.9 else np.array([0.0, 1.0, 0.0])
    u = _normalize(np.cross(ref, n))
    v = _normalize(np.cross(n, u))
    return u, v, n
