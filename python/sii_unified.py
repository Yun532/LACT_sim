#!/usr/bin/env python3
"""恒星强度干涉的最小可复用模拟组件。

模块明确分成三个时间尺度：光学相干与源可见度、用于校准响应的短波形、
以及用显式误差模拟小时级观测的 ``|V|^2`` 统计量。没有必要把整夜的
高采样率波形全部放入内存；短波形负责验证 SPE、时延、ADC 和 SiPM，
长曝光则直接使用相应的充分统计量。
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, replace
from itertools import combinations
from pathlib import Path
import json
import hashlib
import math

import numpy as np
import pandas as pd

try:
    from scipy.special import j1
except ImportError:  # pragma: no cover
    j1 = None


C_M_S = 299_792_458.0
H_J_S = 6.626_070_15e-34
JY_W_M2_HZ = 1.0e-26
MAS_TO_RAD = math.pi / (180.0 * 3600.0 * 1000.0)
# 固定赤道坐标场景采用平均恒星日；曝光仍以SI秒计，不把24小时当作恒星日。
SIDEREAL_DAY_S = 86164.0905


@dataclass(frozen=True)
class Instrument:
    """单像素仪器参数；可由 ``from_repository`` 跟随 main 配置更新。"""

    wavelength_nm: float = 400.0
    optical_width_nm: float = 2.0
    effective_area_m2: float = 24.576860
    throughput: float = 0.20
    electronics_bandwidth_hz: float = 312.5e6
    adc_sample_rate_hz: float = 625.0e6
    # main 没有提供 ADC 位数和满量程时，0 表示保留模拟波形、不量化。
    adc_bits: int = 0
    adc_full_scale_mv: float = 0.0
    detected_nsb_rate_hz: float = 0.573135e6
    # 以下未在 main 中标定的随机效应默认关闭；有实测值后直接替换。
    electronic_noise_rms_mv: float = 0.0
    excess_noise_factor: float = 1.016142
    optical_timing_in_effective_bandwidth: bool = False
    polarization_factor: float = 0.5
    spectral_shape_factor: float = 1.0
    microcells_per_pixel: int = 270_336
    microcell_recovery_time_ns: float = 0.0
    intrinsic_time_jitter_ns: float = 0.0
    telescope_gain_calibration_rms: float = 0.0
    per_night_gain_rms: float = 0.0
    baseline_zero_point_rms: float = 0.0
    residual_timing_rms_ns: float = 0.0
    transparency_fractional_rms: float = 0.0
    nsb_fractional_rms: float = 0.0
    # main 当前没有这些单像素非理想项；接口保留，默认严格为零。
    sipm_crosstalk_probability: float = 0.0
    sipm_afterpulse_probability: float = 0.0
    dark_count_rate_hz: float = 0.0
    filter_angular_response_path: str | None = None
    sii_bandpass_path: str | None = None
    spe_template_path: str | None = None
    charge_samples_path: str | None = None
    microcell_device_path: str | None = None
    optical_timing_kernel_path: str | None = None
    # 窄带积分节点及HBT功率权重；空元组表示显式单色场景。
    visibility_wavelength_nm: tuple = ()
    visibility_spectral_weights: tuple = ()
    parameter_source: str = "built-in defaults"

    @property
    def wavelength_m(self) -> float:
        return self.wavelength_nm * 1.0e-9

    @property
    def optical_width_m(self) -> float:
        return self.optical_width_nm * 1.0e-9

    @property
    def optical_bandwidth_hz(self) -> float:
        return C_M_S * self.optical_width_m / self.wavelength_m**2

    @property
    def sample_width_ns(self) -> float:
        return 1.0e9 / self.adc_sample_rate_hz

    @property
    def coherence_area_s(self) -> float:
        return (self.polarization_factor * self.spectral_shape_factor
                * self.wavelength_m**2
                / (C_M_S * self.optical_width_m))

    @classmethod
    def from_repository(cls, repo_root, response_config=None,
                        source_transmission_scale=0.7836336, **overrides):
        """从 main 使用的 cfg/CSV 读取探测器响应，显式参数仍可覆盖。

        ``source_transmission_scale`` 表示大气及尚未单列的通光损失。当前
        数值是固定场景假设，不是实测大气透过率；镜面、滤光片、PDE、NSB、
        SPE、采样间隔或微单元配置在 main 中更新后会在下次调用时重新读取。
        """
        root = Path(repo_root).resolve()
        if response_config is None:
            response_config = (
                root / "configs" / "examples"
                / "corsika_lact_pylast_root_only_measured_waveform.cfg")
        config_path = Path(response_config).resolve()
        from config_io import (expand_component_config,
                               read_key_value_config,
                               resolve_workspace_path)

        cfg, component_paths = expand_component_config(config_path)
        base = cls()
        wavelength_nm = float(overrides.get("wavelength_nm", base.wavelength_nm))
        width_nm = float(overrides.get("optical_width_nm", base.optical_width_nm))
        if not (np.isfinite(wavelength_nm) and np.isfinite(width_nm)
                and 0 < width_nm < 2*wavelength_nm):
            raise ValueError("SII wavelength and passband width must be positive and finite")

        def workspace_path(value):
            return resolve_workspace_path(config_path, value)

        mirror_path = workspace_path(cfg["efficiency.mirror_reflectivity"])
        filter_path = workspace_path(cfg["efficiency.filter_transmission"])
        pde_value = cfg["sipm.pde"]
        pde_path = workspace_path(pde_value)
        nsb_spectrum_path = Path(cfg["nsb.spectrum_csv"])
        effective_area = float(cfg["nsb.effective_area_m2"])
        collector = float(cfg["nsb.collector_mean_transmission"])
        focal_length = float(cfg.get("telescope.focal_length_m", 8.0))

        camera_csv = workspace_path(cfg["camera.csv_path"])
        camera_data = np.genfromtxt(camera_csv, delimiter=",", names=True)
        pixel_size = float(np.atleast_1d(camera_data["size_m"])[0])

        wave_mirror, mirror = _read_two_column_curve(mirror_path)
        wave_filter, filt = _read_two_column_curve(filter_path)
        wave_pde, pde = _read_two_column_curve(pde_path)
        throughput = (
            np.interp(wavelength_nm, wave_mirror, mirror)
            * np.interp(wavelength_nm, wave_filter, filt)
            * np.interp(wavelength_nm, wave_pde, pde)
            * collector * float(source_transmission_scale))

        wave_nsb, nsb_flux = _read_two_column_curve(nsb_spectrum_path)
        lower, upper = wavelength_nm-width_nm/2, wavelength_nm+width_nm/2
        bandpass_path = overrides.get("sii_bandpass_path")
        bandpass = (None if bandpass_path is None else
                    _read_two_column_curve(workspace_path(bandpass_path)))
        knots = [wave_mirror, wave_filter, wave_pde, wave_nsb,
                 np.linspace(lower, upper, 1025)]
        if bandpass is not None:
            knots.append(bandpass[0])
        wavelength = np.unique(np.concatenate(knots))
        wavelength = wavelength[(wavelength >= lower) & (wavelength <= upper)]
        integration_weights = None
        if bandpass is not None:
            # 在实测曲线的每个线性区间积分，窄峰不会被固定全带宽节点漏掉。
            edges = np.unique(np.concatenate((wave_mirror, wave_filter,
                wave_pde, wave_nsb, bandpass[0], [lower, upper])))
            edges = edges[(edges >= lower) & (edges <= upper)]
            nodes, weights = np.polynomial.legendre.leggauss(5)
            half = np.diff(edges)/2
            wavelength = ((edges[:-1]+half)[:, None]
                          + half[:, None]*nodes).ravel()
            integration_weights = (half[:, None]*weights).ravel()
        def integrate(values):
            return (np.trapezoid(values, wavelength) if integration_weights is None
                    else np.dot(integration_weights, values))
        response = np.ones_like(wavelength)
        for x, y in ((wave_mirror, mirror), (wave_filter, filt), (wave_pde, pde)):
            response *= np.interp(wavelength, x, y, left=0., right=0.)
        if bandpass is not None:
            if np.any((bandpass[1] < 0) | (bandpass[1] > 1)):
                raise ValueError("SII passband transmission must lie in [0, 1]")
            response *= np.interp(wavelength, *bandpass, left=0., right=0.)
        photon_integral = integrate(response/wavelength)
        if photon_integral <= 0:
            raise ValueError("SII passband has no detected throughput")
        central_response = (np.interp(wavelength_nm, wave_mirror, mirror)
                            * np.interp(wavelength_nm, wave_filter, filt)
                            * np.interp(wavelength_nm, wave_pde, pde))
        if central_response <= 0:
            raise ValueError("central wavelength has no detected throughput")
        band_average_ratio = photon_integral/(width_nm/wavelength_nm)/central_response
        spectral_shape = (width_nm*integrate(response**2)
                          / (wavelength_nm**2*photon_integral**2))
        # HBT相关面积按探测光谱密度平方加权；平坦AB谱下换成波长后正比于响应平方。
        nodes, quadrature = np.polynomial.legendre.leggauss(5)
        spectral_nodes = wavelength_nm+width_nm/2*nodes
        spectral_weights = quadrature*np.interp(spectral_nodes, wavelength, response)**2
        if integration_weights is not None:
            spectral_nodes = wavelength.copy()
            spectral_weights = integration_weights*response**2
            active = spectral_weights > 0
            spectral_nodes, spectral_weights = spectral_nodes[active], spectral_weights[active]
        if not np.isfinite(spectral_weights.sum()) or spectral_weights.sum() <= 0:
            raise ValueError("SII visibility quadrature has no finite positive response")
        spectral_weights /= spectral_weights.sum()
        detected_spectrum = np.interp(wavelength, wave_nsb, nsb_flux,
                                     left=0., right=0.)*response
        # main 的诊断工具把结果写成 p.e./ns；这里的公共接口统一使用 Hz。
        nsb_rate = (integrate(detected_spectrum)
                    * effective_area * (pixel_size/focal_length)**2 * collector)

        electronics_path = component_paths["electronics"]
        electronics = read_key_value_config(electronics_path)
        electronics_dir = electronics_path.parent
        template_path = electronics_dir / electronics["single_pe.template"]
        charge_path = electronics_dir / electronics[
            "single_pe.charge_fluctuation.samples"]
        charge = load_empirical_charge_factors(charge_path)
        excess_noise = float(np.sqrt(np.mean(charge**2)))
        sample_width_ns = float(electronics.get("sampling.width_ns", 4.0))

        device_path = electronics_dir / electronics["microcell.device"]
        device = read_key_value_config(device_path)
        microcells = math.prod(int(device[key]) for key in (
            "channel_columns", "channel_rows",
            "microcell_columns_per_channel", "microcell_rows_per_channel"))
        recovery_ns = (float(electronics.get("microcell.recovery_time_ns", 0.0))
                       if electronics.get("microcell.recovery_enabled", "false").lower() == "true"
                       else 0.0)
        # 优先使用由当前 main 实测光学配置生成的完整响应；没有运行过该
        # 校准时才回退到旧的轴上理想误差时间核，并由 notebook 明确标注。
        full_timing_path = (root / "configs" / "optics"
                            / "lact2_measured_single_pixel_400nm.csv")
        fallback_timing_path = (root / "configs" / "optics"
                                / "lact_1229_onaxis_timing_kernel.csv")
        timing_path = (full_timing_path if (full_timing_path.exists()
                                             and np.isclose(wavelength_nm, 400.0))
                       else fallback_timing_path)
        full_provenance_path = full_timing_path.with_suffix(".provenance.json")
        if (full_timing_path.exists() and full_provenance_path.exists()
                and np.isclose(wavelength_nm, 400.0)):
            full_response = json.loads(
                full_provenance_path.read_text(encoding="utf-8"))
            central_effective_area = float(
                full_response["central_pixel_effective_detection_area_m2"])
            # 光追已包含镜面、遮挡、PSF、相机、集光器、滤光片和PDE；
            # source_transmission_scale 只继续表示大气及未单列的上游损失。
            throughput = (central_effective_area/effective_area
                          * float(source_transmission_scale))

        values = {
            "effective_area_m2": effective_area,
            "throughput": float(throughput*band_average_ratio),
            "spectral_shape_factor": float(spectral_shape),
            "visibility_wavelength_nm": tuple(spectral_nodes),
            "visibility_spectral_weights": tuple(spectral_weights),
            "adc_sample_rate_hz": 1.0e9/sample_width_ns,
            # main 当前只给采样间隔、没有独立标定的模拟前端带宽；相关器
            # 可用数字频带以fs/2标记，实际响应由SPE决定；这不是实测前端带宽。
            "electronics_bandwidth_hz": 0.5e9/sample_width_ns,
            "adc_bits": int(electronics.get("adc.bits", base.adc_bits)),
            "adc_full_scale_mv": float(electronics.get(
                "adc.full_scale_mv", base.adc_full_scale_mv)),
            "electronic_noise_rms_mv": float(electronics.get(
                "electronics.noise_rms_mv", base.electronic_noise_rms_mv)),
            "sipm_crosstalk_probability": float(electronics.get(
                "sipm.crosstalk_probability",
                base.sipm_crosstalk_probability)),
            "sipm_afterpulse_probability": float(electronics.get(
                "sipm.afterpulse_probability",
                base.sipm_afterpulse_probability)),
            "dark_count_rate_hz": float(electronics.get(
                "sipm.dark_count_rate_hz", base.dark_count_rate_hz)),
            "detected_nsb_rate_hz": float(nsb_rate),
            "excess_noise_factor": excess_noise,
            "microcells_per_pixel": microcells,
            "microcell_recovery_time_ns": recovery_ns,
            "intrinsic_time_jitter_ns": (float(electronics.get("single_pe.time_jitter.sigma_ns", 0.0))
                if electronics.get("single_pe.time_jitter.enabled", "false").lower() == "true" else 0.0),
            "spe_template_path": str(template_path.resolve()),
            "charge_samples_path": str(charge_path.resolve()),
            "microcell_device_path": str(device_path.resolve()),
            "optical_timing_kernel_path": (
                str(timing_path.resolve()) if timing_path.exists() else None),
            "filter_angular_response_path": (
                str(workspace_path(cfg["efficiency.filter_angular_response"]).resolve())
                if cfg.get("efficiency.filter_angular_response") else None),
            "parameter_source": str(config_path),
        }
        values.update(overrides)
        if bandpass_path is not None:
            values["sii_bandpass_path"] = str(workspace_path(bandpass_path).resolve())
        return replace(base, **values)


@dataclass(frozen=True)
class BinarySource:
    """双星亮度、角分离、位置角、流量比和两颗星的角直径。"""

    ab_magnitude: float = 2.0
    separation_mas: float = 0.20
    position_angle_deg: float = 35.0
    flux_ratio_secondary_to_primary: float = 0.55
    primary_diameter_mas: float = 0.060
    secondary_diameter_mas: float = 0.040


@dataclass(frozen=True)
class EllipseSource:
    """均匀椭圆源；位置角从北向东计算。"""

    ab_magnitude: float = 2.0
    major_diameter_mas: float = 0.24
    minor_diameter_mas: float = 0.12
    position_angle_deg: float = 35.0


@dataclass(frozen=True)
class TransitSource:
    """均匀恒星盘减去完全位于盘内的不透明行星圆盘。"""

    ab_magnitude: float = 2.0
    stellar_diameter_mas: float = 0.30
    planet_diameter_mas: float = 0.08
    planet_east_offset_mas: float = 0.06
    planet_north_offset_mas: float = 0.03


@dataclass(frozen=True)
class Observation:
    """一次完整观测的几何和积分设置。"""

    site_lat_deg: float = 29.36
    source_dec_deg: float = 22.0
    hours_per_night: float = 6.0
    nights: int = 1
    segment_s: float = 1200.0
    visibility_subsamples_per_segment: int = 9


@dataclass
class PipelineResult:
    """一键流程的输出；重建可关闭，因此允许为 ``None``。"""

    uvw: pd.DataFrame
    measurements: pd.DataFrame
    reconstruction: object | None
    metadata: dict


@dataclass
class WaveformGLSCalibration:
    """由短波形标定出的完整相关峰和单块噪声协方差。"""

    lags_ns: np.ndarray
    null_mean: np.ndarray
    peak_per_visibility2: np.ndarray
    covariance_per_block: np.ndarray
    block_duration_s: float
    star_rate_hz: float
    background_rate_hz: float
    calibration_visibility2: float
    hbt_pair_rate_scale: float
    null_records: int
    signal_records: int
    covariance_shrinkage: float
    instrument_signature: str | None = None
    response_relative_uncertainty: float = 0.0
    sigma_relative_uncertainty: float = 0.0
    phase_model: str = ""

    def covariance(self, exposure_s: float) -> np.ndarray:
        """按独立数据块平均的 ``1/T`` 定律返回指定曝光协方差。"""
        if exposure_s <= 0.0:
            raise ValueError("exposure_s must be positive")
        return self.covariance_per_block * self.block_duration_s / exposure_s


def _read_two_column_curve(path) -> tuple[np.ndarray, np.ndarray]:
    """读取 main 中通用的两列光谱 CSV。"""
    data = np.genfromtxt(path, delimiter=",", names=True)
    names = list(data.dtype.names or ())
    if len(names) < 2:
        raise ValueError(f"{path} needs two numeric columns")
    x = np.atleast_1d(np.asarray(data[names[0]], float))
    y = np.atleast_1d(np.asarray(data[names[1]], float))
    finite = np.isfinite(x) & np.isfinite(y)
    order = np.argsort(x[finite])
    return x[finite][order], y[finite][order]


def source_direction_enu(hour_angle_rad: float, dec_rad: float,
                         lat_rad: float) -> np.ndarray:
    """返回 ENU 坐标中的源方向单位向量。

    时角H=LST-RA，西向为正；后面的天球切向基底与此约定一致。
    """
    vector = np.array([
        -math.cos(dec_rad) * math.sin(hour_angle_rad),
        math.sin(dec_rad) * math.cos(lat_rad)
        - math.cos(dec_rad) * math.cos(hour_angle_rad) * math.sin(lat_rad),
        math.sin(dec_rad) * math.sin(lat_rad)
        + math.cos(dec_rad) * math.cos(hour_angle_rad) * math.cos(lat_rad),
    ])
    return vector / np.linalg.norm(vector)


def celestial_tangent_axes_enu(hour_angle_rad: float, dec_rad: float,
                                lat_rad: float) -> tuple[np.ndarray, ...]:
    """返回固定 RA/Dec 天球基底在 ENU 中的三个单位轴。

    ``u`` 是当前时角符号下的天球东西向切向轴，``v`` 是赤纬增加方向，
    ``w`` 指向源。它们由源方向对时角和赤纬的导数得到，不会像
    ``local_up × source`` 那样随视差角额外旋转。
    """
    hour = float(hour_angle_rad)
    dec = float(dec_rad)
    lat = float(lat_rad)
    sh, ch = math.sin(hour), math.cos(hour)
    sd, cd = math.sin(dec), math.cos(dec)
    sl, cl = math.sin(lat), math.cos(lat)
    source = source_direction_enu(hour, dec, lat)
    u_axis = np.array([ch, -sl * sh, cl * sh], dtype=float)
    v_axis = np.array([
        sd * sh,
        cd * cl + sd * ch * sl,
        cd * sl - sd * ch * cl,
    ], dtype=float)
    # 显式归一化抵抗浮点误差，并让 w 与源方向完全一致。
    u_axis /= np.linalg.norm(u_axis)
    v_axis /= np.linalg.norm(v_axis)
    return u_axis, v_axis, source


def uvw_from_enu(baseline_enu_m, hour_angle_rad: float, dec_rad: float,
                 lat_rad: float) -> np.ndarray:
    """把 ENU 基线投影到固定 RA/Dec 天球的 ``(u,v,w)`` 坐标。

    ``w`` 是几何传播方向，几何时延为 ``w/c``；交换基线端点会严格得到
    ``(-u,-v,-w)``。该定义与双星位置角使用的 RA/Dec 切平面保持一致。
    """
    u_axis, v_axis, w_axis = celestial_tangent_axes_enu(
        hour_angle_rad, dec_rad, lat_rad)
    baseline = np.asarray(baseline_enu_m, dtype=float)
    return np.array([baseline @ u_axis, baseline @ v_axis,
                     baseline @ w_axis])


def uniform_disk_visibility(q_lambda, diameter_mas: float) -> np.ndarray:
    """计算均匀圆盘的实复可见度。"""
    if j1 is None:
        raise RuntimeError("SciPy is required for uniform-disk visibility")
    x = np.pi * diameter_mas * MAS_TO_RAD * np.asarray(q_lambda, float)
    output = np.ones_like(x)
    nonzero = np.abs(x) > 1.0e-12
    output[nonzero] = 2.0 * j1(x[nonzero]) / x[nonzero]
    return output


def binary_visibility(u_lambda, v_lambda,
                      source: BinarySource = BinarySource()) -> np.ndarray:
    """计算两颗均匀圆盘双星的归一化复可见度。"""
    u = np.asarray(u_lambda, dtype=float)
    v = np.asarray(v_lambda, dtype=float)
    q = np.hypot(u, v)
    pa = math.radians(source.position_angle_deg)
    dx = 0.5 * source.separation_mas * math.sin(pa) * MAS_TO_RAD
    dy = 0.5 * source.separation_mas * math.cos(pa) * MAS_TO_RAD
    first = (uniform_disk_visibility(q, source.primary_diameter_mas)
             * np.exp(-2j * np.pi * (-u * dx - v * dy)))
    second = (uniform_disk_visibility(q, source.secondary_diameter_mas)
              * np.exp(-2j * np.pi * (u * dx + v * dy)))
    ratio = source.flux_ratio_secondary_to_primary
    return (first + ratio * second) / (1.0 + ratio)


def ellipse_visibility(u_lambda, v_lambda,
                       source: EllipseSource = EllipseSource()) -> np.ndarray:
    """均匀椭圆的归一化复可见度（椭圆中心固定在相位中心）。"""
    u = np.asarray(u_lambda, dtype=float)
    v = np.asarray(v_lambda, dtype=float)
    pa = math.radians(source.position_angle_deg)
    along_major = u*math.sin(pa) + v*math.cos(pa)
    along_minor = u*math.cos(pa) - v*math.sin(pa)
    x = np.pi*MAS_TO_RAD*np.hypot(
        source.major_diameter_mas*along_major,
        source.minor_diameter_mas*along_minor)
    output = np.ones_like(x)
    nonzero = np.abs(x) > 1.0e-12
    output[nonzero] = 2.0*j1(x[nonzero])/x[nonzero]
    return output.astype(complex)


def transit_visibility(u_lambda, v_lambda,
                       source: TransitSource = TransitSource()) -> np.ndarray:
    """均匀恒星盘被盘内不透明圆行星遮挡时的归一化复可见度。"""
    star_radius = 0.5*source.stellar_diameter_mas
    planet_radius = 0.5*source.planet_diameter_mas
    offset = math.hypot(
        source.planet_east_offset_mas, source.planet_north_offset_mas)
    if planet_radius <= 0.0 or planet_radius+offset > star_radius:
        raise ValueError("transit planet must lie completely inside the stellar disk")
    u = np.asarray(u_lambda, dtype=float)
    v = np.asarray(v_lambda, dtype=float)
    q = np.hypot(u, v)
    blocked_fraction = (source.planet_diameter_mas
                        / source.stellar_diameter_mas)**2
    phase = np.exp(-2j*np.pi*MAS_TO_RAD*(
        u*source.planet_east_offset_mas
        + v*source.planet_north_offset_mas))
    star = uniform_disk_visibility(q, source.stellar_diameter_mas)
    planet = uniform_disk_visibility(q, source.planet_diameter_mas)*phase
    return (star-blocked_fraction*planet)/(1.0-blocked_fraction)


def source_visibility(u_lambda, v_lambda, source,
                      source_case: str = "binary") -> np.ndarray:
    """把所有内置源模型映射为同一个归一化复可见度接口。"""
    if source_case == "binary":
        if not isinstance(source, BinarySource):
            raise TypeError("binary source_case needs BinarySource")
        return binary_visibility(u_lambda, v_lambda, source)
    if source_case == "single_disk":
        if not isinstance(source, BinarySource):
            raise TypeError("single_disk source_case needs BinarySource")
        return uniform_disk_visibility(
            np.hypot(u_lambda, v_lambda), source.primary_diameter_mas)
    if source_case == "ellipse":
        if not isinstance(source, EllipseSource):
            raise TypeError("ellipse source_case needs EllipseSource")
        return ellipse_visibility(u_lambda, v_lambda, source)
    if source_case == "transit":
        if not isinstance(source, TransitSource):
            raise TypeError("transit source_case needs TransitSource")
        return transit_visibility(u_lambda, v_lambda, source)
    raise ValueError("source_case must be binary, single_disk, ellipse or transit")


def ab_photon_spectral_density(magnitude: float, wavelength_m: float) -> float:
    """把 AB 星等换算为每平方米、每赫兹的光子谱密度。"""
    f_nu = 3631.0 * JY_W_M2_HZ * 10.0 ** (-0.4 * magnitude)
    photon_energy = H_J_S * C_M_S / wavelength_m
    return f_nu / photon_energy


def detected_star_rate_hz(magnitude: float,
                          instrument: Instrument = Instrument()) -> float:
    """由 AB 星等、有效面积、通光效率和光学带宽得到探测光子率。"""
    return (instrument.effective_area_m2 * instrument.throughput
            * ab_photon_spectral_density(magnitude, instrument.wavelength_m)
            * instrument.optical_bandwidth_hz)


def unit_visibility_snr(magnitude: float, integration_s: float,
                        instrument: Instrument = Instrument(),
                        spectral_channels: int = 1,
                        nsb_rate_hz: float | None = None) -> float:
    """计算未分辨等口径双镜的散粒噪声近似 SNR。

    使用和短波形相关对完全相同的有效相干面积 ``coherence_area_s``；
    其中已经包含偏振与有限滤光片谱形造成的相关对比度稀释。
    """
    nsb = (instrument.detected_nsb_rate_hz if nsb_rate_hz is None
           else nsb_rate_hz)
    star = detected_star_rate_hz(magnitude, instrument)
    total = star + nsb + instrument.dark_count_rate_hz
    one_channel = (
        star**2 / total * instrument.coherence_area_s
        # B_eff 是 0..Nyquist 的单边等效带宽。实值波形在时长 T 内
        # 有 2*B_eff*T 个独立实自由度；白噪声极限 B=fs/2 时正好
        # 回到 N=fs*T 个独立 ADC 样点。
        * math.sqrt(2.0 * instrument.electronics_bandwidth_hz * integration_s)
        / instrument.excess_noise_factor**2
    )
    return one_channel * math.sqrt(spectral_channels)


def mean_recovery_fraction(rate_hz: float, microcells: int,
                           recovery_time_ns: float) -> float:
    """计算均匀随机照明下大量微单元的平均电荷比例。

    每个微单元看到的 Poisson 率为 ``rate_hz/microcells``。对指数间隔分布平均
    ``1-exp(-dt/tau)`` 后，结果正好是 ``1/(1 + rate_per_cell*tau)``。
    """
    if rate_hz < 0 or microcells <= 0 or recovery_time_ns < 0:
        raise ValueError("rate, microcell count, and recovery time are invalid")
    occupancy = rate_hz * recovery_time_ns * 1.0e-9 / microcells
    return 1.0 / (1.0 + occupancy)


def apply_exponential_microcell_recovery(times_ns, cell_ids,
                                         recovery_time_ns: float) -> np.ndarray:
    """对按时间排序的微单元击中序列计算指数恢复电荷比例。"""
    times = np.asarray(times_ns, dtype=float)
    cells = np.asarray(cell_ids, dtype=np.int64)
    if times.shape != cells.shape or recovery_time_ns < 0:
        raise ValueError("times/cells or recovery time are invalid")
    # 恢复时间尚未实测时取0：所有光电子权重为1，不把SPE波形尾部当恢复时间。
    if recovery_time_ns == 0:
        return np.ones_like(times)
    order = np.argsort(times, kind="stable")
    fractions = np.empty_like(times)
    last: dict[int, float] = {}
    for index in order:
        cell = int(cells[index])
        previous = last.get(cell)
        if previous is None:
            fraction = 1.0
        else:
            dt = max(0.0, float(times[index] - previous))
            fraction = -math.expm1(-dt / recovery_time_ns)
        fractions[index] = fraction
        if fraction > 0.0:
            last[cell] = float(times[index])
    return fractions


def load_measured_spe_template(path) -> tuple[np.ndarray, np.ndarray]:
    """读取仓库中的两列实测 SPE 时间/幅度 CSV。"""
    data = np.genfromtxt(Path(path), delimiter=",", names=True)
    names = list(data.dtype.names or ())
    if len(names) < 2:
        raise ValueError(f"{path} does not contain two numeric columns")
    time_ns = np.asarray(data[names[0]], dtype=float)
    amplitude_mv = np.asarray(data[names[1]], dtype=float)
    finite = np.isfinite(time_ns) & np.isfinite(amplitude_mv)
    return time_ns[finite], amplitude_mv[finite]


def load_empirical_charge_factors(
        path, column: str = "charge_factor_mean_one") -> np.ndarray:
    """读取正值实测 SPE 电荷，并严格归一化到均值为一。"""
    data = np.genfromtxt(Path(path), delimiter=",", names=True, dtype=None,
                         encoding="utf-8")
    names = list(data.dtype.names or ())
    if column not in names:
        raise ValueError(f"{path} is missing charge column {column}")
    values = np.asarray(data[column], dtype=float)
    values = values[np.isfinite(values) & (values > 0.0)]
    if values.size == 0:
        raise ValueError(f"{path} has no positive charge factors")
    return values / np.mean(values)


def load_optical_timing_mixture(path) -> dict[str, np.ndarray | float]:
    """读取按镜面分面汇总的光学到达时间混合分布。

    文件中的均值包含任意公共传播时间；抽样和传递计算会去掉这个公共时间，
    因为强度相关峰只受路径展宽影响。
    """
    data = np.genfromtxt(Path(path), delimiter=",", names=True)
    required = {"weight", "mean_time_ns", "std_time_ns"}
    names = set(data.dtype.names or ())
    if not required.issubset(names):
        raise ValueError(f"{path} is missing timing-mixture columns {required-names}")
    weights = np.atleast_1d(np.asarray(data["weight"], dtype=float))
    means = np.atleast_1d(np.asarray(data["mean_time_ns"], dtype=float))
    sigmas = np.atleast_1d(np.asarray(data["std_time_ns"], dtype=float))
    finite = (np.isfinite(weights) & np.isfinite(means) & np.isfinite(sigmas)
              & (weights > 0.0) & (sigmas >= 0.0))
    weights, means, sigmas = weights[finite], means[finite], sigmas[finite]
    if weights.size == 0:
        raise ValueError(f"{path} contains no valid timing components")
    weights = weights / weights.sum()
    global_mean = float(weights @ means)
    centered_means = means - global_mean
    global_variance = float(weights @ (sigmas**2 + centered_means**2))
    return {
        "weights": weights,
        "mean_delay_ns": centered_means,
        "std_delay_ns": sigmas,
        "absolute_mean_time_ns": global_mean,
        "rms_spread_ns": math.sqrt(max(0.0, global_variance)),
    }


def matched_effective_bandwidth_hz(
        instrument: Instrument, source_ab_magnitude: float = 2.0,
        fine_dt_ns: float = 0.01) -> float:
    """由光学时间核、实测 SPE、采样和噪声计算匹配相关带宽。

    ``Instrument.from_repository`` 中的带宽只是 ADC Nyquist 上限；长曝光
    SNR 应使用本函数返回的等效带宽。SPE 在线性散粒噪声极限下同时滤波
    信号和光子噪声，因而相消；存在附加电子噪声或量化噪声时会通过散粒
    噪声占比进入积分。
    """
    if fine_dt_ns <= 0.0:
        raise ValueError("fine_dt_ns must be positive")
    if not instrument.spe_template_path:
        raise ValueError("instrument has no measured SPE template")
    if not instrument.optical_timing_kernel_path:
        raise ValueError("instrument has no optical timing kernel")

    spe_t_ns, spe_mv = load_measured_spe_template(
        instrument.spe_template_path)
    fine_time_ns = np.arange(
        spe_t_ns.min(), spe_t_ns.max() + fine_dt_ns/2.0, fine_dt_ns)
    fine_pulse = np.interp(fine_time_ns, spe_t_ns, spe_mv)
    fft_length = 1 << int(np.ceil(np.log2(len(fine_pulse)*16)))
    full_frequency_hz = np.fft.rfftfreq(fft_length, fine_dt_ns*1.0e-9)
    full_pulse_spectrum = (
        np.fft.rfft(fine_pulse, n=fft_length)*fine_dt_ns*1.0e-9)
    frequency_limit_hz = min(
        instrument.electronics_bandwidth_hz,
        instrument.adc_sample_rate_hz/2.0)
    keep = full_frequency_hz <= frequency_limit_hz
    frequency_hz = full_frequency_hz[keep]
    sampled_pulse = (
        full_pulse_spectrum[keep]
        * np.sinc(frequency_hz/instrument.adc_sample_rate_hz))

    timing = load_optical_timing_mixture(
        instrument.optical_timing_kernel_path)
    weights = np.asarray(timing["weights"], dtype=float)
    means_ns = np.asarray(timing["mean_delay_ns"], dtype=float)
    sigmas_ns = np.asarray(timing["std_delay_ns"], dtype=float)
    frequency_per_ns = frequency_hz*1.0e-9
    optical_transfer = (
        np.exp(-2j*np.pi*frequency_per_ns[:, None]*means_ns[None, :])
        * np.exp(-2*np.pi**2*frequency_per_ns[:, None]**2
                 * sigmas_ns[None, :]**2)) @ weights
    optical_transfer *= np.exp(-2*np.pi**2*frequency_per_ns**2
                               * instrument.intrinsic_time_jitter_ns**2)

    star_rate_hz = detected_star_rate_hz(source_ab_magnitude, instrument)
    total_rate_hz = (star_rate_hz + instrument.detected_nsb_rate_hz
                     + instrument.dark_count_rate_hz)
    shot_psd = (2.0*total_rate_hz*instrument.excess_noise_factor**2
                * np.abs(sampled_pulse)**2)
    adc_step_mv = (
        0.0 if instrument.adc_bits <= 0
        else instrument.adc_full_scale_mv/2**instrument.adc_bits)
    sample_noise_variance_mv2 = (
        instrument.electronic_noise_rms_mv**2 + adc_step_mv**2/12.0)
    additive_noise_psd = (
        2.0*sample_noise_variance_mv2/instrument.adc_sample_rate_hz)
    shot_fraction = np.divide(
        shot_psd, shot_psd + additive_noise_psd,
        out=np.zeros_like(shot_psd),
        where=(shot_psd + additive_noise_psd) > 0.0)
    matched_weight = np.abs(optical_transfer)**4*shot_fraction**2
    return float(np.trapezoid(matched_weight, frequency_hz))


def with_matched_effective_bandwidth(
        instrument: Instrument, source_ab_magnitude: float = 2.0,
        fine_dt_ns: float = 0.01) -> Instrument:
    """返回已把光学时间核计入相关带宽的仪器参数。"""
    bandwidth_hz = matched_effective_bandwidth_hz(
        instrument, source_ab_magnitude, fine_dt_ns)
    return replace(
        instrument,
        electronics_bandwidth_hz=bandwidth_hz,
        optical_timing_in_effective_bandwidth=True)


def sample_optical_delays_ns(rng: np.random.Generator, count: int,
                             mixture) -> np.ndarray:
    """按混合分布抽取已居中的逐光电子光路延迟。"""
    if count < 0:
        raise ValueError("count must be non-negative")
    weights = np.asarray(mixture["weights"], dtype=float)
    means = np.asarray(mixture["mean_delay_ns"], dtype=float)
    sigmas = np.asarray(mixture["std_delay_ns"], dtype=float)
    components = rng.choice(len(weights), size=count, p=weights)
    return rng.normal(means[components], sigmas[components])


def optical_timing_transfer_efficiency(
        mixture, bandwidth_hz: float, samples: int = 4097) -> float:
    """计算光学到达时间展宽在电子带宽内保留的零延迟 HBT 对比度。

    对相同望远镜，若时间特征函数为 H(f)，矩形电子带宽内保留的零延迟
    交叉谱是 ``mean(|H(f)|^2, 0..B)``。校正回原始可见度平方时要除以该值，
    因而其倒数会放大不确定度。
    """
    if bandwidth_hz <= 0.0 or samples < 3:
        raise ValueError("bandwidth and sample count must be positive")
    frequency_hz = np.linspace(0.0, bandwidth_hz, samples)
    frequency_per_ns = frequency_hz * 1.0e-9
    weights = np.asarray(mixture["weights"], dtype=float)
    means = np.asarray(mixture["mean_delay_ns"], dtype=float)
    sigmas = np.asarray(mixture["std_delay_ns"], dtype=float)
    phase = np.exp(-2j * np.pi * frequency_per_ns[:, None] * means[None, :])
    gaussian = np.exp(
        -2.0 * np.pi**2 * frequency_per_ns[:, None]**2 * sigmas[None, :]**2)
    transfer = (phase * gaussian) @ weights
    return float(np.trapezoid(np.abs(transfer)**2, frequency_hz) / bandwidth_hz)


def residual_delay_response(mixture, delay_ns, bandwidth_hz: float,
                            samples: int = 4097) -> np.ndarray:
    """计算未校正双镜时延下的相对零延迟相关响应。"""
    delays = np.asarray(delay_ns, dtype=float)
    frequency_hz = np.linspace(0.0, bandwidth_hz, samples)
    frequency_per_ns = frequency_hz * 1.0e-9
    weights = np.asarray(mixture["weights"], dtype=float)
    means = np.asarray(mixture["mean_delay_ns"], dtype=float)
    sigmas = np.asarray(mixture["std_delay_ns"], dtype=float)
    transfer = (
        np.exp(-2j * np.pi * frequency_per_ns[:, None] * means[None, :])
        * np.exp(-2.0 * np.pi**2 * frequency_per_ns[:, None]**2
                 * sigmas[None, :]**2)
    ) @ weights
    power = np.abs(transfer)**2
    normalization = np.trapezoid(power, frequency_hz)
    response = np.array([
        np.trapezoid(power * np.cos(2*np.pi*frequency_hz*tau*1.0e-9),
                     frequency_hz) / normalization
        for tau in delays.ravel()
    ]).reshape(delays.shape)
    return response


def convolve_pe_times(times_ns, amplitudes, sample_times_ns,
                      template_time_ns, template_amplitude_mv) -> np.ndarray:
    """在 ADC 时间网格上叠加每个平移后的标定 SPE 模板。"""
    output = np.zeros_like(np.asarray(sample_times_ns, dtype=float))
    for time_ns, amplitude in zip(times_ns, amplitudes):
        output += float(amplitude) * np.interp(
            output * 0.0 + sample_times_ns - time_ns,
            template_time_ns, template_amplitude_mv, left=0.0, right=0.0)
    return output


def digitize_adc(waveform_mv, bits: int, full_scale_mv: float) -> np.ndarray:
    """按对称满量程量化；``bits=full_scale=0`` 表示关闭 ADC 量化。"""
    if bits == 0 and full_scale_mv == 0.0:
        return np.asarray(waveform_mv, float).copy()
    if bits < 2 or full_scale_mv <= 0:
        raise ValueError("invalid ADC specification")
    half = full_scale_mv / 2.0
    levels = 2**bits
    step = full_scale_mv / levels
    clipped = np.clip(np.asarray(waveform_mv, float), -half, half - step)
    codes = np.rint((clipped + half) / step)
    return codes * step - half


def simulate_short_pair_waveforms(
        duration_ns: float,
        star_rate_hz: float,
        nsb_rate_hz: float,
        visibility2: float,
        instrument: Instrument,
        template_time_ns,
        template_amplitude_mv,
        delay_ns: float = 0.0,
        recovery_time_ns: float | None = None,
        charge_factors=None,
        optical_timing_mixture=None,
        padding_ns: float | None = None,
        seed: int = 20260824) -> dict[str, np.ndarray | float]:
    """模拟一段代表性的双望远镜 ADC 记录。

    恒星相关计数使用共享 Poisson 分量，使每个 ADC bin 的协方差等于
    ``r_star^2 * coherence_area * dt * |V|^2``；NSB 和剩余星光相互独立。
    这样得到正确的二阶矩，但不假装在数值上解析飞秒光场。
    """
    if duration_ns <= 0.0:
        raise ValueError("duration_ns must be positive")
    rng = np.random.default_rng(seed)
    dt_ns = instrument.sample_width_ns
    sample_edges = np.arange(0.0, duration_ns + dt_ns, dt_ns)
    centers = 0.5 * (sample_edges[:-1] + sample_edges[1:])
    # 模板以 p.e. 到达时刻为零点。必须在有效窗口外继续生成光子，否则
    # 靠近边界的 170 ns SPE 尾部会被人为截断。
    template_padding = max(abs(float(np.min(template_time_ns))),
                           abs(float(np.max(template_time_ns))))
    timing_padding = (0.0 if optical_timing_mixture is None else
                      5.0*float(optical_timing_mixture["rms_spread_ns"]))
    padding = (template_padding + abs(delay_ns) + timing_padding
               if padding_ns is None else float(padding_ns))
    if padding < template_padding:
        raise ValueError("padding_ns is shorter than the SPE template support")
    event_edges = np.arange(-padding, duration_ns + padding + dt_ns, dt_ns)
    event_bins = len(event_edges) - 1
    dt_s = dt_ns * 1.0e-9
    shared_mean = (star_rate_hz**2 * instrument.coherence_area_s * dt_s
                   * max(0.0, float(visibility2)))
    star_mean = star_rate_hz * dt_s
    nsb_mean = nsb_rate_hz * dt_s
    dark_mean = instrument.dark_count_rate_hz * dt_s
    if shared_mean > star_mean:
        raise ValueError("shared thermal component exceeds the star count")
    common = rng.poisson(shared_mean, event_bins)
    counts_a = common + rng.poisson(star_mean - shared_mean, event_bins)
    counts_b = common + rng.poisson(star_mean - shared_mean, event_bins)
    counts_a += rng.poisson(nsb_mean, event_bins)
    counts_b += rng.poisson(nsb_mean, event_bins)
    if dark_mean > 0.0:
        counts_a += rng.poisson(dark_mean, event_bins)
        counts_b += rng.poisson(dark_mean, event_bins)

    empirical_charge = (None if charge_factors is None
                        else np.asarray(charge_factors, dtype=float))
    if empirical_charge is not None and (
            empirical_charge.size == 0
            or np.any(~np.isfinite(empirical_charge))
            or np.any(empirical_charge <= 0.0)):
        raise ValueError("charge_factors must be finite and positive")

    def expand(counts, time_shift_ns):
        indices = np.repeat(np.arange(event_bins), counts)
        times = (event_edges[indices] + rng.random(len(indices))*dt_ns
                 + time_shift_ns)
        if optical_timing_mixture is not None:
            times += sample_optical_delays_ns(
                rng, len(times), optical_timing_mixture)
        cells = rng.integers(0, instrument.microcells_per_pixel, len(times))
        tau = (instrument.microcell_recovery_time_ns if recovery_time_ns is None
               else recovery_time_ns)
        recovery = apply_exponential_microcell_recovery(times, cells, tau)
        if empirical_charge is None:
            charge = np.ones(len(times), dtype=float)
        else:
            charge = rng.choice(empirical_charge, size=len(times), replace=True)
        return times, recovery, charge, recovery * charge

    times_a, recovery_a, charge_a, amplitudes_a = expand(counts_a, 0.0)
    times_b, recovery_b, charge_b, amplitudes_b = expand(counts_b, delay_ns)
    waveform_a = convolve_pe_times(
        times_a, amplitudes_a, centers, template_time_ns,
        template_amplitude_mv)
    waveform_b = convolve_pe_times(
        times_b, amplitudes_b, centers, template_time_ns,
        template_amplitude_mv)
    waveform_a += rng.normal(
        0.0, instrument.electronic_noise_rms_mv, len(centers))
    waveform_b += rng.normal(
        0.0, instrument.electronic_noise_rms_mv, len(centers))
    adc_a = digitize_adc(waveform_a, instrument.adc_bits,
                         instrument.adc_full_scale_mv)
    adc_b = digitize_adc(waveform_b, instrument.adc_bits,
                         instrument.adc_full_scale_mv)
    return {
        "sample_time_ns": centers,
        "adc_a_mv": adc_a,
        "adc_b_mv": adc_b,
        "analog_a_mv": waveform_a,
        "analog_b_mv": waveform_b,
        "pe_times_a_ns": times_a,
        "pe_times_b_ns": times_b,
        "recovery_a": recovery_a,
        "recovery_b": recovery_b,
        "charge_a": charge_a,
        "charge_b": charge_b,
        "shared_count_mean_per_bin": shared_mean,
        "analysis_duration_ns": float(duration_ns),
        "simulated_padding_each_side_ns": float(padding),
        "optical_timing_rms_ns": (
            0.0 if optical_timing_mixture is None
            else float(optical_timing_mixture["rms_spread_ns"])),
        "expected_sample_contrast": (
            instrument.coherence_area_s / dt_s * visibility2),
    }


def normalized_cross_correlation(left, right) -> tuple[np.ndarray, np.ndarray]:
    """计算两条去均值波形的归一化互相关及其整数采样滞后。"""
    left = np.asarray(left, float) - np.mean(left)
    right = np.asarray(right, float) - np.mean(right)
    scale = np.std(left) * np.std(right) * len(left)
    correlation = np.correlate(left, right, mode="full") / scale
    lags = np.arange(-len(left) + 1, len(left))
    return lags, correlation


def hbt_correlated_pair_rate_hz(star_rate_a_hz: float,
                                star_rate_b_hz: float,
                                coherence_area_s: float,
                                visibility2: float) -> float:
    """热光 HBT 超额光子对率：``R_pair=R1 R2 tau_eff |V|^2``。"""
    values = (star_rate_a_hz, star_rate_b_hz, coherence_area_s, visibility2)
    if any(not np.isfinite(value) or value < 0.0 for value in values):
        raise ValueError("HBT rates, coherence area, and visibility must be finite and >= 0")
    if visibility2 > 1.0:
        raise ValueError("squared visibility cannot exceed one")
    return star_rate_a_hz * star_rate_b_hz * coherence_area_s * visibility2


def _sample_active_sensor_positions(rng, count: int, device_path):
    """按 main 的 S17351 通道几何均匀抽取有效微单元位置。"""
    if count == 0:
        return np.empty(0), np.empty(0)
    if not device_path:
        return rng.uniform(-0.006, 0.006, count), rng.uniform(-0.006, 0.006, count)
    from config_io import read_key_value_config
    device = read_key_value_config(device_path)
    columns, rows = int(device["channel_columns"]), int(device["channel_rows"])
    width, height = float(device["channel_size_x_m"]), float(device["channel_size_y_m"])
    gap_x, gap_y = float(device["channel_gap_x_m"]), float(device["channel_gap_y_m"])
    total_x = columns*width + (columns-1)*gap_x
    total_y = rows*height + (rows-1)*gap_y
    column = rng.integers(0, columns, count)
    row = rng.integers(0, rows, count)
    x = -total_x/2 + column*(width+gap_x) + rng.random(count)*width
    y = -total_y/2 + row*(height+gap_y) + rng.random(count)*height
    return x, y


def simulate_hbt_primary_pe(
        rng: np.random.Generator,
        duration_ns: float,
        star_rate_hz: float | tuple[float, float],
        nsb_rate_hz: float | tuple[float, float],
        visibility2: float,
        instrument: Instrument = Instrument(),
        optical_timing_mixture=None,
        geometric_delay_ns: float = 0.0,
        padding_ns: float = 200.0,
        event_id: int = 1,
        pixel_id: int = 0,
        hbt_pair_rate_scale: float = 1.0) -> tuple[pd.DataFrame, dict]:
    """生成可直接交给 main ``run_camera_electronics`` 的双镜 p.e. 流。

    在光学相干时间远短于电子学分辨率时，热光的二阶超额相关可等价为
    稀疏相关光子对，其率为 ``R1*R2*tau_eff*|V|^2``。其余恒星光子和
    NSB 分别是独立 Poisson 流；DC 时间核独立作用于每一个恒星光子。
    """
    if duration_ns <= 0.0 or padding_ns < 0.0:
        raise ValueError("duration_ns must be > 0 and padding_ns must be >= 0")

    def pair(value):
        if np.isscalar(value):
            return float(value), float(value)
        if len(value) != 2:
            raise ValueError("rate must be a scalar or a two-element pair")
        return float(value[0]), float(value[1])

    star = pair(star_rate_hz)
    nsb = pair(nsb_rate_hz)
    if min(*star, *nsb) < 0.0:
        raise ValueError("photon rates must be >= 0")
    if not np.isfinite(hbt_pair_rate_scale) or hbt_pair_rate_scale <= 0.0:
        raise ValueError("hbt_pair_rate_scale must be finite and positive")
    physical_pair_rate = hbt_correlated_pair_rate_hz(
        star[0], star[1], instrument.coherence_area_s, visibility2)
    # 标定时可注入已知倍率的相关光子对，再把响应除回该倍率。总星光率
    # 保持不变，因此不会人为降低散粒噪声或提高最终物理曝光的 SNR。
    pair_rate = physical_pair_rate * hbt_pair_rate_scale
    if pair_rate > min(star):
        raise ValueError("correlated-pair approximation exceeds a marginal star rate")

    start_ns, end_ns = -padding_ns, duration_ns + padding_ns
    simulated_s = (end_ns-start_ns)*1.0e-9
    pair_count = int(rng.poisson(pair_rate*simulated_s))
    pair_centers = rng.uniform(start_ns, end_ns, pair_count)
    next_pair_id = np.arange(pair_count, dtype=np.int64)
    frames = []
    for telescope_id in (0, 1):
        star_single_count = int(rng.poisson((star[telescope_id]-pair_rate)*simulated_s))
        nsb_count = int(rng.poisson(nsb[telescope_id]*simulated_s))
        star_times = rng.uniform(start_ns, end_ns, star_single_count)
        pair_times = pair_centers.copy()
        if optical_timing_mixture is not None:
            star_times += sample_optical_delays_ns(
                rng, star_single_count, optical_timing_mixture)
            pair_times += sample_optical_delays_ns(
                rng, pair_count, optical_timing_mixture)
        if telescope_id == 1:
            star_times += geometric_delay_ns
            pair_times += geometric_delay_ns
        nsb_times = rng.uniform(start_ns, end_ns, nsb_count)
        times = np.concatenate((star_times, pair_times, nsb_times))
        origins = np.concatenate((
            np.full(star_single_count+pair_count, "cherenkov", object),
            np.full(nsb_count, "nsb", object)))
        pair_ids = np.concatenate((
            np.full(star_single_count, -1, np.int64), next_pair_id,
            np.full(nsb_count, -1, np.int64)))
        sensor_x, sensor_y = _sample_active_sensor_positions(
            rng, len(times), instrument.microcell_device_path)
        frames.append(pd.DataFrame({
            "event_id": event_id, "telescope_id": telescope_id,
            "pixel_id": pixel_id, "time_ns": times,
            "sensor_x_m": sensor_x, "sensor_y_m": sensor_y,
            "primary_pe": 1, "wavelength_nm": instrument.wavelength_nm,
            "origin": origins, "hbt_pair_id": pair_ids,
        }))
    hits = pd.concat(frames, ignore_index=True).sort_values(
        ["telescope_id", "time_ns"], ignore_index=True)
    metadata = {
        "duration_ns": float(duration_ns),
        "padding_ns": float(padding_ns),
        "star_rate_a_hz": star[0], "star_rate_b_hz": star[1],
        "nsb_rate_a_hz": nsb[0], "nsb_rate_b_hz": nsb[1],
        "visibility2": float(visibility2),
        "coherence_area_ps": instrument.coherence_area_s*1.0e12,
        "hbt_pair_rate_hz": physical_pair_rate,
        "injected_hbt_pair_rate_hz": pair_rate,
        "hbt_pair_rate_scale": float(hbt_pair_rate_scale),
        "expected_hbt_pairs_in_analysis_window": (
            physical_pair_rate*duration_ns*1.0e-9),
        "generated_hbt_pairs_with_padding": pair_count,
        "optical_timing_rms_ns": (0.0 if optical_timing_mixture is None else
                                  float(optical_timing_mixture["rms_spread_ns"])),
    }
    hits.attrs.update(metadata)
    return hits, metadata


def write_main_primary_pe_csv(hits: pd.DataFrame, path) -> Path:
    """写出 main ``run_camera_electronics`` 原生逐 p.e. CSV。"""
    columns = ["event_id", "telescope_id", "pixel_id", "time_ns",
               "sensor_x_m", "sensor_y_m", "primary_pe",
               "wavelength_nm", "origin"]
    missing = set(columns)-set(hits.columns)
    if missing:
        raise ValueError(f"primary p.e. table is missing {sorted(missing)}")
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    hits.to_csv(output, columns=columns, index=False)
    return output


def make_fast_spe_template(rise_ns: float = 0.15,
                           fall_ns: float = 0.80,
                           support_ns: float = 8.0,
                           step_ns: float = 0.01) -> tuple[np.ndarray, np.ndarray]:
    """生成峰值归一的快速差指数 SPE；参数可替换为实测前端响应。"""
    if not (0.0 < rise_ns < fall_ns and support_ns > 0.0 and step_ns > 0.0):
        raise ValueError("fast SPE needs 0 < rise < fall and positive support/step")
    time_ns = np.arange(0.0, support_ns+step_ns/2, step_ns)
    amplitude = np.exp(-time_ns/fall_ns)-np.exp(-time_ns/rise_ns)
    amplitude /= amplitude.max()
    return time_ns, amplitude


def _validate_waveform_instrument(instrument):
    for name in ("sipm_crosstalk_probability", "sipm_afterpulse_probability"):
        if getattr(instrument, name) != 0.0:
            raise NotImplementedError(f"waveform renderer does not simulate {name}")
    if not np.isclose(instrument.electronics_bandwidth_hz,
                      instrument.adc_sample_rate_hz/2, rtol=1e-12, atol=0):
        raise ValueError("waveform mode requires Nyquist electronics_bandwidth_hz; "
                         "use an SPE template with the desired hardware response")
    if instrument.optical_timing_in_effective_bandwidth:
        raise ValueError("an analytic effective bandwidth is not a waveform instrument")


def waveform_instrument_signature(instrument, template=None):
    """用数值参数和响应内容绑定标定；移动文件不会失效，内容改变必须重标定。"""
    values = asdict(instrument)
    values.pop("parameter_source")
    values.pop("detected_nsb_rate_hz")  # 实际光子率另外核对。
    for name in list(values):
        if name.endswith("_path"):
            path = values[name]
            values[name] = (None if path is None else
                            hashlib.sha256(Path(path).read_bytes()).hexdigest())
    if template is None and instrument.spe_template_path:
        template = load_measured_spe_template(instrument.spe_template_path)
    values.pop("spe_template_path")
    values["spe_template"] = (None if template is None else [
        hashlib.sha256(np.asarray(a, dtype="<f8").tobytes()).hexdigest()
        for a in template])
    values["renderer_version"] = "continuous-spe-v1"
    return hashlib.sha256(json.dumps(values, sort_keys=True, allow_nan=False).encode()).hexdigest()


def render_pe_waveform(
        rng: np.random.Generator, pe_times_ns, duration_ns: float,
        instrument: Instrument = Instrument(), sample_width_ns: float | None = None,
        template: tuple[np.ndarray, np.ndarray] | None = None,
        electronic_noise_rms_mv: float | None = None) -> dict:
    """直接采样平移后的连续SPE叠加，保留光电子到达时刻的亚采样间隔信息。"""
    _validate_waveform_instrument(instrument)
    times = np.asarray(pe_times_ns, dtype=float)
    if instrument.intrinsic_time_jitter_ns < 0:
        raise ValueError("intrinsic timing jitter must be nonnegative")
    if instrument.intrinsic_time_jitter_ns > 0:
        times = times+rng.normal(0., instrument.intrinsic_time_jitter_ns, len(times))
    dt = instrument.sample_width_ns if sample_width_ns is None else float(sample_width_ns)
    if (not np.isfinite(duration_ns) or not np.isfinite(dt) or duration_ns <= 0 or dt <= 0
            or times.ndim != 1 or not np.all(np.isfinite(times))):
        raise ValueError("duration, sample width and photon times must be finite and valid")
    if template is None:
        if not instrument.spe_template_path:
            raise ValueError("instrument has no measured SPE template")
        template = load_measured_spe_template(instrument.spe_template_path)
    template_time, template_amplitude = [np.asarray(value, float) for value in template]
    if (template_time.ndim != 1 or len(template_time) < 2
            or template_time.shape != template_amplitude.shape
            or not np.all(np.isfinite(template_time+template_amplitude))
            or not np.all(np.diff(template_time) > 0)):
        raise ValueError("SPE template needs finite amplitudes and strictly increasing times")
    centers = np.arange(dt/2, duration_ns, dt)
    samples = len(centers)
    if samples == 0:
        raise ValueError("duration contains no ADC sample center")
    cells = (rng.integers(0, instrument.microcells_per_pixel, len(times))
             if instrument.microcell_recovery_time_ns > 0 else np.zeros(len(times), dtype=int))
    recovery = apply_exponential_microcell_recovery(
        times, cells, instrument.microcell_recovery_time_ns)
    if instrument.charge_samples_path:
        empirical = load_empirical_charge_factors(instrument.charge_samples_path)
        charge = rng.choice(empirical, len(times))
    else:
        charge = np.ones(len(times))
    weights = recovery*charge
    analog = np.zeros(samples)
    offsets = np.arange(int(math.ceil(np.ptp(template_time)/dt))+1)
    # 只计算各个SPE支持区间内的样本，并限制临时数组的内存占用。
    chunk = max(1, 1_000_000//len(offsets))
    for begin in range(0, len(times), chunk):
        arrival = times[begin:begin+chunk]
        first = np.ceil((arrival+template_time[0])/dt-0.5).astype(np.int64)
        indices = first[:, None]+offsets
        relative = (indices+0.5)*dt-arrival[:, None]
        amplitudes = np.interp(relative, template_time, template_amplitude,
                               left=0., right=0.)*weights[begin:begin+chunk, None]
        valid = (indices >= 0) & (indices < samples)
        analog += np.bincount(indices[valid], weights=amplitudes[valid], minlength=samples)
    noise_rms = (instrument.electronic_noise_rms_mv
                 if electronic_noise_rms_mv is None else electronic_noise_rms_mv)
    if not np.isfinite(noise_rms) or noise_rms < 0:
        raise ValueError("electronic noise RMS must be finite and nonnegative")
    noisy = analog + rng.normal(0.0, noise_rms, samples) if noise_rms > 0 else analog.copy()
    return {
        "sample_time_ns": centers, "analog_mv": analog, "noisy_mv": noisy,
        "adc_mv": digitize_adc(noisy, instrument.adc_bits, instrument.adc_full_scale_mv),
        "recovery_fraction": recovery, "charge_factor": charge,
        "sample_width_ns": dt, "pe_count": len(times),
    }



def waveform_cross_correlation(left, right, sample_width_ns: float,
                               max_lag_ns: float) -> tuple[np.ndarray, np.ndarray]:
    """FFT 计算有限波形的无偏归一互相关。"""
    left = np.asarray(left, float)
    right = np.asarray(right, float)
    if left.shape != right.shape or left.ndim != 1 or left.size < 2:
        raise ValueError("waveforms must be equal non-trivial one-dimensional arrays")
    if sample_width_ns <= 0.0 or max_lag_ns < 0.0:
        raise ValueError("sample width must be positive and max lag non-negative")
    x, y = left-left.mean(), right-right.mean()
    try:
        from scipy.signal import correlate
        raw = correlate(x, y, mode="full", method="fft")
    except ImportError:  # pragma: no cover
        raw = np.correlate(x, y, mode="full")
    lags = np.arange(-len(x)+1, len(x))
    overlap = len(x)-np.abs(lags)
    scale = np.std(x)*np.std(y)
    correlation = raw/overlap/scale if scale > 0.0 else np.zeros_like(raw)
    keep = np.abs(lags*sample_width_ns) <= max_lag_ns
    return lags[keep]*sample_width_ns, correlation[keep]


def calibrate_waveform_gls(
        lags_ns, null_correlations, signal_correlations,
        block_duration_s: float, star_rate_hz: float,
        background_rate_hz: float, calibration_visibility2: float = 1.0,
        hbt_pair_rate_scale: float = 1.0,
        covariance_shrinkage: float = 0.05) -> WaveformGLSCalibration:
    """由无信号和已知注入记录标定完整峰 GLS 估计器。

    ``signal_correlations`` 的注入强度为
    ``calibration_visibility2*hbt_pair_rate_scale``。增强只用于更快得到
    平均响应模板，函数会除回该倍率；噪声协方差始终来自物理的零信号记录。
    """
    lags = np.asarray(lags_ns, float)
    null = np.asarray(null_correlations, float)
    signal = np.asarray(signal_correlations, float)
    if (null.ndim != 2 or signal.ndim != 2 or null.shape[1:] != lags.shape
            or signal.shape[1:] != lags.shape):
        raise ValueError("correlation records must have shape (records, lags)")
    if null.shape[0] < 3 or signal.shape[0] < 1:
        raise ValueError("calibration needs at least 3 null and 1 signal records")
    scale = calibration_visibility2*hbt_pair_rate_scale
    if (block_duration_s <= 0.0 or star_rate_hz < 0.0
            or background_rate_hz < 0.0 or scale <= 0.0
            or not 0.0 <= covariance_shrinkage <= 1.0):
        raise ValueError("invalid waveform GLS calibration parameters")

    null_mean = null.mean(axis=0)
    peak = (signal.mean(axis=0)-null_mean)/scale
    covariance = np.atleast_2d(np.cov(null, rowvar=False, ddof=1))
    diagonal = np.diag(np.diag(covariance))
    covariance = ((1.0-covariance_shrinkage)*covariance
                  + covariance_shrinkage*diagonal)
    if np.any(np.diag(covariance) <= 0.0) or not np.any(np.abs(peak) > 0.0):
        raise ValueError("degenerate waveform GLS calibration")
    return WaveformGLSCalibration(
        lags_ns=lags, null_mean=null_mean,
        peak_per_visibility2=peak, covariance_per_block=covariance,
        block_duration_s=float(block_duration_s),
        star_rate_hz=float(star_rate_hz),
        background_rate_hz=float(background_rate_hz),
        calibration_visibility2=float(calibration_visibility2),
        hbt_pair_rate_scale=float(hbt_pair_rate_scale),
        null_records=int(null.shape[0]), signal_records=int(signal.shape[0]),
        covariance_shrinkage=float(covariance_shrinkage))


def waveform_gls_weights(
        calibration: WaveformGLSCalibration,
        exposure_s: float) -> tuple[np.ndarray, float]:
    """返回GLS权重及条件sigma_P；有限标定精度另见标定diagnostics。"""
    covariance = calibration.covariance(exposure_s)
    peak = calibration.peak_per_visibility2
    solved = np.linalg.solve(covariance, peak)
    information = float(peak @ solved)
    if not np.isfinite(information) or information <= 0.0:
        raise ValueError("waveform GLS covariance has no positive information")
    return solved/information, 1.0/math.sqrt(information)


def estimate_visibility2_gls(
        correlations, calibration: WaveformGLSCalibration,
        exposure_s: float | None = None) -> tuple[np.ndarray | float, float]:
    """从一条或多条完整互相关峰估计 ``P=|V|^2``。"""
    exposure = (calibration.block_duration_s if exposure_s is None
                else float(exposure_s))
    weight, sigma = waveform_gls_weights(calibration, exposure)
    values = np.asarray(correlations, float)
    if values.shape[-1] != len(calibration.lags_ns):
        raise ValueError("correlation lag grid does not match calibration")
    estimate = (values-calibration.null_mean) @ weight
    return (float(estimate) if np.ndim(estimate) == 0 else estimate), sigma


def simulate_waveform_gls_calibration(
        instrument: Instrument, source_ab_magnitude: float = 2.0,
        block_duration_ns: float = 20_000.0, max_lag_ns: float = 200.0,
        null_records: int = 1024, signal_records: int = 256,
        calibration_visibility2: float = 1.0,
        hbt_pair_rate_scale: float = 10_000.0,
        covariance_shrinkage: float | str = "auto",
        nsb_rate_hz: float | None = None,
        spe_template=None,
        seed: int = 20260824) -> tuple[WaveformGLSCalibration, dict]:
    """用 main 同参数的单像素短波形生成完整峰 GLS 标定。

    输出峰模板和零信号协方差；后续长曝光只需按 ``1/T`` 缩放协方差，
    不需要在内存中保存整段20分钟波形。
    """
    if block_duration_ns <= 0.0 or max_lag_ns < 0.0:
        raise ValueError("block duration and maximum lag are invalid")
    _validate_waveform_instrument(instrument)
    if not np.isclose(block_duration_ns/instrument.sample_width_ns,
                      round(block_duration_ns/instrument.sample_width_ns), atol=1e-10, rtol=0):
        raise ValueError("calibration block duration must contain an integer number of samples")
    if null_records < 16 or signal_records < 16:
        raise ValueError("calibration needs at least 16 null and 16 signal records")
    if spe_template is None and not instrument.spe_template_path:
        raise ValueError("instrument has no measured SPE template")

    rng = np.random.default_rng(seed)
    star_rate = detected_star_rate_hz(source_ab_magnitude, instrument)
    nsb_rate = (instrument.detected_nsb_rate_hz if nsb_rate_hz is None
                else float(nsb_rate_hz))
    background_rate = nsb_rate+instrument.dark_count_rate_hz
    template = (load_measured_spe_template(instrument.spe_template_path)
                if spe_template is None else spe_template)
    timing = (None if not instrument.optical_timing_kernel_path else
              load_optical_timing_mixture(instrument.optical_timing_kernel_path))
    timing_padding = 0.0 if timing is None else 5.0*timing["rms_spread_ns"]
    padding_ns = (max(abs(template[0].min()), abs(template[0].max()))+timing_padding
                  +5*instrument.intrinsic_time_jitter_ns)

    def one_record(visibility2, pair_scale):
        hits, _ = simulate_hbt_primary_pe(
            rng, block_duration_ns, star_rate, background_rate, visibility2,
            instrument, optical_timing_mixture=timing, padding_ns=padding_ns,
            hbt_pair_rate_scale=pair_scale)
        waveforms = []
        for telescope_id in (0, 1):
            times = hits.loc[hits.telescope_id == telescope_id, "time_ns"]
            rendered = render_pe_waveform(
                rng, times.to_numpy(float), block_duration_ns, instrument,
                template=template)
            waveforms.append(rendered["adc_mv"])
        return waveform_cross_correlation(
            waveforms[0], waveforms[1], instrument.sample_width_ns,
            max_lag_ns)

    null = []
    for _ in range(null_records):
        lags, correlation = one_record(0.0, 1.0)
        null.append(correlation)
    signal = []
    for _ in range(signal_records):
        signal_lags, correlation = one_record(
            calibration_visibility2, hbt_pair_rate_scale)
        if not np.array_equal(signal_lags, lags):
            raise RuntimeError("waveform correlation lag grids changed")
        signal.append(correlation)

    null = np.asarray(null)
    signal = np.asarray(signal)
    null_split = null_records//2
    signal_split = signal_records//2
    null_scale_start = 3*null_records//4
    signal_test_start = 3*signal_records//4
    candidates = (1e-4, 1e-3, 1e-2, 5e-2) if covariance_shrinkage == "auto" else (float(covariance_shrinkage),)
    if any(not np.isfinite(x) or not 0 <= x <= 1 for x in candidates):
        raise ValueError("covariance_shrinkage must be auto or a number in [0, 1]")
    calibration = calibrate_waveform_gls(
        lags, null[:null_split], signal[:signal_split], block_duration_ns*1.0e-9,
        star_rate, background_rate, calibration_visibility2,
        hbt_pair_rate_scale, candidates[0])

    # 对平稳波形，两个相关滞后点的协方差只依赖滞后之差。沿每条对角线
    # 平均可大幅减少有限标定样本造成的协方差噪声，再轻微向对角阵收缩。
    raw_covariance = np.cov(null[:null_split], rowvar=False, ddof=1)
    lag_count = len(lags)
    toeplitz_values = np.array([
        np.mean(np.diag(raw_covariance, offset))
        for offset in range(lag_count)])
    stationary_covariance = np.fromfunction(
        lambda i, j: toeplitz_values[np.abs(i-j).astype(int)],
        (lag_count, lag_count), dtype=int)
    calibration.null_mean = np.zeros(lag_count)
    calibration.peak_per_visibility2 = 0.5*(
        calibration.peak_per_visibility2+calibration.peak_per_visibility2[::-1])
    # 线性响应的连续模板保留亚采样形状；独立波形仍决定协方差、增益和误差。
    if (spe_template is None and nsb_rate == instrument.detected_nsb_rate_hz
            and not any((instrument.adc_bits, instrument.microcell_recovery_time_ns,
                         instrument.sipm_crosstalk_probability, instrument.sipm_afterpulse_probability))):
        from sii_validation import analytic_waveform_calibration
        reference, _ = analytic_waveform_calibration(instrument, source_ab_magnitude,
            block_duration_ns=block_duration_ns, max_lag_ns=max_lag_ns)
        calibration.peak_per_visibility2 = reference.peak_per_visibility2.copy()
        calibration.phase_model = "linear_spe"
    diagonal = np.diag(np.diag(stationary_covariance))
    injected_visibility2 = calibration_visibility2*hbt_pair_rate_scale
    selection = []
    for shrinkage in candidates:
        covariance = (1-shrinkage)*stationary_covariance+shrinkage*diagonal
        eigenvalue, eigenvector = np.linalg.eigh(covariance)
        eigenvalue = np.maximum(eigenvalue, eigenvalue.max()*1e-10)
        calibration.covariance_per_block = (eigenvector*eigenvalue) @ eigenvector.T
        tune_null, _ = estimate_visibility2_gls(null[null_split:null_scale_start], calibration)
        gain_values, _ = estimate_visibility2_gls(signal[signal_split:signal_test_start], calibration)
        gain = float(gain_values.mean()/injected_visibility2)
        score = float(tune_null.std(ddof=1)/gain) if gain > 0 else float("inf")
        selection.append((score, shrinkage, calibration.covariance_per_block.copy(), gain))
    score, shrinkage, covariance, gain = min(selection, key=lambda row: row[0])
    if not np.isfinite(score):
        raise ValueError("calibration injection has no positive held-out response")
    calibration.covariance_per_block = covariance
    calibration.covariance_shrinkage = shrinkage
    gain_values, _ = estimate_visibility2_gls(signal[signal_split:signal_test_start], calibration)
    calibration.response_relative_uncertainty = float(
        gain_values.std(ddof=1)/np.sqrt(len(gain_values))/abs(gain_values.mean()))
    calibration.peak_per_visibility2 *= gain
    calibration.null_records = null_records
    calibration.signal_records = signal_records
    calibration.instrument_signature = waveform_instrument_signature(instrument, template)

    # 用未参加正则化选择的记录确定固定权重的噪声尺度。
    validation_null, predicted_sigma = estimate_visibility2_gls(null[null_scale_start:], calibration)
    empirical_sigma = float(np.std(validation_null, ddof=1))
    covariance_scale = (empirical_sigma/predicted_sigma)**2
    calibration.covariance_per_block *= covariance_scale
    calibration.sigma_relative_uncertainty = 1/math.sqrt(2*(len(validation_null)-1))
    validation_signal, _ = estimate_visibility2_gls(signal[signal_test_start:], calibration)
    response_closure = float(np.mean(validation_signal)/injected_visibility2)
    _, sigma_block = waveform_gls_weights(calibration, calibration.block_duration_s)
    diagnostics = {
        "null_correlations": null,
        "signal_correlations": signal,
        "mean_signal_correlation": np.mean(signal, axis=0),
        "sigma_visibility2_per_block": sigma_block,
        "validation_null_mean_visibility2": float(np.mean(validation_null)),
        "validation_null_std_visibility2": empirical_sigma,
        "validation_response_closure": response_closure,
        "covariance_validation_scale": float(covariance_scale),
        "null_training_records": null_split,
        "signal_training_records": signal_split,
        # 正态近似下，有限尺度标定本身的相对标准误差；不是额外的独立UV噪声。
        "null_scale_records": null_records-null_scale_start,
        "null_scale_start": null_scale_start,
        "signal_test_start": signal_test_start,
        "selected_covariance_shrinkage": shrinkage,
        "shrinkage_selection": [{"shrinkage": row[1], "tuning_sigma": row[0]} for row in selection],
        "response_gain_correction": gain,
        "response_relative_calibration_error": calibration.response_relative_uncertainty,
        "sigma_relative_calibration_error_approx": calibration.sigma_relative_uncertainty,
        "padding_ns": float(padding_ns),
        "source_ab_magnitude": float(source_ab_magnitude),
        "physical_hbt_pair_rate_hz": hbt_correlated_pair_rate_hz(
            star_rate, star_rate, instrument.coherence_area_s, 1.0),
    }
    return calibration, diagnostics


def combine_visibility_estimates(values, sigmas) -> tuple[float, float]:
    """用逆方差合并多个独立 ``P`` 估计。"""
    estimates = np.asarray(values, float)
    uncertainty = np.asarray(sigmas, float)
    if (estimates.shape != uncertainty.shape or estimates.size == 0
            or np.any(~np.isfinite(estimates))
            or np.any(~np.isfinite(uncertainty))
            or np.any(uncertainty <= 0.0)):
        raise ValueError("visibility estimates and sigmas are invalid")
    inverse_variance = 1.0/uncertainty**2
    return (float(np.sum(estimates*inverse_variance)/np.sum(inverse_variance)),
            float(1.0/math.sqrt(np.sum(inverse_variance))))


def sample_waveform_gls_visibility2(
        visibility2, calibration: WaveformGLSCalibration,
        exposure_s: float,
        rng: np.random.Generator | None = None) -> tuple[np.ndarray, float]:
    """抽样完整相关峰的线性GLS投影，得到长曝光 ``P``。

    这是短波形到长曝光的充分统计量：均值为 ``C0+P*k``，协方差由
    标定块按 ``T_block/T`` 缩放。没有生成无法存储的连续20分钟数组。
    """
    truth = np.asarray(visibility2, float)
    if np.any(~np.isfinite(truth)):
        raise ValueError("visibility2 must be finite")
    generator = np.random.default_rng() if rng is None else rng
    _, sigma = waveform_gls_weights(calibration, exposure_s)
    # 同一标定用于整批UV点：增益误差只抽一次，不能当作每条基线的独立噪声。
    gain_sigma = calibration.response_relative_uncertainty
    gain = 1.+generator.normal(0., gain_sigma) if gain_sigma > 0 else 1.
    # 多维高斯相关函数的线性GLS投影严格等价于这个一维高斯抽样。
    return gain*truth + generator.normal(0.0, sigma, size=truth.shape), sigma


def _load_layout(layout) -> pd.DataFrame:
    if isinstance(layout, (str, Path)) and Path(layout).suffix.lower() == '.input':
        from sii_layout import read_corsika_layout
        frame = read_corsika_layout(layout)
    else:
        frame = pd.read_csv(layout) if isinstance(layout, (str, Path)) else layout.copy()
    aliases = {"position_x_m": "east_m", "position_y_m": "north_m",
               "position_z_m": "up_m"}
    frame = frame.rename(columns=aliases)
    required = {"east_m", "north_m", "up_m"}
    if not required.issubset(frame):
        raise ValueError(f"layout is missing {sorted(required-set(frame.columns))}")
    if not np.all(np.isfinite(frame[['east_m', 'north_m', 'up_m']].to_numpy(float))):
        raise ValueError('望远镜ENU坐标必须为有限数值')
    if "telescope_id" not in frame:
        frame["telescope_id"] = np.arange(1, len(frame)+1)
    if "name" not in frame:
        frame["name"] = frame.telescope_id.map(lambda value: f"TEL.{value}")
    if len(frame) < 2 or frame.telescope_id.duplicated().any():
        raise ValueError("layout needs at least two unique telescopes")
    return frame.reset_index(drop=True)


def generate_uvw(layout, observation: Observation,
                 instrument: Instrument = Instrument()) -> pd.DataFrame:
    """生成一次完整观测的所有基线和时间段，而不是单个 UV 点。

    对每对望远镜 ``B_ij = r_j-r_i``，使用
    ``(u,v,w)=(B·u_hat, B·v_hat, B·s_hat)``，几何时延为 ``w/c``。
    """
    telescopes = _load_layout(layout)
    if observation.hours_per_night <= 0 or observation.segment_s <= 0:
        raise ValueError("observation duration and segment_s must be positive")
    if observation.nights <= 0:
        raise ValueError("nights must be positive")
    segment_count = int(round(
        observation.hours_per_night*3600.0/observation.segment_s))
    if not np.isclose(segment_count*observation.segment_s,
                      observation.hours_per_night*3600., rtol=1e-12, atol=1e-9):
        raise ValueError("observation duration must be an integer number of segments")
    if segment_count < 1:
        raise ValueError("observation has no integration segment")
    half_segment_h = observation.segment_s/7200.0
    hour_angles_h = np.linspace(
        -observation.hours_per_night/2 + half_segment_h,
        observation.hours_per_night/2 - half_segment_h,
        segment_count)*86400.0/SIDEREAL_DAY_S
    lat = math.radians(observation.site_lat_deg)
    dec = math.radians(observation.source_dec_deg)
    # 时角H=LST-RA，西向为正；地平线以下不能继续累计可见光有效曝光。
    for hour in (-observation.hours_per_night/2, 0., observation.hours_per_night/2):
        if source_direction_enu(hour*3600.*2*math.pi/SIDEREAL_DAY_S, dec, lat)[2] <= 0:
            raise ValueError("source is below the horizon during the observation")
    rows = []
    for i, j in combinations(range(len(telescopes)), 2):
        first, second = telescopes.iloc[i], telescopes.iloc[j]
        baseline = (
            second[["east_m", "north_m", "up_m"]].to_numpy(float)
            - first[["east_m", "north_m", "up_m"]].to_numpy(float))
        for segment, hour_angle_h in enumerate(hour_angles_h):
            hour_angle = hour_angle_h*math.pi/12.0
            u_m, v_m, w_m = uvw_from_enu(baseline, hour_angle, dec, lat)
            rows.append({
                "telescope_i": str(first["name"]),
                "telescope_j": str(second["name"]),
                "telescope_i_index": i, "telescope_j_index": j,
                "segment": segment, "hour_angle_h": hour_angle_h,
                "baseline_east_m": baseline[0],
                "baseline_north_m": baseline[1],
                "baseline_up_m": baseline[2],
                "baseline_m": float(np.linalg.norm(baseline)),
                "u_m": u_m, "v_m": v_m, "w_m": w_m,
                "u_lambda": u_m/instrument.wavelength_m,
                "v_lambda": v_m/instrument.wavelength_m,
                "projected_baseline_m": math.hypot(u_m, v_m),
                "geometric_delay_ns": w_m/C_M_S*1.0e9,
            })
    return pd.DataFrame(rows)


def segment_uv_samples(uvw, observation: Observation,
                       instrument: Instrument = Instrument()):
    """只由基线和时间生成积分子采样；模拟和重建共用，不接触天体真值。"""
    frame = uvw
    required = {"hour_angle_h", "baseline_east_m", "baseline_north_m",
                "baseline_up_m"}
    if not required.issubset(frame):
        missing = sorted(required-set(frame.columns))
        raise ValueError(
            f"uvw is missing segment-average columns {missing}")
    samples = observation.visibility_subsamples_per_segment
    if not isinstance(samples, int) or samples < 1:
        raise ValueError("visibility_subsamples_per_segment must be a positive integer")

    baseline = frame[["baseline_east_m", "baseline_north_m",
                      "baseline_up_m"]].to_numpy(float)
    center_hour = frame.hour_angle_h.to_numpy(float)
    offsets_h = ((np.arange(samples)+0.5)/samples-0.5
                 )*observation.segment_s/3600.0*86400.0/SIDEREAL_DAY_S
    lat = math.radians(observation.site_lat_deg)
    dec = math.radians(observation.source_dec_deg)
    sl, cl = math.sin(lat), math.cos(lat)
    sd, cd = math.sin(dec), math.cos(dec)
    us, vs = [], []
    for offset_h in offsets_h:
        hour = (center_hour+offset_h)*math.pi/12.0
        sh, ch = np.sin(hour), np.cos(hour)
        u_m = (baseline[:, 0]*ch - baseline[:, 1]*sl*sh
               + baseline[:, 2]*cl*sh)
        v_m = (baseline[:, 0]*sd*sh
               + baseline[:, 1]*(cd*cl+sd*ch*sl)
               + baseline[:, 2]*(cd*sl-sd*ch*cl))
        u, v = u_m/instrument.wavelength_m, v_m/instrument.wavelength_m
        us.append(u)
        vs.append(v)
    u, v = np.stack(us, axis=1), np.stack(vs, axis=1)
    wavelength = np.asarray(instrument.visibility_wavelength_nm or (instrument.wavelength_nm,))
    scale = instrument.wavelength_nm/wavelength
    return ((u[..., None]*scale).reshape(len(frame), -1),
            (v[..., None]*scale).reshape(len(frame), -1))


def segment_sampling_weights(observation, instrument):
    """时间中点等权，波长按HBT相关面积加权；生成和重建必须使用同一权重。"""
    spectral = np.asarray(instrument.visibility_spectral_weights or (1.,), float)
    if (spectral.size != len(instrument.visibility_wavelength_nm or (instrument.wavelength_nm,))
            or np.any(spectral < 0) or not np.isclose(spectral.sum(), 1.)):
        raise ValueError("invalid spectral quadrature weights")
    count = observation.visibility_subsamples_per_segment
    return np.tile(spectral/count, count)


def segment_averaged_visibility2(
        uvw, source, observation: Observation,
        instrument: Instrument = Instrument(), source_case: str = "binary"):
    """计算每个积分段沿地球自转轨迹平均后的 ``|V|²`` 期望值。"""
    u, v = segment_uv_samples(uvw, observation, instrument)
    return np.abs(source_visibility(u, v, source, source_case))**2 @ segment_sampling_weights(observation, instrument)


def _electronics_correlation_efficiency(instrument, total_rate_hz):
    """由主程序 SPE 模板估算加性电子噪声造成的相关效率。"""
    if not instrument.spe_template_path:
        return 1.0
    template_t, template_v = load_measured_spe_template(
        instrument.spe_template_path)
    charge_second_moment = instrument.excess_noise_factor**2
    shot_variance = (total_rate_hz/1.0e9 * charge_second_moment
                     * np.trapezoid(template_v**2, template_t))
    adc_step = (0.0 if instrument.adc_bits <= 0
                else instrument.adc_full_scale_mv/2**instrument.adc_bits)
    additive = instrument.electronic_noise_rms_mv**2 + adc_step**2/12.0
    return float(shot_variance/(shot_variance+additive))


def simulate_uv_observation(
        uvw, source, observation: Observation,
        instrument: Instrument = Instrument(), seed: int = 20260824,
        source_case: str = "binary", nsb_multiplier: float = 1.0,
        electronics_case: str = "reference", estimator: str = "analytic",
        waveform_calibration: WaveformGLSCalibration | None = None,
        ) -> tuple[pd.DataFrame, dict]:
    """把完整 UVW 表转换成带误差的长曝光 ``|V|²`` 测量。

    主报告显式选waveform_gls：使用匹配当前仪器和光子率的短波形标定，
    按积分时间缩放协方差并抽样P；不是重新产生数小时ADC波形。
    为兼容保留的analytic分支才另外抽样每镜每段的Poisson计数及共享状态，
    并使用解析SNR计算sigma。两条估计器不能混称；GLS未标定的系统误差会报错。
    """
    frame = uvw.copy().reset_index(drop=True)
    sample_u, sample_v = segment_uv_samples(frame, observation, instrument)
    frame["uv_samples_u"] = list(map(tuple, sample_u))
    frame["uv_samples_v"] = list(map(tuple, sample_v))
    frame["uv_samples_weight"] = [tuple(segment_sampling_weights(observation, instrument))]*len(frame)
    required = {"u_lambda", "v_lambda", "segment",
                "telescope_i_index", "telescope_j_index"}
    if not required.issubset(frame):
        raise ValueError(f"uvw is missing {sorted(required-set(frame.columns))}")
    if nsb_multiplier < 0:
        raise ValueError("nsb_multiplier must be non-negative")

    u = frame.u_lambda.to_numpy(float)
    v = frame.v_lambda.to_numpy(float)
    center_truth = np.abs(source_visibility(u, v, source, source_case))**2
    truth = segment_averaged_visibility2(
        frame, source, observation, instrument, source_case)
    frame["visibility2_center"] = center_truth
    frame["visibility2_true"] = truth
    frame["segment_time_smearing_delta"] = truth-center_truth
    frame["baseline_integration_s"] = observation.segment_s*observation.nights

    rng = np.random.default_rng(seed)
    row_i = frame.telescope_i_index.to_numpy(int)
    row_j = frame.telescope_j_index.to_numpy(int)
    row_segment = frame.segment.to_numpy(int)
    telescope_count = int(max(row_i.max(), row_j.max())+1)
    segment_count = int(row_segment.max()+1)
    pairs = list(zip(row_i, row_j))
    pair_zero = {pair: rng.normal(0.0, instrument.baseline_zero_point_rms)
                 for pair in set(pairs)}
    static_gain = rng.normal(
        0.0, instrument.telescope_gain_calibration_rms, telescope_count)

    def ar1_lognormal(rms, rho=0.85):
        innovation = rng.normal(size=(segment_count, telescope_count))
        state = np.zeros_like(innovation)
        state[0] = innovation[0]
        for index in range(1, segment_count):
            state[index] = (rho*state[index]
                            + math.sqrt(1-rho**2)*innovation[index])
        return np.exp(rms*state-0.5*rms**2)

    nsb_rate = nsb_multiplier*instrument.detected_nsb_rate_hz
    dark_rate = instrument.dark_count_rate_hz
    star_rate = detected_star_rate_hz(source.ab_magnitude, instrument)
    if estimator == "waveform_gls":
        if waveform_calibration is None:
            raise ValueError("waveform_gls estimator needs waveform_calibration")
        _validate_waveform_instrument(instrument)
        if waveform_calibration.instrument_signature != waveform_instrument_signature(instrument):
            raise ValueError("waveform calibration instrument/response signature mismatch; recalibrate")
        expected_background = nsb_rate+dark_rate
        if (not np.isclose(waveform_calibration.star_rate_hz, star_rate,
                           rtol=1e-9)
                or not np.isclose(waveform_calibration.background_rate_hz,
                                  expected_background, rtol=1e-9)):
            raise ValueError(
                "waveform calibration rates do not match this source/background")
        unsupported = {
            "telescope_gain_calibration_rms": instrument.telescope_gain_calibration_rms,
            "per_night_gain_rms": instrument.per_night_gain_rms,
            "baseline_zero_point_rms": instrument.baseline_zero_point_rms,
            "residual_timing_rms_ns": instrument.residual_timing_rms_ns,
            "transparency_fractional_rms": instrument.transparency_fractional_rms,
            "nsb_fractional_rms": instrument.nsb_fractional_rms,
        }
        enabled = [name for name, value in unsupported.items() if value != 0.0]
        if enabled:
            raise NotImplementedError(
                "waveform_gls systematics are not calibrated: " + ", ".join(enabled))
        exposure_s = observation.segment_s*observation.nights
        from sii_performance import phase_template_bank, tracked_segment_precision
        bank = phase_template_bank(instrument, source.ab_magnitude,
            block_duration_ns=waveform_calibration.block_duration_s*1e9,
            calibration=waveform_calibration)
        sigma_stat = tracked_segment_precision(frame, observation, instrument, bank)
        gain = 1+rng.normal(0., waveform_calibration.response_relative_uncertainty)
        measured = gain*truth+rng.normal(0., sigma_stat)
        sigma_value = float(np.median(sigma_stat))
        frame["visibility2_measured"] = measured
        frame["sigma_visibility2_stat"] = sigma_stat
        frame["sigma_visibility2"] = sigma_stat
        frame["sigma_visibility2_calibration"] = (
            np.abs(measured)*waveform_calibration.response_relative_uncertainty)
        frame["calibration_relative_sigma"] = waveform_calibration.response_relative_uncertainty
        frame["calibration_id"] = hashlib.sha256(
            waveform_calibration.peak_per_visibility2.tobytes()).hexdigest()
        endpoint = observation.hours_per_night*3600.*math.pi/SIDEREAL_DAY_S
        altitude = [math.asin(source_direction_enu(
            sign*endpoint, math.radians(observation.source_dec_deg),
            math.radians(observation.site_lat_deg))[2]) for sign in (-1, 1)]
        metadata = {
            "estimator": "waveform_gls",
            "delay_processing": "integer ADC shifts with continuous residual-phase templates; equal-time variance integration",
            "uncertainty_model": "independent statistical errors plus one shared Gaussian calibration gain, profiled in reconstruction",
            "calibration_response_relative_uncertainty": waveform_calibration.response_relative_uncertainty,
            "calibration_sigma_relative_uncertainty": waveform_calibration.sigma_relative_uncertainty,
            "hours_per_night": observation.hours_per_night,
            "nights": observation.nights,
            "total_integration_hours": observation.hours_per_night*observation.nights,
            "integration_s_per_uv_measurement": exposure_s,
            "visibility_subsamples_per_segment": (
                observation.visibility_subsamples_per_segment),
            "max_abs_segment_time_smearing": float(
                np.max(np.abs(truth-center_truth))),
            "source_ab_magnitude": source.ab_magnitude,
            "detected_star_rate_MHz": star_rate/1.0e6,
            "source_case": source_case,
            "nsb_multiplier": nsb_multiplier,
            "nsb_rate_MHz": nsb_rate/1.0e6,
            "dark_count_rate_MHz": dark_rate/1.0e6,
            "electronics_case": "waveform_gls_calibrated",
            "waveform_block_duration_us": (
                waveform_calibration.block_duration_s*1.0e6),
            "waveform_lag_samples": len(waveform_calibration.lags_ns),
            "waveform_null_records": waveform_calibration.null_records,
            "waveform_signal_records": waveform_calibration.signal_records,
            "waveform_hbt_pair_rate_scale": (
                waveform_calibration.hbt_pair_rate_scale),
            "minimum_altitude_deg": math.degrees(min(altitude)),
            "uv_measurements": len(frame),
            "sigma_visibility2_stat": float(sigma_value),
            "sigma_visibility2_total": float(sigma_value),
            "instrument_parameter_source": instrument.parameter_source,
        }
        frame.attrs.update(metadata)
        return frame, metadata
    if estimator != "analytic":
        raise ValueError("estimator must be analytic or waveform_gls")
    unsupported = ('telescope_gain_calibration_rms', 'per_night_gain_rms',
                   'residual_timing_rms_ns', 'transparency_fractional_rms', 'nsb_fractional_rms')
    enabled = [name for name in unsupported if getattr(instrument, name) != 0]
    if enabled:
        raise NotImplementedError('analytic drift likelihood is not calibrated: '+', '.join(enabled))

    if electronics_case == "ideal":
        charge_rms, electronics_efficiency = 1.0, 1.0
    elif electronics_case == "reference":
        charge_rms = instrument.excess_noise_factor
        electronics_efficiency = _electronics_correlation_efficiency(
            instrument, star_rate+nsb_rate+dark_rate)
        if instrument.optical_timing_in_effective_bandwidth:
            electronics_efficiency = 1.0  # 匹配积分已包含电子/量化噪声，不再扣一次。
    else:
        raise ValueError("electronics_case must be reference or ideal")

    if instrument.optical_timing_kernel_path:
        timing = load_optical_timing_mixture(
            instrument.optical_timing_kernel_path)
        optical_efficiency = (
            1.0 if instrument.optical_timing_in_effective_bandwidth
            else optical_timing_transfer_efficiency(
                timing, instrument.electronics_bandwidth_hz))
        timing_grid = np.linspace(-3.0, 3.0, 2401)
        timing_response = residual_delay_response(
            timing, timing_grid, instrument.electronics_bandwidth_hz)
    else:
        optical_efficiency = 1.0
        timing_grid, timing_response = np.array([-1.0, 1.0]), np.ones(2)

    nightly_values, nightly_variances = [], []
    photon_relative_rms, timing_factors, gain_factors = [], [], []
    for _night in range(observation.nights):
        transparency = ar1_lognormal(instrument.transparency_fractional_rms)
        nsb_factor = ar1_lognormal(instrument.nsb_fractional_rms)
        star_counts = rng.poisson(
            star_rate*transparency*observation.segment_s)
        nsb_counts = rng.poisson(
            nsb_rate*nsb_factor*observation.segment_s)
        # 暗计数与恒星/NSB一样按“时间段×望远镜”独立抽样；同一格计数
        # 再由该望远镜在该时间段参与的全部基线共享。
        dark_counts = rng.poisson(
            dark_rate*observation.segment_s,
            size=(segment_count, telescope_count))
        observed_star = star_counts/observation.segment_s
        observed_total = (
            star_counts + nsb_counts + dark_counts
        ) / observation.segment_s
        photon_relative_rms.append(float(np.median(
            1.0/np.sqrt(np.maximum(
                star_counts + nsb_counts + dark_counts, 1)))))

        nightly_gain = rng.normal(0.0, instrument.per_night_gain_rms,
                                  telescope_count)
        clock_ns = rng.normal(0.0, instrument.residual_timing_rms_ns,
                              telescope_count)
        delay_factor = np.interp(
            clock_ns[row_i]-clock_ns[row_j], timing_grid, timing_response,
            left=timing_response[0], right=timing_response[-1])
        pair_gain = ((1+static_gain[row_i])*(1+nightly_gain[row_i])
                     *(1+static_gain[row_j])*(1+nightly_gain[row_j]))
        timing_factors.extend(delay_factor)
        gain_factors.extend(pair_gain)

        star_i = observed_star[row_segment, row_i]
        star_j = observed_star[row_segment, row_j]
        total_i = observed_total[row_segment, row_i]
        total_j = observed_total[row_segment, row_j]
        unit_snr = (
            star_i*star_j/np.sqrt(total_i*total_j)
            * instrument.coherence_area_s
            * np.sqrt(2.0 * instrument.electronics_bandwidth_hz
                      * observation.segment_s)
            / charge_rms**2)
        sigma = 1.0/(unit_snr*electronics_efficiency*optical_efficiency)
        expected = (truth*pair_gain*delay_factor
                    + np.asarray([pair_zero[pair] for pair in pairs]))
        nightly_values.append(expected+rng.normal(0.0, sigma))
        nightly_variances.append(sigma**2)

    inverse_variance = 1.0/np.asarray(nightly_variances)
    summed_weight = inverse_variance.sum(axis=0)
    measured = (np.asarray(nightly_values)*inverse_variance).sum(0)/summed_weight
    sigma_stat = np.sqrt(1.0/summed_weight)
    sigma_shared = np.sqrt(
        (math.sqrt(2)*instrument.telescope_gain_calibration_rms*truth)**2
        + instrument.baseline_zero_point_rms**2)
    sigma_total = np.sqrt(sigma_stat**2+sigma_shared**2)
    frame["visibility2_measured"] = measured
    frame["sigma_visibility2_stat"] = sigma_stat
    frame["sigma_visibility2"] = sigma_total
    frame["baseline_zero_point_sigma"] = instrument.baseline_zero_point_rms

    endpoint = observation.hours_per_night*3600.*math.pi/SIDEREAL_DAY_S
    altitude = [math.asin(source_direction_enu(
        sign*endpoint, math.radians(observation.source_dec_deg),
        math.radians(observation.site_lat_deg))[2]) for sign in (-1, 1)]
    metadata = {
        "estimator": "analytic",
        "hours_per_night": observation.hours_per_night,
        "nights": observation.nights,
        "total_integration_hours": observation.hours_per_night*observation.nights,
        "integration_s_per_uv_measurement": (
            observation.segment_s*observation.nights),
        "visibility_subsamples_per_segment": (
            observation.visibility_subsamples_per_segment),
        "max_abs_segment_time_smearing": float(
            np.max(np.abs(truth-center_truth))),
        "source_ab_magnitude": source.ab_magnitude,
        "detected_star_rate_MHz": star_rate/1.0e6,
        "source_case": source_case,
        "nsb_multiplier": nsb_multiplier,
        "nsb_rate_MHz": nsb_rate/1.0e6,
        "dark_count_rate_MHz": dark_rate/1.0e6,
        "electronics_case": electronics_case,
        "electronics_correlation_efficiency": electronics_efficiency,
        "optical_timing_efficiency": optical_efficiency,
        "median_timing_attenuation": float(np.median(timing_factors)),
        "rms_pair_gain_error": float(np.std(gain_factors)),
        "median_photon_count_relative_rms": float(np.median(photon_relative_rms)),
        "minimum_altitude_deg": math.degrees(min(altitude)),
        "uv_measurements": len(frame),
        "sigma_visibility2_stat": float(np.median(sigma_stat)),
        "sigma_visibility2_total": float(np.median(sigma_total)),
        "instrument_parameter_source": instrument.parameter_source,
    }
    frame.attrs.update(metadata)
    return frame, metadata


def prepare_reconstruction_uv(measurements, cell_mlambda=120.0, legacy=False):
    """逆方差合并观测，并保留每个时间子采样在预测量中的相同权重。

    统计噪声独立；共享标定增益通过独立先验传入，不参与逆方差合并。
    其他已知相关统计误差应由UvData.covariance直接传入重建器。
    cell_mlambda=None不合并，仅用于大小受控的独立检查。
    """
    from sii_reconstruction import UvData
    data = measurements.copy()
    gain_sigma = 0.
    if "calibration_relative_sigma" in data:
        if data.calibration_relative_sigma.nunique() != 1:
            raise ValueError("one reconstruction must use one shared calibration prior")
        if "calibration_id" in data and data.calibration_id.nunique() != 1:
            raise ValueError("multiple independent calibrations need separate likelihood groups")
        gain_sigma = float(data.calibration_relative_sigma.iloc[0])
        if not np.isfinite(gain_sigma) or gain_sigma < 0:
            raise ValueError("calibration prior must be finite and nonnegative")
        if legacy and gain_sigma > 0:
            raise ValueError("legacy reconstruction does not support shared calibration uncertainty")
    required = ["u_lambda", "v_lambda", "visibility2_measured", "sigma_visibility2"]
    if (len(data) == 0 or not np.all(np.isfinite(data[required].to_numpy(float)))
            or np.any(data.sigma_visibility2 <= 0)):
        raise ValueError("UV coordinates, measurements and positive sigma must be finite")
    if cell_mlambda is None:
        data["ku"], data["kv"] = np.arange(len(data)), 0
    else:
        if cell_mlambda <= 0:
            raise ValueError("cell_mlambda must be positive or None")
        data["ku"] = np.round(data.u_lambda/1.0e6/cell_mlambda).astype(int)
        data["kv"] = np.round(data.v_lambda/1.0e6/cell_mlambda).astype(int)
    shared_zero = ("baseline_zero_point_sigma" in data
                   and np.any(data.baseline_zero_point_sigma.to_numpy() != 0))
    if shared_zero:
        if legacy:
            raise ValueError("legacy reconstruction does not support shared baseline offsets")
        if not {"sigma_visibility2_stat", "telescope_i", "telescope_j"}.issubset(data):
            raise ValueError("shared baseline offsets require statistical sigma and telescope IDs")
        if (not np.all(np.isfinite(data.baseline_zero_point_sigma))
                or np.any(data.baseline_zero_point_sigma < 0)
                or not np.all(np.isfinite(data.sigma_visibility2_stat))
                or np.any(data.sigma_visibility2_stat <= 0)):
            raise ValueError("invalid shared or independent baseline uncertainty")
    independent_sigma = data.sigma_visibility2_stat if shared_zero else data.sigma_visibility2
    data["inverse_variance"] = 1.0/independent_sigma**2
    data["weighted_value"] = (
        data.visibility2_measured*data.inverse_variance)
    grouped = data.groupby(["ku", "kv"], as_index=False).agg(
        u_lambda=("u_lambda", "mean"), v_lambda=("v_lambda", "mean"),
        weighted_value=("weighted_value", "sum"),
        inverse_variance=("inverse_variance", "sum"),
        multiplicity=("visibility2_measured", "size"))
    grouped["visibility2"] = grouped.weighted_value/grouped.inverse_variance
    grouped["sigma"] = np.sqrt(1.0/grouped.inverse_variance)
    sampling = None
    if legacy:
        grouped = pd.concat([grouped, pd.DataFrame([{
            "u_lambda": 0.0, "v_lambda": 0.0, "visibility2": 1.0,
            "sigma": 0.01, "multiplicity": 1}])], ignore_index=True)
        weight = 1.0/grouped.sigma.to_numpy(float)**2
        weight = np.clip(weight/np.mean(weight), 0.1, 10.0)
        weight /= weight.mean()
    else:
        weight = 1.0/grouped.sigma.to_numpy(float)**2
        groups = data.groupby(["ku", "kv"], sort=True).ngroup().to_numpy()
        if "uv_samples_u" in data and "uv_samples_v" in data:
            us = [np.asarray(x, float) for x in data.uv_samples_u]
            vs = [np.asarray(x, float) for x in data.uv_samples_v]
        else:
            us = [np.array([x]) for x in data.u_lambda]
            vs = [np.array([x]) for x in data.v_lambda]
        lengths = np.array([len(x) for x in us])
        if np.any(lengths == 0) or any(len(v) != len(u) for u, v in zip(us, vs)):
            raise ValueError("invalid per-measurement UV samples")
        time_weights = ([np.asarray(x, float) for x in data.uv_samples_weight]
                        if "uv_samples_weight" in data else
                        [np.full(n, 1.0/n) for n in lengths])
        if any(len(w) != n or not np.isclose(w.sum(), 1) or np.any(w < 0)
               for w, n in zip(time_weights, lengths)):
            raise ValueError("time weights must be nonnegative and sum to one")
        fractions = data.inverse_variance.to_numpy()/grouped.inverse_variance.to_numpy()[groups]
        sampling = (np.concatenate(us), np.concatenate(vs), np.repeat(groups, lengths),
                    np.concatenate([w*f for w, f in zip(time_weights, fractions)]))
    covariance = None
    if shared_zero:
        # A先合并独立测量；AB保留同一基线跨时间共享的零点，不按点数缩小。
        baseline_ids = pd.MultiIndex.from_frame(data[["telescope_i", "telescope_j"]])
        baseline_index, unique_baselines = pd.factorize(baseline_ids)
        design = np.zeros((len(grouped), len(unique_baselines)))
        np.add.at(design, (groups, baseline_index),
                  fractions*data.baseline_zero_point_sigma.to_numpy())
        covariance = np.diag(grouped.sigma.to_numpy()**2)+design@design.T
        grouped["sigma"] = np.sqrt(np.diag(covariance))
        weight = 1/grouped.sigma.to_numpy()**2
    uv = UvData(
        u_lambda=grouped.u_lambda.to_numpy(float),
        v_lambda=grouped.v_lambda.to_numpy(float),
        visibility_abs2=grouped.visibility2.to_numpy(float),
        sigma=grouped.sigma.to_numpy(float), weight=weight,
        multiplicity=grouped.multiplicity.to_numpy(int),
        input_rows=len(data), finite_rows=len(data),
        physical_violations=int(((grouped.visibility2 < 0)
                                 | (grouped.visibility2 > 1)).sum()), sampling=sampling,
        calibration_relative_sigma=gain_sigma, covariance=covariance)
    return uv


def reconstruct_uv(measurements, cell_mlambda=120.0, **kwargs):
    """独立重建入口：只读测量、误差和采样几何；默认使用绝对高斯似然。"""
    from sii_reconstruction import reconstruct_uv_data
    legacy = kwargs.get("likelihood") == "legacy"
    if legacy:
        kwargs.setdefault("smoothness", 0.020)
        kwargs.setdefault("huber_delta", 0.15)
    uv = prepare_reconstruction_uv(measurements, cell_mlambda, legacy=legacy)
    return reconstruct_uv_data(uv, **kwargs)


def run_sii_pipeline(
        layout, source=BinarySource(),
        observation: Observation = Observation(),
        instrument: Instrument = Instrument(), seed: int = 20260824,
        do_reconstruction: bool = True, reconstruction_kwargs=None,
        **simulation_kwargs) -> PipelineResult:
    """最小的一键入口：模型 → UVW → 测量 → 可选独立重建。"""
    uvw = generate_uvw(layout, observation, instrument)
    measurements, metadata = simulate_uv_observation(
        uvw, source, observation, instrument, seed=seed, **simulation_kwargs)
    reconstruction = None
    if do_reconstruction:
        defaults = {
            "grid_size": 28, "fov_mas": 0.70,
            "support_radius_mas": 0.32, "starts": 3,
            "max_iter": 1600, "smoothness": "cv",
            "huber_delta": 2.5, "seed": seed+1,
            "peak_minimum_separation_mas": 0.10,
        }
        defaults.update(reconstruction_kwargs or {})
        reconstruction = reconstruct_uv(measurements, **defaults)
    return PipelineResult(uvw, measurements, reconstruction, metadata)
