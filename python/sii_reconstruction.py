#!/usr/bin/env python3
"""由离散采样的 SII ``|V(u,v)|²`` 做非参数图像重建。

优化器只读取观测量和不确定度，不读取模拟真值。真值参数只允许在拟合完成后
计算闭合指标，避免把答案泄漏给重建算法。
"""

from __future__ import annotations

import csv
import json
import math
import sys
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np
from scipy.signal import fftconvolve
from scipy.linalg import solve_triangular

try:
    from scipy.optimize import minimize
except ImportError:  # pragma: no cover - 重建需要SciPy。
    minimize = None

try:
    import matplotlib

    # 不替换Notebook已经启用的绘图后端；脚本尚未导入pyplot时才选择Agg。
    if "matplotlib.pyplot" not in sys.modules:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:  # pragma: no cover - CSV重建仍可使用。
    plt = None


MAS_TO_RAD = math.pi / (180.0 * 3600.0 * 1000.0)


@dataclass
class UvData:
    """合并重复采样后的 UV 数据、误差权重和输入质量统计。"""

    u_lambda: np.ndarray
    v_lambda: np.ndarray
    visibility_abs2: np.ndarray
    sigma: np.ndarray
    weight: np.ndarray
    multiplicity: np.ndarray
    input_rows: int
    finite_rows: int
    physical_violations: int
    # 每个实际子采样的(u,v,所属输出格,归一化平均权重)，不含源真值。
    sampling: tuple | None = None
    covariance: np.ndarray | None = None
    # 这是全体测量共享的标定先验，不随合并UV点数按平方根缩小。
    calibration_relative_sigma: float = 0.0


@dataclass
class ReconstructionResult:
    """重建图、脏自相关、可见度预测和诊断指标。"""

    theta_mas: np.ndarray
    image: np.ndarray
    dirty_autocorrelation: np.ndarray
    uv: UvData
    predicted_visibility_abs2: np.ndarray
    metrics: dict


def _parse_float(row, key):
    try:
        return float(row.get(key, "nan"))
    except (TypeError, ValueError):
        return math.nan


def write_uv_data(path, uv):
    """无损保存推断输入：积分节点、共享先验和协方差均随测量一起保存。"""
    payload = {name: getattr(uv, name) for name in (
        "u_lambda", "v_lambda", "visibility_abs2", "sigma", "weight",
        "multiplicity", "input_rows", "finite_rows", "physical_violations",
        "calibration_relative_sigma")}
    payload["schema_version"] = 1
    if uv.sampling is not None:
        payload.update(zip(("sample_u", "sample_v", "sample_group", "sample_weight"), uv.sampling))
    if uv.covariance is not None:
        payload["covariance"] = uv.covariance
    # 文件句柄避免numpy自动补后缀导致调用者读错文件。
    with Path(path).open("wb") as handle:
        np.savez_compressed(handle, **payload)


def read_uv_data(path):
    """读取完整UV数据；禁止pickle，且不以默认值补掉缺失的物理字段。"""
    with np.load(path, allow_pickle=False) as archive:
        if int(archive["schema_version"]) != 1:
            raise ValueError("unsupported UV data schema")
        values = {name: archive[name].copy() for name in (
            "u_lambda", "v_lambda", "visibility_abs2", "sigma", "weight", "multiplicity")}
        for name in ("input_rows", "finite_rows", "physical_violations"):
            values[name] = int(archive[name])
        values["calibration_relative_sigma"] = float(archive["calibration_relative_sigma"])
        sample_keys = ("sample_u", "sample_v", "sample_group", "sample_weight")
        if any(key in archive for key in sample_keys):
            values["sampling"] = tuple(archive[key].copy() for key in sample_keys)
        if "covariance" in archive:
            values["covariance"] = archive["covariance"].copy()
    uv = UvData(**values)
    n = len(uv.u_lambda)
    if (n == 0 or any(np.shape(values[key]) != (n,) for key in (
            "v_lambda", "visibility_abs2", "sigma", "weight", "multiplicity"))
            or any(not np.all(np.isfinite(values[key])) for key in (
                "u_lambda", "v_lambda", "visibility_abs2", "sigma", "weight"))
            or np.any(uv.sigma <= 0) or np.any(uv.weight <= 0)
            or not np.isfinite(uv.calibration_relative_sigma)
            or uv.calibration_relative_sigma < 0):
        raise ValueError("invalid complete UV measurements or uncertainties")
    if uv.sampling is not None:
        u, v, group, weight = uv.sampling
        if (any(x.ndim != 1 or len(x) != len(u) for x in (u, v, group, weight))
                or not all(np.all(np.isfinite(x)) for x in (u, v, group, weight))
                or not np.issubdtype(group.dtype, np.integer)
                or np.any(group < 0) or np.any(group >= n) or np.any(weight < 0)
                or not np.allclose(np.bincount(group, weights=weight, minlength=n), 1)):
            raise ValueError("invalid UV integration sampling")
    if uv.covariance is not None:
        cov = uv.covariance
        if (cov.shape != (n, n) or not np.all(np.isfinite(cov))
                or not np.allclose(cov, cov.T)
                or not np.allclose(np.diag(cov), uv.sigma**2)):
            raise ValueError("invalid UV covariance")
        np.linalg.cholesky(cov)
    return uv


def read_uv_measurements(
        path,
        value_column,
        duplicate_tolerance_lambda=1e-3,
        sigma_column="auto"):
    """读取并合并重复 UV 样本；不会插值或虚构未观测的 UV 点。"""
    path = Path(path)
    if path.suffix.lower() == ".npz":
        return read_uv_data(path)
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if rows and ({"uv_samples_u", "uv_samples_v", "calibration_relative_sigma",
                  "baseline_zero_point_sigma"} & set(rows[0])):
        raise ValueError("CSV contains integrated/shared-error observations; use write_uv_data/read_uv_data NPZ to preserve the likelihood")
    required = {"u_lambda", "v_lambda", value_column}
    if not rows or not required.issubset(rows[0]):
        missing = ", ".join(sorted(required - (set(rows[0]) if rows else set())))
        raise ValueError(f"{path} is missing required columns: {missing}")
    if duplicate_tolerance_lambda <= 0.0:
        raise ValueError("duplicate_tolerance_lambda must be positive")
    if sigma_column == "auto":
        candidate = f"{value_column}_sem"
        sigma_column = candidate if candidate in rows[0] else None
    elif sigma_column and sigma_column not in rows[0]:
        raise ValueError(f"{path} is missing uncertainty column: {sigma_column}")

    grouped = {}
    finite_rows = 0
    physical_violations = 0
    for row in rows:
        u = _parse_float(row, "u_lambda")
        v = _parse_float(row, "v_lambda")
        value = _parse_float(row, value_column)
        if not (math.isfinite(u) and math.isfinite(v) and math.isfinite(value)):
            continue
        finite_rows += 1
        if value < 0.0 or value > 1.0:
            physical_violations += 1
        key = (
            int(round(u / duplicate_tolerance_lambda)),
            int(round(v / duplicate_tolerance_lambda)),
        )
        entry = grouped.setdefault(key, [0.0, 0.0, 0.0, 0, 0.0, 0.0, 0])
        entry[0] += u
        entry[1] += v
        entry[2] += value
        entry[3] += 1
        sigma = _parse_float(row, sigma_column) if sigma_column else math.nan
        if sigma_column and not (math.isfinite(sigma) and sigma > 0.0):
            raise ValueError("every measured UV row must have a finite positive uncertainty")
        if math.isfinite(sigma) and sigma > 0.0:
            inverse_variance = 1.0 / (sigma * sigma)
            entry[4] += inverse_variance
            entry[5] += inverse_variance * value
            entry[6] += 1
    if not grouped:
        raise ValueError(f"{path} contains no finite uv measurements")

    values = []
    for (sum_u, sum_v, sum_value, count,
         inverse_variance, weighted_value, sigma_count) in grouped.values():
        if sigma_count > 0 and inverse_variance > 0.0:
            value = weighted_value / inverse_variance
            sigma = math.sqrt(1.0 / inverse_variance)
        else:
            value = sum_value / count
            sigma = math.nan
        values.append((sum_u / count, sum_v / count, value, count, sigma))
    values.sort(key=lambda item: (math.hypot(item[0], item[1]), item[0], item[1]))
    array = np.asarray(values, dtype=float)
    sigma = array[:, 4]
    finite_sigma = sigma[np.isfinite(sigma) & (sigma > 0.0)]
    if finite_sigma.size:
        if finite_sigma.size != len(sigma):
            raise ValueError("some UV measurements lack a finite positive uncertainty")
        weight = 1.0 / sigma**2
    else:
        weight = np.ones(len(array), dtype=float)
    return UvData(
        u_lambda=array[:, 0],
        v_lambda=array[:, 1],
        visibility_abs2=array[:, 2],
        sigma=sigma,
        weight=weight,
        multiplicity=array[:, 3].astype(int),
        input_rows=len(rows),
        finite_rows=finite_rows,
        physical_violations=physical_violations,
    )


def _image_grid(grid_size, fov_mas, support_radius_mas):
    if grid_size < 12:
        raise ValueError("grid_size must be at least 12")
    if fov_mas <= 0.0:
        raise ValueError("fov_mas must be positive")
    if not (0.0 < support_radius_mas <= fov_mas / math.sqrt(2.0)):
        raise ValueError("support_radius_mas must be positive and fit inside the image")
    theta = np.linspace(-0.5 * fov_mas, 0.5 * fov_mas, grid_size)
    xx, yy = np.meshgrid(theta, theta)
    support = xx * xx + yy * yy <= support_radius_mas * support_radius_mas
    return theta, xx, yy, support


def _fourier_matrix(uv, xx, yy, support):
    theta_x = xx[support] * MAS_TO_RAD
    theta_y = yy[support] * MAS_TO_RAD
    phase = -2j * math.pi * (
        uv.u_lambda[:, None] * theta_x[None, :]
        + uv.v_lambda[:, None] * theta_y[None, :]
    )
    return np.exp(phase)


def _smoothness_value_gradient(image, strength):
    if strength <= 0.0:
        return 0.0, np.zeros_like(image)
    dx = image[:, 1:] - image[:, :-1]
    dy = image[1:, :] - image[:-1, :]
    value = 0.5 * strength * (np.sum(dx * dx) + np.sum(dy * dy))
    gradient = np.zeros_like(image)
    gradient[:, :-1] -= strength * dx
    gradient[:, 1:] += strength * dx
    gradient[:-1, :] -= strength * dy
    gradient[1:, :] += strength * dy
    return float(value), gradient


def _huber(residual, delta):
    absolute = np.abs(residual)
    quadratic = absolute <= delta
    value = np.where(
        quadratic,
        0.5 * residual * residual,
        delta * (absolute - 0.5 * delta),
    )
    derivative = np.where(quadratic, residual, delta * np.sign(residual))
    return value, derivative


def _softmax(z):
    shifted = z - np.max(z)
    value = np.exp(shifted)
    return value / np.sum(value)


def power_sampling_kernel(uv, grid_size, fov_mas):
    """精确预计算平均功率算子，避免以平均UV坐标替代平均|V|²。

    Wiener–Khinchin: P(q)=sum_d A_I(d) cos(2πq·d)。A_I是像素图自相关。
    先平均cos核，再作用于A_I，与逐个子采样算|V|²后平均代数等价。
    这里没有补相位、插值UV、或对测量做自相关。
    """
    lag = np.arange(1-grid_size, grid_size)*fov_mas/(grid_size-1)*MAS_TO_RAD
    count = len(uv.u_lambda)
    if count*len(lag)**2*8 > 1_000_000_000:
        raise ValueError("exact power kernel exceeds 1 GB; use UV grouping or a smaller diagnostic subset")
    if uv.sampling is None:
        u, v = uv.u_lambda, uv.v_lambda
        group, fraction = np.arange(count), np.ones(count)
    else:
        u, v, group, fraction = uv.sampling
    u, v, fraction = map(lambda x: np.asarray(x, float), (u, v, fraction))
    original_group = np.asarray(group)
    group = np.asarray(group, int)
    if np.any(original_group != group):
        raise ValueError("sampling group indices must be integers")
    if not (u.shape == v.shape == group.shape == fraction.shape):
        raise ValueError("sampling arrays must have equal shapes")
    if (not np.all(np.isfinite(u+v+fraction)) or np.any(fraction < 0)
            or np.any(group < 0) or np.any(group >= count)):
        raise ValueError("invalid sampling coordinates or averaging weights")
    if not np.allclose(np.bincount(group, weights=fraction, minlength=count), 1.0):
        raise ValueError("sampling weights must sum to one for every measurement")
    order = np.argsort(group, kind="stable")
    boundaries = np.r_[0, np.cumsum(np.bincount(group, minlength=count))]
    kernel = np.empty((count, len(lag)**2))
    for index in range(count):
        selected = order[boundaries[index]:boundaries[index+1]]
        # 分离u、v指数，只需两个窄矩阵；积的实部就是余弦平均。
        ex = np.exp(-2j*np.pi*u[selected, None]*lag)
        ey = np.exp(-2j*np.pi*v[selected, None]*lag)
        kernel[index] = ((ey*fraction[selected, None]).T @ ex).real.ravel()
    return kernel


def power_from_image(kernel, image):
    """给定天图返回与实际测量使用完全相同平均方式的P。"""
    autocorrelation = fftconvolve(image, image[::-1, ::-1], mode="full")
    return kernel @ autocorrelation.ravel()


def _power_gradient(kernel, influence, image):
    lag_gradient = (kernel.T @ influence).reshape((2*len(image)-1,)*2)
    # 实图自相关为中心对称量，消除浮点级非对称后取精确伴随导数。
    lag_gradient = 0.5*(lag_gradient+lag_gradient[::-1, ::-1])
    return 2*fftconvolve(lag_gradient, image, mode="valid")


def profile_calibration_gain(uv, prediction):
    """对共享增益g作解析剖面拟合：最小化统计卡方与(g-1)^2/s_g^2之和。"""
    scale = uv.calibration_relative_sigma
    if scale <= 0:
        return 1., 0.
    base = replace(uv, calibration_relative_sigma=0.)
    _, precision_prediction = statistical_loss(base, prediction)
    information = float(prediction @ precision_prediction)+1/scale**2
    gain = (float(uv.visibility_abs2 @ precision_prediction)+1/scale**2)/information
    return gain, 1/np.sqrt(information)


def statistical_loss(uv, residual, likelihood="gaussian", huber_delta=2.5):
    """返回绝对统计损失及对P残差的导数；不归一化/截断逆方差。"""
    if likelihood == "noiseless":
        return 0.5*float(residual @ residual), residual
    if not np.isfinite(uv.calibration_relative_sigma) or uv.calibration_relative_sigma < 0:
        raise ValueError("calibration prior must be finite and nonnegative")
    if uv.calibration_relative_sigma > 0:
        if likelihood != "gaussian":
            raise ValueError("shared calibration profiling requires Gaussian likelihood")
        prediction = np.asarray(residual)+uv.visibility_abs2
        gain, _ = profile_calibration_gain(uv, prediction)
        value, gradient = statistical_loss(replace(uv, calibration_relative_sigma=0.),
                                           gain*prediction-uv.visibility_abs2)
        value += .5*((gain-1)/uv.calibration_relative_sigma)**2
        # 包络定理：在最优g处，对图像求导不需要额外的dg/dI项。
        return value, gain*gradient
    if not np.all(np.isfinite(uv.sigma) & (uv.sigma > 0)):
        raise ValueError("statistical reconstruction requires finite positive sigma; use likelihood='noiseless' explicitly for ideal data")
    if uv.covariance is None:
        standardized = residual/uv.sigma
        unwhiten_gradient = lambda x: x/uv.sigma
    else:
        covariance = np.asarray(uv.covariance, float)
        if (covariance.shape != (len(residual),)*2
                or not np.all(np.isfinite(covariance))
                or not np.allclose(covariance, covariance.T, rtol=1e-10,
                                   atol=1e-14*np.max(np.abs(covariance)))
                or not np.allclose(np.diag(covariance), uv.sigma**2, rtol=1e-8, atol=0)):
            raise ValueError("UV covariance must be symmetric and have diagonal sigma²")
        chol = np.linalg.cholesky(covariance)
        standardized = solve_triangular(chol, residual, lower=True)
        unwhiten_gradient = lambda x: solve_triangular(chol.T, x, lower=False)
    if likelihood == "gaussian":
        return 0.5*float(standardized @ standardized), unwhiten_gradient(standardized)
    if likelihood == "huber":
        values, influence = _huber(standardized, huber_delta)
        return float(values.sum()), unwhiten_gradient(influence)
    raise ValueError("likelihood must be gaussian, huber, noiseless or legacy")


def _initial_image(xx, yy, support, rng, start_index, fov_mas):
    """使用集中、弥散和随机平滑初值；不把双星形态预先编码到所有随机起点。"""
    from scipy.ndimage import gaussian_filter
    if start_index % 3 == 0:
        image = np.exp(-(xx**2+yy**2)/(2*(.12*fov_mas)**2))
    elif start_index % 3 == 1:
        image = np.ones_like(xx)
    else:
        field = gaussian_filter(rng.normal(size=xx.shape), max(1., len(xx)/12.))
        image = np.exp(field/max(float(field.std()), 1e-12))
    return np.log(image[support]+1e-8)+rng.normal(0., .15, support.sum())



def _dirty_autocorrelation(uv, xx, yy):
    theta_x = xx.ravel() * MAS_TO_RAD
    theta_y = yy.ravel() * MAS_TO_RAD
    phase = 2.0 * math.pi * (
        uv.u_lambda[:, None] * theta_x[None, :]
        + uv.v_lambda[:, None] * theta_y[None, :]
    )
    dirty = np.mean(uv.visibility_abs2[:, None] * np.cos(phase), axis=0)
    return dirty.reshape(xx.shape)


def _shift_without_wrap(image, shift_y, shift_x):
    output = np.zeros_like(image)
    src_y0 = max(0, -shift_y)
    src_y1 = min(image.shape[0], image.shape[0] - shift_y)
    src_x0 = max(0, -shift_x)
    src_x1 = min(image.shape[1], image.shape[1] - shift_x)
    dst_y0 = src_y0 + shift_y
    dst_y1 = src_y1 + shift_y
    dst_x0 = src_x0 + shift_x
    dst_x1 = src_x1 + shift_x
    if src_y1 > src_y0 and src_x1 > src_x0:
        output[dst_y0:dst_y1, dst_x0:dst_x1] = image[src_y0:src_y1, src_x0:src_x1]
    return output


def _center_image(image):
    total = float(np.sum(image))
    indices = np.arange(image.shape[0], dtype=float)
    cy = float(np.sum(image * indices[:, None]) / total)
    cx = float(np.sum(image * indices[None, :]) / total)
    target = 0.5 * (image.shape[0] - 1)
    shift_y = int(round(target - cy))
    shift_x = int(round(target - cx))
    shifted = _shift_without_wrap(image, shift_y, shift_x)
    if np.sum(shifted) > 0.0:
        shifted /= np.sum(shifted)
    return shifted, (cy, cx), (shift_y, shift_x)


def _peak_diagnostic(image, theta_mas, minimum_separation_mas):
    """只报告实际局部极大值；单峰图不强行挑出第二个亮像素冒充双星。"""
    from scipy.ndimage import maximum_filter, label
    yy, xx = np.meshgrid(theta_mas, theta_mas, indexing="ij")
    maxima = (image == maximum_filter(image, size=3)) & (image > .05*image.max())
    regions, count = label(maxima)
    candidates = []
    for region in range(1, count+1):
        indices = np.flatnonzero(regions.ravel() == region)
        candidates.append(indices[np.argmax(image.ravel()[indices])])
    order = sorted(candidates, key=lambda index: image.ravel()[index], reverse=True)
    peaks = []
    for flat_index in order:
        row, col = np.unravel_index(flat_index, image.shape)
        x = float(xx[row, col])
        y = float(yy[row, col])
        if all(math.hypot(x - px, y - py) >= minimum_separation_mas
               for px, py, _ in peaks):
            peaks.append((x, y, float(image[row, col])))
        if len(peaks) == 2:
            break
    result = {"peaks": [
        {"theta_x_mas": p[0], "theta_y_mas": p[1], "brightness": p[2]}
        for p in peaks
    ]}
    if len(peaks) == 2:
        dx = peaks[1][0] - peaks[0][0]
        dy = peaks[1][1] - peaks[0][1]
        result.update({
            "two_peak_separation_mas": math.hypot(dx, dy),
            "two_peak_position_angle_deg": math.degrees(math.atan2(dx, dy)) % 180.0,
            "two_peak_brightness_ratio": min(peaks[0][2], peaks[1][2])
            / max(peaks[0][2], peaks[1][2]),
        })
    return result


def _truth_grid(path, theta_mas):
    with Path(path).open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    xs = sorted({_parse_float(row, "theta_x_mas") for row in rows})
    ys = sorted({_parse_float(row, "theta_y_mas") for row in rows})
    if not xs or not ys or len(rows) != len(xs) * len(ys):
        raise ValueError("truth sky image must be a complete rectangular grid")
    x_index = {value: index for index, value in enumerate(xs)}
    y_index = {value: index for index, value in enumerate(ys)}
    source = np.zeros((len(ys), len(xs)), dtype=float)
    for row in rows:
        x = _parse_float(row, "theta_x_mas")
        y = _parse_float(row, "theta_y_mas")
        source[y_index[y], x_index[x]] = max(0.0, _parse_float(row, "intensity"))
    along_x = np.vstack([
        np.interp(theta_mas, xs, source_row, left=0.0, right=0.0)
        for source_row in source
    ])
    target = np.vstack([
        np.interp(theta_mas, ys, along_x[:, column], left=0.0, right=0.0)
        for column in range(along_x.shape[1])
    ]).T
    if np.sum(target) <= 0.0:
        raise ValueError("truth sky image has zero flux in the reconstruction field")
    return target / np.sum(target)


def _best_truth_alignment(image, truth):
    best = None
    for mirrored, candidate in ((False, image), (True, image[::-1, ::-1])):
        correlation = np.fft.ifft2(
            np.fft.fft2(truth) * np.conj(np.fft.fft2(candidate))
        ).real
        row, col = np.unravel_index(np.argmax(correlation), correlation.shape)
        shift_y = int(row if row <= image.shape[0] // 2 else row - image.shape[0])
        shift_x = int(col if col <= image.shape[1] // 2 else col - image.shape[1])
        # 配准仅用于拟合后评价；禁止周期卷绕把视场边缘通量搬到另一侧。
        aligned = _shift_without_wrap(candidate, shift_y, shift_x)
        a = aligned.ravel()
        b = truth.ravel()
        coefficient = float(np.corrcoef(a, b)[0, 1])
        nrmse = float(np.sqrt(np.mean((a - b) ** 2)) / (np.sqrt(np.mean(b ** 2)) + 1e-15))
        score = coefficient
        if best is None or score > best[0]:
            best = (score, aligned, {
                "truth_correlation": coefficient,
                "truth_nrmse": nrmse,
                "truth_uses_180_degree_mirror": mirrored,
                "truth_alignment_shift_pixels": [shift_y, shift_x],
                "truth_alignment_retained_flux": float(aligned.sum()/candidate.sum()),
            })
    return best[1], best[2]


def reconstruct_uv_data(
        uv,
        grid_size=48,
        fov_mas=4.0,
        support_radius_mas=None,
        starts=4,
        max_iter=8000,
        smoothness="cv",
        huber_delta=2.5,
        seed=12345,
        peak_minimum_separation_mas=0.35,
        truth_sky_image=None,
        likelihood="gaussian",
        smoothness_candidates=(0.0, 1e-5, 1e-3, 1e-1),
        cv_seed=314159,
        optimizer_ftol=None,
        parameterization=None,
        _power_kernel=None):
    """在正值和有限支撑约束下，以多起点 L-BFGS-B 拟合 ``|V|²``。

    默认直接优化有界非负像素再归一化，避免softmax使暗像素梯度过小；
    parameterization="softmax"保留旧参数化，legacy默认仍用它。两者物理目标相同。
    默认使用绝对Gaussian误差，
    可选Huber对标准化残差降低异常点影响。smoothness="cv"只用留出测量选先验，
    数值强度对应视场归一化坐标上的亮度密度梯度；旧定义仅在legacy模式保留。
    绝对平移和180°中心反演是幅度干涉的固有简并。
    """
    if minimize is None:
        raise RuntimeError("reconstruct-uv requires scipy.optimize")
    if starts <= 0 or max_iter <= 0:
        raise ValueError("starts and max_iter must be positive")
    if parameterization is None:
        parameterization = "softmax" if likelihood == "legacy" else "flux"
    if parameterization not in ("softmax", "flux"):
        raise ValueError("parameterization must be softmax or flux")
    if (len(uv.u_lambda) == 0 or any(len(x) != len(uv.u_lambda) for x in
            (uv.v_lambda, uv.visibility_abs2, uv.sigma, uv.weight, uv.multiplicity))
            or not np.all(np.isfinite(uv.visibility_abs2))):
        raise ValueError("UV arrays must be nonempty, equally sized, with finite measurements")
    if (smoothness != "cv" and smoothness < 0.0) or huber_delta <= 0.0:
        raise ValueError("smoothness must be non-negative and huber_delta positive")
    if support_radius_mas is None:
        support_radius_mas = 0.47 * fov_mas
    theta, xx, yy, support = _image_grid(grid_size, fov_mas, support_radius_mas)
    legacy = likelihood == "legacy"
    # 新目标是卡方总和，不是旧平均损失。保留可调停止容差，并在报告中做严格容差对照。
    ftol = (1e-12 if parameterization == "flux" or likelihood in ("legacy", "noiseless")
            else 1e-9) if optimizer_ftol is None else optimizer_ftol
    if not np.isfinite(ftol) or ftol <= 0:
        raise ValueError("optimizer_ftol must be finite and positive")
    if legacy:
        if smoothness == "cv":
            raise ValueError("legacy comparison requires an explicit old smoothness")
        fourier = _fourier_matrix(uv, xx, yy, support)
        kernel = None
    else:
        kernel = (power_sampling_kernel(uv, grid_size, fov_mas)
                  if _power_kernel is None else _power_kernel)
        statistical_loss(uv, np.zeros(len(uv.u_lambda)), likelihood, huber_delta)
    cv_results = []
    if smoothness == "cv":
        if uv.covariance is not None or likelihood == "noiseless":
            raise ValueError("CV currently requires independent noisy UV measurements")
        if len(uv.u_lambda) < 10:
            raise ValueError("CV needs at least 10 UV measurements; specify smoothness explicitly")
        # 一次固定80/20留出，只读测量；不读真图，也不按图像相似度选参数。
        order = np.random.default_rng(cv_seed).permutation(len(uv.u_lambda))
        validation, training = order[:max(2, len(order)//5)], order[max(2, len(order)//5):]
        def subset(indices):
            return replace(uv, u_lambda=uv.u_lambda[indices], v_lambda=uv.v_lambda[indices],
                visibility_abs2=uv.visibility_abs2[indices], sigma=uv.sigma[indices],
                weight=uv.weight[indices], multiplicity=uv.multiplicity[indices], sampling=None)
        train_uv, validation_uv = subset(training), subset(validation)
        for alpha in smoothness_candidates:
            fit = reconstruct_uv_data(train_uv, grid_size=grid_size, fov_mas=fov_mas,
                support_radius_mas=support_radius_mas, starts=min(starts, 2),
                max_iter=max_iter, smoothness=float(alpha), huber_delta=huber_delta,
                seed=seed, likelihood=likelihood, optimizer_ftol=ftol,
                parameterization=parameterization, _power_kernel=kernel[training])
            validation_prediction = power_from_image(kernel[validation], fit.image)
            train_prediction = power_from_image(kernel[training], fit.image)
            gain, gain_sigma = profile_calibration_gain(train_uv, train_prediction)
            residual = gain*validation_prediction-validation_uv.visibility_abs2
            # 验证集只接收训练集的增益后验；不允许验证集自行重新标定。
            standardized = residual/validation_uv.sigma
            factor = gain_sigma*validation_prediction/validation_uv.sigma
            denominator = 1.+factor @ factor
            score = .5*(standardized @ standardized
                         -(standardized @ factor)**2/denominator+np.log(denominator))
            cv_results.append({"smoothness": float(alpha), "validation_chi2_per_point":
                2*score/len(validation), "training_converged": fit.metrics["optimizer_success"]})
        smoothness = min(cv_results, key=lambda row: row["validation_chi2_per_point"])["smoothness"]
    observed = uv.visibility_abs2
    rng = np.random.default_rng(seed)
    trials = []

    def image_objective(supported_image):
        """先在真实像素流量上求导，再通过参数化传回优化器。"""
        full_image = np.zeros_like(xx)
        full_image[support] = supported_image
        if legacy:
            field = fourier @ supported_image
            predicted = np.abs(field) ** 2
            residual = predicted - observed
            loss, influence = _huber(residual, huber_delta)
            weight_sum = float(np.sum(uv.weight))
            value = float(np.sum(uv.weight * loss) / weight_sum)
            gradient_x = 2.0 * np.real(
                fourier.T @ (uv.weight * influence * np.conj(field))) / weight_sum
        else:
            predicted = power_from_image(kernel, full_image)
            value, influence = statistical_loss(uv, predicted-observed, likelihood, huber_delta)
            gradient_x = _power_gradient(kernel, influence, full_image)[support]
        # 在视场归一化坐标上惩罚亮度密度梯度；(N-1)^4补偿像素流量随网格变化。
        strength = smoothness if legacy else smoothness*(grid_size-1)**4
        regularization, regularization_gradient = _smoothness_value_gradient(
            full_image, strength)
        value += regularization
        gradient_x += regularization_gradient[support]
        return value, gradient_x

    def objective(z):
        if parameterization == "flux" and z.sum() <= 0:
            # 全零向量没有定义归一化图像；拒绝这一步，不虚构预测量。
            return np.inf, -np.ones_like(z)
        pixels = _softmax(z) if parameterization == "softmax" else z/z.sum()
        value, gradient = image_objective(pixels)
        tangent = gradient-np.dot(gradient, pixels)
        return value, tangent*(pixels if parameterization == "softmax" else 1/z.sum())

    for start_index in range(starts):
        initial = _initial_image(xx, yy, support, rng, start_index, fov_mas)
        if parameterization == "flux":
            initial = _softmax(initial)
        fit = minimize(
            objective,
            initial,
            method="L-BFGS-B",
            jac=True,
            bounds=[(0.0, 1.0)]*len(initial) if parameterization == "flux" else None,
            options={
                "maxiter": max_iter,
                "ftol": ftol,
                "gtol": 1e-8,
                "maxcor": 20,
            },
        )
        trials.append(fit)
    selected_index = int(np.argmin([trial.fun for trial in trials]))
    selected = trials[selected_index]
    image = np.zeros_like(xx)
    image[support] = (_softmax(selected.x) if parameterization == "softmax"
                      else selected.x/selected.x.sum())
    _, pixel_gradient = image_objective(image[support])
    # 单纯形上的一阶驻点检查：把少量流量移到梯度最小像素，能否降低目标。
    # 这是局部必要条件，不是非凸问题的全局最优证明。
    stationarity_gap = float(image[support] @ pixel_gradient-pixel_gradient.min())
    centered, centroid_before, center_shift = _center_image(image)
    # 中心平移若会裁掉通量或移出支撑，不改变优化后的解。禁止裁剪后重新归一化改变P。
    if legacy:
        image = centered
    else:
        candidate = _shift_without_wrap(image, *center_shift)
        if (np.isclose(candidate.sum(), image.sum(), rtol=0, atol=1e-13)
                and np.sum(candidate[~support]) < 1e-13):
            image = candidate
        else:
            center_shift = (0, 0)

    # 选定允许的平移后重新计算可见度，保证保存的图像与预测一致。
    supported_image = image[support]
    supported_image /= np.sum(supported_image)
    predicted = (np.abs(fourier @ supported_image) ** 2 if legacy
                 else power_from_image(kernel, image))
    residual = predicted - observed
    calibration_gain, calibration_gain_sigma = profile_calibration_gain(uv, predicted)
    metrics = {
        "algorithm": f"positive_{parameterization}_multistart_lbfgsb",
        "parameterization": parameterization,
        "calibration_gain": float(calibration_gain),
        "calibration_gain_posterior_sigma": float(calibration_gain_sigma),
        "chi2_statistical": (float("nan") if legacy or likelihood == "noiseless" else
            2*statistical_loss(replace(uv, calibration_relative_sigma=0.),
                              calibration_gain*predicted-observed, "gaussian")[0]),
        "chi2_calibration_prior": ((calibration_gain-1)/uv.calibration_relative_sigma)**2
            if uv.calibration_relative_sigma > 0 else 0.,
        "simplex_stationarity_gap": stationarity_gap,
        "stationarity_gap_per_objective": stationarity_gap/max(1., abs(float(selected.fun))),
        "stationarity_passed": stationarity_gap <= 1e-4*max(1., abs(float(selected.fun))),
        "objective": float(selected.fun),
        "fit_rmse": float(np.sqrt(np.mean(residual * residual))),
        "weighted_fit_rmse": float(np.sqrt(
            np.sum(uv.weight * residual * residual) / np.sum(uv.weight)
        )),
        "fit_mae": float(np.mean(np.abs(residual))),
        "selected_start": selected_index,
        "starts": starts,
        "optimizer_success": bool(selected.success),
        "optimizer_status": int(selected.status),
        "optimizer_message": str(selected.message),
        "optimizer_iterations": int(selected.nit),
        "optimizer_ftol": float(ftol),
        "grid_size": grid_size,
        "fov_mas": fov_mas,
        "pixel_scale_mas": float(theta[1] - theta[0]),
        "support_radius_mas": support_radius_mas,
        "smoothness": smoothness,
        "likelihood": likelihood,
        "smoothness_selection": cv_results,
        "per_start_objective": [float(t.fun) for t in trials],
        "per_start_success": [bool(t.success) for t in trials],
        "per_start_iterations": [int(t.nit) for t in trials],
        "chi2": (float("nan") if legacy or likelihood == "noiseless" else
                 2*statistical_loss(uv, residual, "gaussian")[0]),
        "forward_model": "representative_uv" if legacy else "averaged_power_exact",
        "huber_delta": huber_delta,
        "seed": seed,
        "input_rows": uv.input_rows,
        "finite_rows": uv.finite_rows,
        "unique_uv_samples": len(observed),
        "physical_visibility_violations": uv.physical_violations,
        "uses_uncertainty_weights": bool(np.any(np.isfinite(uv.sigma))),
        "centroid_before_gauge_fix_pixels": [centroid_before[0], centroid_before[1]],
        "centering_shift_pixels": [center_shift[0], center_shift[1]],
        "ambiguities": [
            "absolute image translation is not observable from visibility magnitude",
            "a 180-degree image mirror has the same visibility magnitude",
        ],
    }
    metrics.update(_peak_diagnostic(image, theta, peak_minimum_separation_mas))
    if truth_sky_image:
        truth = _truth_grid(truth_sky_image, theta)
        aligned, truth_metrics = _best_truth_alignment(image, truth)
        metrics.update(truth_metrics)
        metrics["truth_sky_image"] = str(Path(truth_sky_image))
        metrics["truth_image"] = truth
        metrics["truth_aligned_reconstruction"] = aligned
    return ReconstructionResult(
        theta_mas=theta,
        image=image,
        dirty_autocorrelation=_dirty_autocorrelation(uv, xx, yy),
        uv=uv,
        predicted_visibility_abs2=predicted,
        metrics=metrics,
    )


def _serializable_metrics(metrics):
    return {key: value for key, value in metrics.items()
            if not isinstance(value, np.ndarray)}


def _write_image_csv(path, result):
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["theta_x_mas", "theta_y_mas", "brightness"])
        for row, theta_y in enumerate(result.theta_mas):
            for column, theta_x in enumerate(result.theta_mas):
                writer.writerow([theta_x, theta_y, result.image[row, column]])


def _write_visibility_fit(path, result):
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "u_lambda", "v_lambda", "observed_visibility_abs2",
            "predicted_visibility_abs2", "residual", "sigma", "weight",
            "multiplicity",
        ])
        for index in range(len(result.uv.visibility_abs2)):
            observed = result.uv.visibility_abs2[index]
            predicted = result.predicted_visibility_abs2[index]
            writer.writerow([
                result.uv.u_lambda[index], result.uv.v_lambda[index], observed,
                predicted, predicted - observed, result.uv.sigma[index],
                result.uv.weight[index], result.uv.multiplicity[index],
            ])


def _plot_reconstruction(path, result, value_column):
    if plt is None:
        return False
    extent = [
        result.theta_mas[0], result.theta_mas[-1],
        result.theta_mas[0], result.theta_mas[-1],
    ]
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 9.0))
    scatter = axes[0, 0].scatter(
        result.uv.u_lambda / 1e6,
        result.uv.v_lambda / 1e6,
        c=result.uv.visibility_abs2,
        cmap="viridis",
        s=24,
    )
    axes[0, 0].scatter(
        -result.uv.u_lambda / 1e6,
        -result.uv.v_lambda / 1e6,
        c=result.uv.visibility_abs2,
        cmap="viridis",
        s=12,
        alpha=0.35,
    )
    axes[0, 0].set_title(f"Sampled {value_column}")
    axes[0, 0].set_xlabel("u [Mlambda]")
    axes[0, 0].set_ylabel("v [Mlambda]")
    axes[0, 0].set_aspect("equal", adjustable="datalim")
    fig.colorbar(scatter, ax=axes[0, 0], label="visibility squared")

    dirty_limit = max(1e-12, float(np.max(np.abs(result.dirty_autocorrelation))))
    dirty_plot = axes[0, 1].imshow(
        result.dirty_autocorrelation,
        origin="lower",
        extent=extent,
        cmap="coolwarm",
        vmin=-dirty_limit,
        vmax=dirty_limit,
    )
    axes[0, 1].set_title("Dirty autocorrelation from sampled power")
    axes[0, 1].set_xlabel("delta theta_x [mas]")
    axes[0, 1].set_ylabel("delta theta_y [mas]")
    fig.colorbar(dirty_plot, ax=axes[0, 1])

    image_plot = axes[1, 0].imshow(
        result.image,
        origin="lower",
        extent=extent,
        cmap="magma",
    )
    axes[1, 0].set_title("Positive regularized reconstruction")
    axes[1, 0].set_xlabel("theta_x [mas]")
    axes[1, 0].set_ylabel("theta_y [mas]")
    fig.colorbar(image_plot, ax=axes[1, 0], label="normalized brightness / pixel")

    observed = result.uv.visibility_abs2
    predicted = result.predicted_visibility_abs2
    lo = min(float(np.min(observed)), float(np.min(predicted)), 0.0)
    hi = max(float(np.max(observed)), float(np.max(predicted)), 1.0)
    axes[1, 1].scatter(observed, predicted, s=19, alpha=0.75)
    axes[1, 1].plot([lo, hi], [lo, hi], color="0.25", linewidth=1.0)
    axes[1, 1].set_xlim(lo, hi)
    axes[1, 1].set_ylim(lo, hi)
    axes[1, 1].set_xlabel("observed visibility squared")
    axes[1, 1].set_ylabel("reconstructed visibility squared")
    axes[1, 1].set_title(f"Forward closure, RMSE={result.metrics['fit_rmse']:.4g}")
    axes[1, 1].grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)
    return True


def _plot_truth_comparison(path, result):
    truth = result.metrics.get("truth_image")
    aligned = result.metrics.get("truth_aligned_reconstruction")
    if plt is None or truth is None or aligned is None:
        return False
    extent = [
        result.theta_mas[0], result.theta_mas[-1],
        result.theta_mas[0], result.theta_mas[-1],
    ]
    fig, axes = plt.subplots(1, 3, figsize=(12.0, 3.8))
    vmax = max(float(np.max(truth)), float(np.max(aligned)))
    for axis, image, title in (
            (axes[0], truth, "Simulation truth (validation only)"),
            (axes[1], aligned, "Reconstruction, ambiguity-aligned"),
            (axes[2], aligned - truth, "Reconstruction - truth")):
        if axis is axes[2]:
            limit = max(1e-15, float(np.max(np.abs(image))))
            artist = axis.imshow(image, origin="lower", extent=extent,
                                 cmap="coolwarm", vmin=-limit, vmax=limit)
        else:
            artist = axis.imshow(image, origin="lower", extent=extent,
                                 cmap="magma", vmin=0.0, vmax=vmax)
        axis.set_title(title)
        axis.set_xlabel("theta_x [mas]")
        axis.set_ylabel("theta_y [mas]")
        fig.colorbar(artist, ax=axis)
    fig.suptitle(
        f"truth correlation={result.metrics['truth_correlation']:.4f}, "
        f"NRMSE={result.metrics['truth_nrmse']:.4f}"
    )
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)
    return True


def reconstruct_uv_csv(
        input_csv,
        output_dir,
        value_column="measured_visibility_abs2",
        duplicate_tolerance_lambda=1e-3,
        sigma_column="auto",
        **kwargs):
    """从标准 UV CSV 重建，并写出图像、拟合表、指标 JSON 和诊断图。"""
    uv = read_uv_measurements(
        input_csv, value_column, duplicate_tolerance_lambda, sigma_column)
    result = reconstruct_uv_data(uv, **kwargs)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    _write_image_csv(output_dir / "reconstruction.csv", result)
    _write_visibility_fit(output_dir / "visibility_fit.csv", result)
    metrics = _serializable_metrics(result.metrics)
    metrics.update({
        "input_csv": str(Path(input_csv)),
        "value_column": value_column,
        "sigma_column": sigma_column,
        "truth_used_by_optimizer": False,
    })
    (output_dir / "reconstruction_metrics.json").write_text(
        json.dumps(metrics, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    _plot_reconstruction(output_dir / "reconstruction.png", result, value_column)
    _plot_truth_comparison(output_dir / "truth_comparison.png", result)
    return result


def reconstruction_self_test():
    """用确定性双星数据执行傅里叶功率重建的最小闭合检查。"""
    rng = np.random.default_rng(90210)
    angles = rng.uniform(0.0, 2.0 * math.pi, 120)
    radii = np.sqrt(rng.uniform(0.0, 1.0, 120)) * 2.0e8
    u = radii * np.cos(angles)
    v = radii * np.sin(angles)
    separation_mas = 1.0
    position_angle = math.radians(32.0)
    dx = separation_mas * math.sin(position_angle) * MAS_TO_RAD
    dy = separation_mas * math.cos(position_angle) * MAS_TO_RAD
    ratio = 0.55
    field = (1.0 + ratio * np.exp(-2j * math.pi * (u * dx + v * dy))) / (1.0 + ratio)
    uv = UvData(
        u_lambda=np.concatenate(([0.0], u)),
        v_lambda=np.concatenate(([0.0], v)),
        visibility_abs2=np.concatenate(([1.0], np.abs(field) ** 2)),
        sigma=np.full(121, np.nan),
        weight=np.ones(121, dtype=float),
        multiplicity=np.ones(121, dtype=int),
        input_rows=121,
        finite_rows=121,
        physical_violations=0,
    )
    result = reconstruct_uv_data(
        uv,
        grid_size=32,
        fov_mas=3.0,
        support_radius_mas=1.4,
        starts=5,
        max_iter=800,
        smoothness=0.02,
        likelihood="legacy",
        huber_delta=0.5,
        seed=321,
        peak_minimum_separation_mas=0.35,
    )
    if result.metrics["fit_rmse"] > 0.05:
        raise RuntimeError(
            f"reconstruction self-test visibility RMSE is too large: "
            f"{result.metrics['fit_rmse']:.6g}"
        )
    recovered = result.metrics.get("two_peak_separation_mas", math.nan)
    if not math.isfinite(recovered) or abs(recovered - separation_mas) > 0.35:
        raise RuntimeError(
            f"reconstruction self-test separation mismatch: {recovered:.6g} mas"
        )
    return result.metrics
