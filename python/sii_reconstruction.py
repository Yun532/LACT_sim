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
from dataclasses import dataclass
from pathlib import Path

import numpy as np

try:
    from scipy.optimize import minimize
except ImportError:  # pragma: no cover - reconstruction requires SciPy.
    minimize = None

try:
    import matplotlib

    # A reusable library must not replace an already active notebook backend.
    # Select Agg only for headless/script imports that have not loaded pyplot.
    if "matplotlib.pyplot" not in sys.modules:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:  # pragma: no cover - CSV reconstruction remains usable.
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


def read_uv_measurements(
        path,
        value_column,
        duplicate_tolerance_lambda=1e-3,
        sigma_column="auto"):
    """读取并合并重复 UV 样本；不会插值或虚构未观测的 UV 点。"""
    path = Path(path)
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
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
        floor = max(1e-6, 0.5 * float(np.percentile(finite_sigma, 20.0)))
        effective = np.where(np.isfinite(sigma) & (sigma > 0.0),
                             np.maximum(sigma, floor), np.nan)
        fallback = float(np.median(finite_sigma))
        effective = np.where(np.isfinite(effective), effective, fallback)
        weight = 1.0 / (effective * effective)
        weight /= np.mean(weight)
        weight = np.clip(weight, 0.1, 10.0)
        weight /= np.mean(weight)
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


def _initial_image(xx, yy, support, rng, start_index, fov_mas):
    if start_index == 0:
        image = np.exp(-(xx * xx + yy * yy) / (2.0 * (0.12 * fov_mas) ** 2))
    else:
        separation = rng.uniform(0.15, 0.45) * fov_mas
        angle = rng.uniform(0.0, math.pi)
        dx = 0.5 * separation * math.sin(angle)
        dy = 0.5 * separation * math.cos(angle)
        sigma = rng.uniform(0.04, 0.11) * fov_mas
        ratio = rng.uniform(0.25, 1.0)
        image = np.exp(-((xx - dx) ** 2 + (yy - dy) ** 2) / (2.0 * sigma ** 2))
        image += ratio * np.exp(
            -((xx + dx) ** 2 + (yy + dy) ** 2) / (2.0 * sigma ** 2)
        )
    supported = image[support] + 1e-8
    noise = rng.normal(0.0, 0.35, size=supported.size)
    return np.log(supported) + noise


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
    yy, xx = np.meshgrid(theta_mas, theta_mas, indexing="ij")
    order = np.argsort(image.ravel())[::-1]
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
        aligned = np.roll(candidate, (shift_y, shift_x), axis=(0, 1))
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
            })
    return best[1], best[2]


def reconstruct_uv_data(
        uv,
        grid_size=48,
        fov_mas=4.0,
        support_radius_mas=None,
        starts=4,
        max_iter=600,
        smoothness=0.2,
        huber_delta=0.18,
        seed=12345,
        peak_minimum_separation_mas=0.35,
        truth_sky_image=None):
    """在正值和有限支撑约束下，以多起点 L-BFGS-B 拟合 ``|V|²``。

    图像通过 softmax 保证非负并归一化；Huber 损失降低异常 UV 点的影响，
    ``smoothness`` 控制相邻像素正则。绝对平移和 180° 镜像是幅度干涉本身
    无法消除的固有简并，而不是优化器错误。
    """
    if minimize is None:
        raise RuntimeError("reconstruct-uv requires scipy.optimize")
    if starts <= 0 or max_iter <= 0:
        raise ValueError("starts and max_iter must be positive")
    if smoothness < 0.0 or huber_delta <= 0.0:
        raise ValueError("smoothness must be non-negative and huber_delta positive")
    if support_radius_mas is None:
        support_radius_mas = 0.47 * fov_mas
    theta, xx, yy, support = _image_grid(grid_size, fov_mas, support_radius_mas)
    fourier = _fourier_matrix(uv, xx, yy, support)
    observed = uv.visibility_abs2
    rng = np.random.default_rng(seed)
    trials = []

    def objective(z):
        supported_image = _softmax(z)
        field = fourier @ supported_image
        predicted = np.abs(field) ** 2
        residual = predicted - observed
        loss, influence = _huber(residual, huber_delta)
        weight_sum = float(np.sum(uv.weight))
        value = float(np.sum(uv.weight * loss) / weight_sum)
        gradient_x = 2.0 * np.real(
            fourier.T @ (uv.weight * influence * np.conj(field))
        ) / weight_sum
        full_image = np.zeros_like(xx)
        full_image[support] = supported_image
        regularization, regularization_gradient = _smoothness_value_gradient(
            full_image, smoothness)
        value += regularization
        gradient_x += regularization_gradient[support]
        gradient_z = supported_image * (
            gradient_x - np.dot(gradient_x, supported_image)
        )
        return value, gradient_z

    for start_index in range(starts):
        initial = _initial_image(xx, yy, support, rng, start_index, fov_mas)
        fit = minimize(
            objective,
            initial,
            method="L-BFGS-B",
            jac=True,
            options={
                "maxiter": max_iter,
                "ftol": 1e-12,
                "gtol": 1e-8,
                "maxcor": 20,
            },
        )
        trials.append(fit)
    selected_index = int(np.argmin([trial.fun for trial in trials]))
    selected = trials[selected_index]
    image = np.zeros_like(xx)
    image[support] = _softmax(selected.x)
    image, centroid_before, center_shift = _center_image(image)

    # Recalculate the visibility after choosing the translation gauge.
    supported_image = image[support]
    supported_image /= np.sum(supported_image)
    predicted = np.abs(fourier @ supported_image) ** 2
    residual = predicted - observed
    metrics = {
        "algorithm": "positive_softmax_multistart_lbfgsb",
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
        "grid_size": grid_size,
        "fov_mas": fov_mas,
        "pixel_scale_mas": float(theta[1] - theta[0]),
        "support_radius_mas": support_radius_mas,
        "smoothness": smoothness,
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
