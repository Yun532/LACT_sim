#!/usr/bin/env python3
"""Build the paper-style LACT stellar intensity-interferometry notebook."""

from pathlib import Path
import textwrap

import nbformat as nbf


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "notebooks" / "lact_sii_paper_simulation.ipynb"


def md(source):
    return nbf.v4.new_markdown_cell(textwrap.dedent(source).strip())


def code(source):
    return nbf.v4.new_code_cell(textwrap.dedent(source).strip())


nb = nbf.v4.new_notebook()
nb["metadata"] = {
    "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
    "language_info": {"name": "python", "version": "3"},
}

nb["cells"] = [
    md(r"""
    # LACT 恒星强度干涉：从天图、(uv) 测量到无相位成像

    **目的。** 这是一条可执行、可审计的论文级闭环：输入天体亮度分布，生成理论复可见度和完整 (uv) 功率面；用 LACT 本地 ENU 阵列与地球自转生成每一对望远镜的 (u,v,w) 和几何时延；按 Hanbury Brown 灵敏度关系生成带不确定度的 (|V|^2)；最后仅用模拟观测值做正值、有限支撑、平滑正则的无相位图像重建。

    ## tl;dr

    - 统一实现保留旧 v1 的多源模型、全阵列相干矩阵思想和不确定度加权重建，并吸收旧 v2 对飞秒光学相干、纳秒电子学及小时级积分必须分层处理的修正。
    - notebook 同时计算仓库 7 镜验证阵列和历史 `layout_0803` 生产输入的 **TEL.1–TEL.32 实际模拟坐标**；不再使用理想化 grouped layout。
    - 默认科学案例是 (m_{AB}=2)、分离 (0.20\,\mathrm{mas}) 的不等亮双星；重建器看不到真值相位，真值只在优化结束后用于平移/180°镜像歧义下的验证。
    - 星等极限和角分辨率是**假设驱动的预测**，不是 LACT 已实测性能；结果表明确列出通光效率、NSB、电子带宽、观测时间和光谱复用数。
    """),
    md(r"""
    ## 1. Context & Methods

    对混沌（热）光，Siegert 关系给出

    \[
    g^{(2)}_{ij}(\tau)-1=\beta\,|\gamma_{ij}|^2\,|g^{(1)}(\tau)|^2,
    \]

    而 van Cittert–Zernike 定理将复相干度 \(\gamma(u,v)\) 与归一化源亮度天图的傅里叶变换相连。强度干涉直接测得的是 \(|\gamma|^2\)，不测傅里叶相位，因此图像存在绝对平移与 180°反演歧义。

    长曝光、shot-noise 主导时，本 notebook 使用与 Hanbury Brown、Le Bohec–Holder、Rou 等工作一致的二阶相关近似：

    \[
    \mathrm{SNR}_{ij}=\frac{r_{\star,i}r_{\star,j}}
    {\sqrt{r_{\mathrm{tot},i}r_{\mathrm{tot},j}}\,\Delta\nu}
    |\gamma_{ij}|^2\sqrt{\frac{\Delta f\,T}{2}}\,\frac{1}{F_{\rm EN}},
    \]

    等口径时化为代码中的 `unit_visibility_snr × |V|²`。背景进入总 shot noise，但不产生跨镜恒星相关。观测模拟不裁剪负的 \(|V|^2\) 估计，避免低 SNR 下的正偏差。

    ### 关键参考

    1. Le Bohec & Holder (2006), *Optical Intensity Interferometry with Atmospheric Cherenkov Telescope Arrays*, DOI: [10.1086/507030](https://doi.org/10.1086/507030).
    2. Dravins et al. (2012), *Stellar Intensity Interferometry: Prospects for sub-milliarcsecond optical imaging*, [arXiv:1207.0808](https://arxiv.org/abs/1207.0808).
    3. Nuñez et al. (2012), *Imaging sub-milliarcsecond stellar features...*, [arXiv:1205.5743](https://arxiv.org/abs/1205.5743).
    4. Rou et al. (2013), *Monte-Carlo simulation of stellar intensity interferometry*, [arXiv:1303.4023](https://arxiv.org/abs/1303.4023).
    5. Zhang et al. (2025), LACT 32 镜/6 m 与 grouped layout 间距研究, [arXiv:2409.14382](https://arxiv.org/abs/2409.14382).
    6. H.E.S.S.、MAGIC、VERITAS 连续读出实现：[arXiv:2312.08015](https://arxiv.org/abs/2312.08015)、[arXiv:2402.04755](https://arxiv.org/abs/2402.04755)、[arXiv:1908.03095](https://arxiv.org/abs/1908.03095).
    """),
    code(r"""
    from __future__ import annotations

    import math
    import os
    import platform
    import sys
    from dataclasses import dataclass, asdict
    from pathlib import Path
    from itertools import combinations

    import matplotlib
    os.environ.setdefault("MPLCONFIGDIR", str(Path.cwd() / ".mplconfig"))
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import scipy
    from scipy.optimize import brentq
    from scipy.ndimage import gaussian_filter
    from scipy.special import j1
    from IPython.display import display, Markdown

    REPO_ROOT = Path.cwd().resolve()
    if not (REPO_ROOT / "python" / "sii_reconstruction.py").exists():
        REPO_ROOT = Path.cwd().resolve().parent
    assert (REPO_ROOT / "python" / "sii_reconstruction.py").exists(), REPO_ROOT
    sys.path.insert(0, str(REPO_ROOT / "python"))
    OUTPUT_DIR = REPO_ROOT / "run_logs" / "sii_paper_notebook"
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    RNG = np.random.default_rng(20260823)
    MAS_TO_RAD = np.deg2rad(1.0 / 3_600_000.0)
    RAD_TO_MAS = 1.0 / MAS_TO_RAD
    C = 299_792_458.0
    H_PLANCK = 6.626_070_15e-34
    SIDEREAL_DAY_S = 86_164.0905
    JY = 1e-26

    BLUE = "#2463A8"
    ORANGE = "#D97A2B"
    GOLD = "#B9972E"
    INK = "#25313C"
    GREY = "#98A2AD"
    plt.rcParams.update({
        "figure.dpi": 120,
        "savefig.dpi": 180,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.titleweight": "semibold",
        "axes.labelcolor": INK,
        "text.color": INK,
        "font.size": 10,
    })

    print({
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "matplotlib": matplotlib.__version__,
        "repo": str(REPO_ROOT),
        "output": str(OUTPUT_DIR),
    })
    """),
    md(r"""
    ### 1.1 可见的参数和假设

    默认采用 400 nm、2 nm 光学带宽、200 MHz 电子相关带宽、主程序给出的 24.576860 m² 遮挡后有效面积、20% 总通光/探测效率和 70.527 MHz/pixel 暗天 NSB 预测。实测 SPE 电荷样本给出 RMS 因子 1.016142。站点纬度使用仓库验证配置的 (29.36^\circ)。10 ns 微单元恢复时间、ADC 和电子噪声是待标定工程假设，下面单独标识。
    """),
    code(r"""
    @dataclass(frozen=True)
    class Instrument:
        wavelength_nm: float = 400.0
        optical_width_nm: float = 2.0
        telescope_diameter_m: float = 6.0
        effective_area_m2: float = 24.576860
        throughput: float = 0.20
        electronics_bandwidth_hz: float = 200e6
        adc_sample_rate_hz: float = 625e6
        detected_nsb_rate_hz: float = 70.527e6
        excess_noise_factor: float = 1.016142
        microcells_per_pixel: int = 270336
        microcell_recovery_time_ns: float = 10.0
        electronic_noise_rms_mv: float = 0.35
        adc_bits: int = 8
        adc_full_scale_mv: float = 200.0
        polarization_factor: float = 0.5
        spectral_shape_factor: float = 0.842
        calibration_floor_visibility2: float = 0.015
        site_lat_deg: float = 29.36
        source_dec_deg: float = 22.0
        observing_hours: float = 6.0
        segment_s: float = 1200.0
        source_ab_magnitude: float = 2.0

        @property
        def wavelength_m(self): return self.wavelength_nm * 1e-9
        @property
        def optical_width_m(self): return self.optical_width_nm * 1e-9
        @property
        def collecting_area_m2(self): return self.effective_area_m2
        @property
        def optical_bandwidth_hz(self):
            return C * self.optical_width_m / self.wavelength_m**2
        @property
        def sample_width_s(self): return 1.0 / self.adc_sample_rate_hz
        @property
        def coherence_area_s(self):
            return (self.polarization_factor * self.spectral_shape_factor
                    * self.wavelength_m**2 / (C * self.optical_width_m))

    @dataclass(frozen=True)
    class BinarySource:
        separation_mas: float = 0.20
        position_angle_deg: float = 35.0  # east of north
        flux_ratio_secondary_to_primary: float = 0.55
        primary_diameter_mas: float = 0.060
        secondary_diameter_mas: float = 0.040

    INST = Instrument()
    SOURCE = BinarySource()
    display(pd.DataFrame([asdict(INST)]).T.rename(columns={0: "value"}))
    display(pd.DataFrame([asdict(SOURCE)]).T.rename(columns={0: "value"}))
    """),
    code(r"""
    parameter_status = pd.DataFrame([
        ["32 telescope coordinates", "available", "layout_0803 run-card export", "used directly"],
        ["mirror area after obstruction", "available", "main NSB model: 24.576860 m²", "used directly"],
        ["mirror/filter/PDE curves", "available", "main configs/efficiency", "collapsed to 20% reference throughput"],
        ["measured SPE template", "available", "537 clean pulses; 84.0349557 mV ns/PE", "used in short waveform"],
        ["SPE charge dispersion", "available", "RMS factor 1.016142", "used in long-exposure ENF"],
        ["dark-sky NSB", "model estimate", "70.527 MHz/pixel under documented SkyCalc conditions", "scenario input"],
        ["S17351 microcell geometry", "available", "8 channels, 270336 cells/pixel, 25 µm pitch", "used directly"],
        ["S17351 recovery time", "not publicly specified", "10 ns provisional; Hamamatsu generic scale ~15 ns", "scan 1/10/30 ns"],
        ["ADC sample rate/full scale/bits", "engineering assumption", "625 MS/s, 200 mV, 8 bit from old v2", "used in short waveform"],
        ["electronic noise", "engineering assumption", "0.35 mV RMS from old v2", "used in short waveform"],
        ["electronics correlation bandwidth", "engineering assumption", "200 MHz", "sensitivity scenario"],
        ["crosstalk/afterpulse/dark rate", "missing", "requires temperature/overvoltage calibration", "not enabled"],
        ["clock drift/residual delay", "missing", "requires timing distribution and calibration runs", "0.2 ns illustrative residual"],
        ["narrow-band filter angular response", "missing", "2 nm top-hat approximation", "systematic caveat"],
    ], columns=["parameter", "status", "source_or_value", "treatment"])
    display(parameter_status)
    parameter_status.to_csv(OUTPUT_DIR / "parameter_availability.csv", index=False)
    """),
    md(r"""
    ## 2. Array data: 仓库阵列与 layout_0803 实际 32 镜坐标

    32 镜数据来自本工作区历史运行日志 `main_userv2_comparison_20260820/main_corsika.log:248–279`，其完整记录了 `layout_0803` 的 `TELESCOPE ... TEL.1`至`TEL.32`卡。同一坐标块在 2026-08-15至21 的 8 份独立本地运行日志中逐行一致（规范化后 SHA-256 前缀 `bc7f1cfa970b2d58`）。日志第 59–60 行声明输入系为 `corsika_nwu_relative`。CORSIKA 卡坐标单位为 cm，转为本 notebook 的 ENU 系时使用

    \[
    E=-y_{\rm west}/100,\qquad N=x_{\rm north}/100,\qquad U=z_{\rm up}/100.
    \]

    导出的 CSV 同时保留原始 NWU 列、ENU 列、CORSIKA 的 1-based 镜号和 LACTsim 的 0-based 索引，以便审计。
    """),
    code(r"""
    def load_branch_layout(path: Path) -> pd.DataFrame:
        frame = pd.read_csv(path)
        return frame.rename(columns={
            "position_x_m": "east_m",
            "position_y_m": "north_m",
            "position_z_m": "up_m",
        })[["telescope_id", "name", "east_m", "north_m", "up_m"]]

    def load_layout_0803_reco32(path: Path) -> pd.DataFrame:
        frame = pd.read_csv(path)
        assert len(frame) == 32
        assert frame.telescope_id.tolist() == list(range(1, 33))
        assert frame.lactsim_index.tolist() == list(range(32))
        assert frame.telescope_id.is_unique and frame.lactsim_index.is_unique
        assert np.allclose(frame.east_m, -frame.corsika_west_cm / 100.0)
        assert np.allclose(frame.north_m, frame.corsika_north_cm / 100.0)
        assert np.allclose(frame.up_m, frame.corsika_up_cm / 100.0)
        return frame

    def baseline_table(layout: pd.DataFrame) -> pd.DataFrame:
        rows = []
        for i, j in combinations(range(len(layout)), 2):
            a, b = layout.iloc[i], layout.iloc[j]
            vec = b[["east_m", "north_m", "up_m"]].to_numpy(float) - a[["east_m", "north_m", "up_m"]].to_numpy(float)
            rows.append({
                "i": int(a.telescope_id), "j": int(b.telescope_id),
                "name_i": a["name"], "name_j": b["name"],
                "east_m": vec[0], "north_m": vec[1], "up_m": vec[2],
                "baseline_m": np.linalg.norm(vec),
            })
        return pd.DataFrame(rows)

    BRANCH_LAYOUT = load_branch_layout(REPO_ROOT / "configs" / "arrays" / "star_interferometry_v2_array.csv")
    LACT32 = load_layout_0803_reco32(
        REPO_ROOT / "configs" / "arrays" / "layout_0803_reco32_coordinates.csv"
    )
    branch_baselines = baseline_table(BRANCH_LAYOUT)
    lact_baselines = baseline_table(LACT32)

    fig, axes = plt.subplots(1, 3, figsize=(13.2, 3.7))
    for axis, layout, title in [
        (axes[0], BRANCH_LAYOUT, "Branch validation layout (7 telescopes)"),
        (axes[1], LACT32, "layout_0803 production coordinates (32 telescopes)"),
    ]:
        axis.scatter(layout.east_m, layout.north_m, s=34, color=BLUE, edgecolor="white", linewidth=0.6)
        axis.set(title=title, xlabel="East [m]", ylabel="North [m]")
        axis.set_aspect("equal", adjustable="box")
        axis.grid(alpha=0.18)
    axes[2].hist(branch_baselines.baseline_m, bins=15, alpha=0.8, color=GOLD, label="7-telescope branch")
    axes[2].hist(lact_baselines.baseline_m, bins=25, alpha=0.70, color=BLUE, label="layout_0803 reco32")
    axes[2].set(title="Physical baseline distribution", xlabel="Baseline length [m]", ylabel="Pair count")
    axes[2].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "array_layouts_and_baselines.png", bbox_inches="tight")
    plt.show()

    layout_summary = pd.DataFrame([
        {"layout": "branch_7", "telescopes": len(BRANCH_LAYOUT), "pairs": len(branch_baselines),
         "min_nonzero_baseline_m": branch_baselines.loc[branch_baselines.baseline_m > 0, "baseline_m"].min(),
         "max_baseline_m": branch_baselines.baseline_m.max()},
        {"layout": "layout_0803_reco32", "telescopes": len(LACT32), "pairs": len(lact_baselines),
         "min_nonzero_baseline_m": lact_baselines.loc[lact_baselines.baseline_m > 0, "baseline_m"].min(),
         "max_baseline_m": lact_baselines.baseline_m.max()},
    ])
    display(layout_summary.round(3))
    layout_summary.to_csv(OUTPUT_DIR / "array_layout_summary.csv", index=False)
    """),
    md(r"""
    ## 3. Data: 理论源天图与完整理论 (uv) 图

    双星的两颗分量均为均匀圆盘。单个角直径为 \(\theta\) 的圆盘复可见度为

    \[
    V_{\rm UD}(q)=\frac{2J_1(\pi\theta q)}{\pi\theta q},\quad q=\sqrt{u^2+v^2}.
    \]

    两分量带上各自位置的傅里叶相位后按流量相加。天图和公式由同一参数对象生成，避免“画的源”和“模拟的可见度”不一致。
    """),
    code(r"""
    def uniform_disk_visibility(q_lambda, diameter_mas):
        x = np.pi * diameter_mas * MAS_TO_RAD * np.asarray(q_lambda)
        out = np.ones_like(x, dtype=float)
        mask = np.abs(x) > 1e-12
        out[mask] = 2.0 * j1(x[mask]) / x[mask]
        return out

    def binary_offsets_rad(source=SOURCE):
        pa = np.deg2rad(source.position_angle_deg)
        dx = 0.5 * source.separation_mas * np.sin(pa) * MAS_TO_RAD
        dy = 0.5 * source.separation_mas * np.cos(pa) * MAS_TO_RAD
        return (-dx, -dy), (dx, dy)

    def binary_visibility(u_lambda, v_lambda, source=SOURCE):
        u = np.asarray(u_lambda, dtype=float)
        v = np.asarray(v_lambda, dtype=float)
        q = np.hypot(u, v)
        (x1, y1), (x2, y2) = binary_offsets_rad(source)
        primary = uniform_disk_visibility(q, source.primary_diameter_mas) * np.exp(-2j * np.pi * (u*x1 + v*y1))
        secondary = uniform_disk_visibility(q, source.secondary_diameter_mas) * np.exp(-2j * np.pi * (u*x2 + v*y2))
        ratio = source.flux_ratio_secondary_to_primary
        return (primary + ratio * secondary) / (1.0 + ratio)

    def render_binary_sky(source=SOURCE, fov_mas=0.70, pixels=192):
        theta = np.linspace(-fov_mas/2, fov_mas/2, pixels)
        xx, yy = np.meshgrid(theta, theta)
        (x1, y1), (x2, y2) = binary_offsets_rad(source)
        x1, y1, x2, y2 = np.array([x1, y1, x2, y2]) * RAD_TO_MAS
        image = ((xx-x1)**2 + (yy-y1)**2 <= (source.primary_diameter_mas/2)**2).astype(float)
        component2 = ((xx-x2)**2 + (yy-y2)**2 <= (source.secondary_diameter_mas/2)**2).astype(float)
        if component2.sum() > 0:
            component2 *= source.flux_ratio_secondary_to_primary * image.sum() / component2.sum()
        image += component2
        image /= image.sum()
        return theta, image

    theta_sky, sky_image = render_binary_sky()
    extent = [theta_sky[0], theta_sky[-1], theta_sky[0], theta_sky[-1]]
    fig, axis = plt.subplots(figsize=(5.0, 4.2))
    artist = axis.imshow(sky_image, origin="lower", extent=extent, cmap="magma")
    axis.set(title="Theoretical source sky image", xlabel=r"$\Delta\alpha\cos\delta$ [mas]", ylabel=r"$\Delta\delta$ [mas]")
    axis.set_aspect("equal")
    fig.colorbar(artist, ax=axis, label="Normalized surface brightness / pixel")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "theoretical_source_sky.png", bbox_inches="tight")
    plt.show()
    """),
    code(r"""
    bmax_m = lact_baselines.baseline_m.max()
    uv_limit = 1.05 * bmax_m / INST.wavelength_m
    uv_axis = np.linspace(-uv_limit, uv_limit, 301)
    uu, vv = np.meshgrid(uv_axis, uv_axis)
    theory_power = np.abs(binary_visibility(uu, vv))**2

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.4))
    im0 = axes[0].imshow(theory_power, origin="lower",
                         extent=np.array([-uv_limit, uv_limit, -uv_limit, uv_limit])/1e6,
                         cmap="viridis", vmin=0, vmax=1)
    axes[0].set(title=r"Theoretical $|V(u,v)|^2$", xlabel=r"$u$ [M$\lambda$]", ylabel=r"$v$ [M$\lambda$]")
    axes[0].set_aspect("equal")
    fig.colorbar(im0, ax=axes[0], label=r"$|V|^2$")
    im1 = axes[1].imshow(np.log10(np.clip(theory_power, 1e-4, 1)), origin="lower",
                         extent=np.array([-uv_limit, uv_limit, -uv_limit, uv_limit])/1e6,
                         cmap="cividis", vmin=-4, vmax=0)
    axes[1].set(title="Same power plane (log scale)", xlabel=r"$u$ [M$\lambda$]", ylabel=r"$v$ [M$\lambda$]")
    axes[1].set_aspect("equal")
    fig.colorbar(im1, ax=axes[1], label=r"$\log_{10}|V|^2$")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "theoretical_uv_power.png", bbox_inches="tight")
    plt.show()

    assert np.isclose(abs(binary_visibility(0.0, 0.0))**2, 1.0, atol=1e-12)
    assert np.allclose(theory_power, theory_power[::-1, ::-1], atol=2e-12)
    """),
    md(r"""
    ## 4. 两台望远镜之间到底怎么算

    对望远镜 (i,j)，先在本地 ENU 坐标中作基线

    \[
    \mathbf B_{ij}=\mathbf r_j-\mathbf r_i.
    \]

    给定站点纬度 \(\phi\)、源赤纬 \(\delta\) 和时角 \(H\)，源方向在 ENU 中为

    \[
    \hat s=(\cos\delta\sin H,\;\sin\delta\cos\phi-\cos\delta\cos H\sin\phi,\;
    \sin\delta\sin\phi+\cos\delta\cos H\cos\phi).
    \]

    再定义天球东向 \(\hat u\propto \hat z\times\hat s\)、天球北向 \(\hat v=\hat s\times\hat u\)，得到

    \[
    (u,v,w)=(\mathbf B\cdot\hat u,\mathbf B\cdot\hat v,\mathbf B\cdot\hat s),
    \quad \tau_g=w/c.
    \]

    这一步必须包含站点纬度；直接把本地 ENU 代入不含纬度的地固坐标矩阵会得到错误的轨迹。
    """),
    code(r"""
    def source_direction_enu(hour_angle_rad, dec_rad, lat_rad):
        return np.array([
            np.cos(dec_rad) * np.sin(hour_angle_rad),
            np.sin(dec_rad) * np.cos(lat_rad) - np.cos(dec_rad) * np.cos(hour_angle_rad) * np.sin(lat_rad),
            np.sin(dec_rad) * np.sin(lat_rad) + np.cos(dec_rad) * np.cos(hour_angle_rad) * np.cos(lat_rad),
        ])

    def uvw_from_enu(baseline_enu_m, hour_angle_rad, dec_rad, lat_rad):
        source = source_direction_enu(hour_angle_rad, dec_rad, lat_rad)
        source = source / np.linalg.norm(source)
        up = np.array([0.0, 0.0, 1.0])
        u_axis = np.cross(up, source)
        if np.linalg.norm(u_axis) < 1e-12:
            u_axis = np.array([1.0, 0.0, 0.0])
        else:
            u_axis /= np.linalg.norm(u_axis)
        v_axis = np.cross(source, u_axis)
        baseline = np.asarray(baseline_enu_m, float)
        return np.array([baseline @ u_axis, baseline @ v_axis, baseline @ source])

    def legacy_v2_uvw_labeled_enu(baseline_enu_m, hour_angle_rad, dec_rad):
        east, north, up = baseline_enu_m
        sh, ch = np.sin(hour_angle_rad), np.cos(hour_angle_rad)
        sd, cd = np.sin(dec_rad), np.cos(dec_rad)
        return np.array([
            sh*east + ch*north,
            -sd*ch*east + sd*sh*north + cd*up,
            cd*ch*east - cd*sh*north + sd*up,
        ])

    lat = np.deg2rad(INST.site_lat_deg)
    dec = np.deg2rad(INST.source_dec_deg)
    hour_angles_h = np.linspace(-INST.observing_hours/2, INST.observing_hours/2,
                                int(INST.observing_hours*3600/INST.segment_s)+1)
    hour_angles_rad = hour_angles_h * np.pi / 12.0

    farthest = lact_baselines.loc[lact_baselines.baseline_m.idxmax()]
    baseline_vec = farthest[["east_m", "north_m", "up_m"]].to_numpy(float)
    pair_track = np.array([uvw_from_enu(baseline_vec, h, dec, lat) for h in hour_angles_rad])
    legacy_track = np.array([legacy_v2_uvw_labeled_enu(baseline_vec, h, dec) for h in hour_angles_rad])
    projected = np.hypot(pair_track[:, 0], pair_track[:, 1])
    legacy_projected = np.hypot(legacy_track[:, 0], legacy_track[:, 1])
    pair_v2 = np.abs(binary_visibility(pair_track[:, 0]/INST.wavelength_m,
                                       pair_track[:, 1]/INST.wavelength_m))**2

    transit = int(np.argmin(np.abs(hour_angles_h)))
    pair_example = pd.DataFrame([{
        "telescope_i": farthest.name_i,
        "telescope_j": farthest.name_j,
        "B_E_m": baseline_vec[0], "B_N_m": baseline_vec[1], "B_U_m": baseline_vec[2],
        "physical_B_m": np.linalg.norm(baseline_vec),
        "transit_u_m": pair_track[transit, 0], "transit_v_m": pair_track[transit, 1],
        "transit_w_m": pair_track[transit, 2],
        "transit_projected_B_m": projected[transit],
        "transit_geometric_delay_ns": pair_track[transit, 2]/C*1e9,
        "transit_visibility2": pair_v2[transit],
    }])
    display(pair_example.T.rename(columns={0: "value"}).round(6))
    pair_example.to_csv(OUTPUT_DIR / "farthest_pair_calculation.csv", index=False)

    fig, axes = plt.subplots(1, 3, figsize=(13.0, 3.6))
    axes[0].plot(pair_track[:,0]/1e3, pair_track[:,1]/1e3, "o-", color=BLUE, label="correct ENU")
    axes[0].plot(legacy_track[:,0]/1e3, legacy_track[:,1]/1e3, "s--", color=ORANGE, label="v2 formula on ENU")
    axes[0].set(title=f"Pair {farthest.name_i}–{farthest.name_j}: uv track", xlabel=r"$u$ [km]", ylabel=r"$v$ [km]")
    axes[0].set_aspect("equal", adjustable="datalim")
    axes[0].legend(frameon=False, fontsize=8)
    axes[1].plot(hour_angles_h, projected, "o-", color=BLUE, label="correct")
    axes[1].plot(hour_angles_h, legacy_projected, "s--", color=ORANGE, label="v2-on-ENU")
    axes[1].set(title="Projected baseline", xlabel="Hour angle [h]", ylabel=r"$B_\perp$ [m]")
    axes[1].legend(frameon=False, fontsize=8)
    axes[2].plot(hour_angles_h, pair_track[:,2]/C*1e6, "o-", color=GOLD, label=r"$w/c$")
    ax2 = axes[2].twinx()
    ax2.plot(hour_angles_h, pair_v2, "s--", color=BLUE, label=r"$|V|^2$")
    axes[2].set(title="Delay and source visibility", xlabel="Hour angle [h]", ylabel="Geometric delay [µs]")
    ax2.set_ylabel(r"True $|V|^2$", color=BLUE)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "two_telescope_calculation.png", bbox_inches="tight")
    plt.show()

    invariant_error = np.max(np.abs(np.sum(pair_track**2, axis=1) - np.dot(baseline_vec, baseline_vec)))
    assert invariant_error < 1e-6
    """),
    md(r"""
    ## 5. Results: 全阵列理论 (uv) 覆盖

    下图把所有 (N(N-1)/2=496) 条物理基线在 6 小时内的地球自转轨迹叠加到理论功率面。图中同时画出共轭点 ((-u,-v))；这不是额外的独立测量，而是实天图傅里叶功率的厄米对称性。
    """),
    code(r"""
    def make_uv_coverage(layout, hour_angles_rad, instrument=INST):
        rows = []
        baselines = baseline_table(layout)
        for _, pair in baselines.iterrows():
            vec = pair[["east_m", "north_m", "up_m"]].to_numpy(float)
            for segment, hour_angle in enumerate(hour_angles_rad):
                u_m, v_m, w_m = uvw_from_enu(vec, hour_angle, dec, lat)
                rows.append({
                    "telescope_i": pair.name_i, "telescope_j": pair.name_j,
                    "segment": segment, "hour_angle_h": hour_angle*12/np.pi,
                    "u_m": u_m, "v_m": v_m, "w_m": w_m,
                    "u_lambda": u_m/instrument.wavelength_m,
                    "v_lambda": v_m/instrument.wavelength_m,
                    "projected_baseline_m": np.hypot(u_m, v_m),
                    "geometric_delay_ns": w_m/C*1e9,
                })
        return pd.DataFrame(rows)

    uv_coverage = make_uv_coverage(LACT32, hour_angles_rad)
    uv_coverage["visibility2_true"] = np.abs(binary_visibility(
        uv_coverage.u_lambda.to_numpy(), uv_coverage.v_lambda.to_numpy()))**2

    fig, axis = plt.subplots(figsize=(7.0, 6.0))
    bg = axis.imshow(theory_power, origin="lower",
                     extent=np.array([-uv_limit, uv_limit, -uv_limit, uv_limit])/1e6,
                     cmap="Greys", vmin=0, vmax=1, alpha=0.38)
    sc = axis.scatter(uv_coverage.u_lambda/1e6, uv_coverage.v_lambda/1e6,
                      c=uv_coverage.hour_angle_h, cmap="coolwarm", s=5, alpha=0.65, linewidth=0)
    axis.scatter(-uv_coverage.u_lambda/1e6, -uv_coverage.v_lambda/1e6,
                 c=uv_coverage.hour_angle_h, cmap="coolwarm", s=3, alpha=0.25, linewidth=0)
    axis.set(title="layout_0803 reco32 earth-rotation uv coverage",
             xlabel=r"$u$ [M$\lambda$]", ylabel=r"$v$ [M$\lambda$]")
    axis.set_aspect("equal", adjustable="box")
    fig.colorbar(sc, ax=axis, label="Hour angle [h]")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "theoretical_uv_with_array_coverage.png", bbox_inches="tight")
    plt.show()
    display(uv_coverage.head())
    print(f"uv rows = {len(uv_coverage):,}; unique physical pairs = {uv_coverage[['telescope_i','telescope_j']].drop_duplicates().shape[0]}")
    pd.DataFrame([{
        "uv_rows": len(uv_coverage),
        "unique_physical_pairs": uv_coverage[["telescope_i", "telescope_j"]].drop_duplicates().shape[0],
        "maximum_projected_baseline_m": uv_coverage.projected_baseline_m.max(),
        "maximum_abs_geometric_delay_ns": uv_coverage.geometric_delay_ns.abs().max(),
    }]).to_csv(OUTPUT_DIR / "uv_coverage_summary.csv", index=False)
    """),
    md(r"""
    ## 6. 从光学相干到电子相关峰

    对非偏振、近似带宽为 (Delta\lambda) 的通带，未分辨相关峰的面积尺度为

    \[
    \tau_{c,\mathrm{area}}=p\,s_{\rm band}\frac{\lambda^2}{c\,\Delta\lambda}.
    \]

    它是飞秒量级，远小于 625 MS/s 的 1.6 ns 采样。电子学看到的零时延峰高度因此被稀释约 (	au_c/\Delta t)，而不是把相干时间人为设成 1 ns。下图用归一化电子相关核画出这条基线在过中天时的理论峰，并展示几何时延校正前后的峰位。
    """),
    code(r"""
    coherence_fs = INST.coherence_area_s * 1e15
    unresolved_contrast = INST.coherence_area_s / INST.sample_width_s
    star_fraction_at_m2 = None

    lag_ns = np.linspace(-6000, 6000, 6001)
    electronics_sigma_ns = 2.0
    kernel = np.exp(-0.5*(lag_ns/electronics_sigma_ns)**2)
    kernel /= np.trapezoid(kernel, lag_ns*1e-9)  # unit area in s^-1
    transit_delay_ns = pair_example.loc[0, "transit_geometric_delay_ns"]
    transit_visibility2 = pair_example.loc[0, "transit_visibility2"]
    shifted_kernel = np.interp(lag_ns-transit_delay_ns, lag_ns, kernel, left=0, right=0)
    peak_corrected = INST.coherence_area_s * transit_visibility2 * kernel
    peak_uncorrected = INST.coherence_area_s * transit_visibility2 * shifted_kernel

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 3.7))
    axes[0].plot(lag_ns, peak_uncorrected, color=ORANGE, label="before geometric-delay correction")
    axes[0].plot(lag_ns, peak_corrected, color=BLUE, label="after correction")
    axes[0].set(title=r"Expected electronic $g^{(2)}(\tau)-1$", xlabel="Lag [ns]", ylabel=r"$g^{(2)}-1$")
    axes[0].legend(frameon=False, fontsize=8)
    zoom = np.abs(lag_ns) < 20
    axes[1].plot(lag_ns[zoom], peak_corrected[zoom], color=BLUE)
    axes[1].axhline(unresolved_contrast*transit_visibility2, color=GREY, linestyle="--",
                    label=r"simple $\tau_c/\Delta t$ scale")
    axes[1].set(title="Corrected peak (zoom)", xlabel="Lag [ns]", ylabel=r"$g^{(2)}-1$")
    axes[1].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "two_telescope_g2_peak.png", bbox_inches="tight")
    plt.show()

    coherence_summary = pd.DataFrame([{
        "coherence_area_fs": coherence_fs,
        "adc_sample_width_ns": INST.sample_width_s*1e9,
        "unresolved_zero_baseline_contrast": unresolved_contrast,
        "pair_transit_visibility2": transit_visibility2,
        "pair_expected_sample_contrast": unresolved_contrast*transit_visibility2,
    }])
    display(coherence_summary.T.rename(columns={0:"value"}))
    """),
    md(r"""
    ### 6.1 170 ns 完整电子学记录与微单元恢复

    这一节回答“完整波形有没有意义”：有，但用途是校准延迟、SPE 形状、ADC、基线噪声、饱和和恢复偏差；它不能单独承担小时级天文灵敏度。170 ns 在 625 MS/s 下只有约 106 个采样，而且单采样预期 HBT 对比度仅约 (10^{-6})，所以单帧相关曲线必然由随机噪声主导。

    输入不是人为的 1 ns 光学波形，而是每个 ADC bin 内满足正确一、二阶矩的光电子计数：两镜恒星流包含一个极小的共享 Poisson 分量，使协方差等于 (r_{\star,1}r_{\star,2}\tau_c\Delta t|V|^2)；星光余项和 NSB 独立。每个 p.e. 随机落到 270336 个微单元，电荷按

    \[
    q(\Delta t)=1-\exp(-\Delta t/\tau_{\rm rec})
    \]

    恢复，然后叠加主程序的实测 SPE 模板、0.35 mV 假设电子噪声并做 8-bit/200 mV ADC 量化。几何时延 (w/c) 先由上节计算并用数字缓冲粗校正；这里仅保留 0.2 ns 的示意残差。
    """),
    code(r"""
    from sii_unified import (
        Instrument as ShortWaveformInstrument,
        detected_star_rate_hz as unified_detected_star_rate_hz,
        load_measured_spe_template,
        mean_recovery_fraction,
        normalized_cross_correlation,
        simulate_short_pair_waveforms,
    )

    short_instrument = ShortWaveformInstrument(
        wavelength_nm=INST.wavelength_nm,
        optical_width_nm=INST.optical_width_nm,
        effective_area_m2=INST.effective_area_m2,
        throughput=INST.throughput,
        electronics_bandwidth_hz=INST.electronics_bandwidth_hz,
        adc_sample_rate_hz=INST.adc_sample_rate_hz,
        adc_bits=INST.adc_bits,
        adc_full_scale_mv=INST.adc_full_scale_mv,
        detected_nsb_rate_hz=INST.detected_nsb_rate_hz,
        electronic_noise_rms_mv=INST.electronic_noise_rms_mv,
        excess_noise_factor=INST.excess_noise_factor,
        polarization_factor=INST.polarization_factor,
        spectral_shape_factor=INST.spectral_shape_factor,
        microcells_per_pixel=INST.microcells_per_pixel,
        microcell_recovery_time_ns=INST.microcell_recovery_time_ns,
    )
    spe_t_ns, spe_mv = load_measured_spe_template(
        REPO_ROOT / "configs" / "electronics" / "parameters" / "spe_template_measured.csv"
    )
    reference_star_rate = unified_detected_star_rate_hz(
        INST.source_ab_magnitude, short_instrument)
    short_record = simulate_short_pair_waveforms(
        duration_ns=170.0,
        star_rate_hz=reference_star_rate,
        nsb_rate_hz=INST.detected_nsb_rate_hz,
        visibility2=float(transit_visibility2),
        instrument=short_instrument,
        template_time_ns=spe_t_ns,
        template_amplitude_mv=spe_mv,
        delay_ns=0.2,
        recovery_time_ns=10.0,
        seed=20260824,
    )

    recovery_rows = []
    total_rate = reference_star_rate + INST.detected_nsb_rate_hz
    for tau_ns in [1.0, 10.0, 30.0]:
        record = simulate_short_pair_waveforms(
            170.0, reference_star_rate, INST.detected_nsb_rate_hz,
            float(transit_visibility2), short_instrument, spe_t_ns, spe_mv,
            delay_ns=0.2, recovery_time_ns=tau_ns, seed=20260824)
        recovery_rows.append({
            "recovery_time_ns": tau_ns,
            "analytic_mean_charge_fraction": mean_recovery_fraction(
                total_rate, INST.microcells_per_pixel, tau_ns),
            "record_mean_fraction_tel_A": np.mean(record["recovery_a"]),
            "record_mean_fraction_tel_B": np.mean(record["recovery_b"]),
            "record_pe_tel_A": len(record["pe_times_a_ns"]),
            "record_pe_tel_B": len(record["pe_times_b_ns"]),
        })
    recovery_table = pd.DataFrame(recovery_rows)
    display(recovery_table)
    recovery_table.to_csv(OUTPUT_DIR / "microcell_recovery_sensitivity.csv", index=False)

    lags, short_corr = normalized_cross_correlation(
        short_record["adc_a_mv"], short_record["adc_b_mv"])
    lag_ns_short = lags * INST.sample_width_s * 1e9
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 3.7))
    axes[0].plot(spe_t_ns, spe_mv, color=GOLD)
    axes[0].set(title="Measured common SPE template", xlabel="Time from peak [ns]", ylabel="Amplitude [mV]")
    axes[1].step(short_record["sample_time_ns"], short_record["adc_a_mv"], where="mid", color=BLUE, label="telescope A")
    axes[1].step(short_record["sample_time_ns"], short_record["adc_b_mv"], where="mid", color=ORANGE, alpha=.8, label="telescope B")
    axes[1].set(title="One 170 ns digitized record", xlabel="Time [ns]", ylabel="ADC-equivalent [mV]")
    axes[1].legend(frameon=False, fontsize=8)
    axes[2].plot(lag_ns_short, short_corr, color=BLUE)
    axes[2].axvline(0, color=GREY, linestyle="--", linewidth=1)
    axes[2].set(title="Single-record correlation: noise dominated", xlabel="Lag [ns]", ylabel="Normalized correlation")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "short_waveform_and_recovery.png", bbox_inches="tight")
    plt.show()

    short_waveform_summary = pd.DataFrame([{
        "duration_ns": 170.0,
        "samples": len(short_record["sample_time_ns"]),
        "star_rate_MHz": reference_star_rate / 1e6,
        "nsb_rate_MHz": INST.detected_nsb_rate_hz / 1e6,
        "expected_sample_contrast": short_record["expected_sample_contrast"],
        "shared_correlated_counts_expected_in_record": (
            short_record["shared_count_mean_per_bin"] * len(short_record["sample_time_ns"])),
        "adc_clipped_fraction_A": np.mean(np.abs(short_record["adc_a_mv"]) >= INST.adc_full_scale_mv/2 - INST.adc_full_scale_mv/2**INST.adc_bits),
    }])
    display(short_waveform_summary.T.rename(columns={0:"value"}))
    short_waveform_summary.to_csv(OUTPUT_DIR / "short_waveform_summary.csv", index=False)
    """),
    md(r"""
    ## 7. 模拟的 (uv) 测量：星等、NSB、电子带宽与误差

    AB 星等先转成每单位频率光子谱密度 (n_\nu=F_\nu/(h\nu))，再得到已探测星光率 (r_\star=A\eta n_\nu\Delta\nu)。每个 20 min 段生成一个不裁剪的高斯近似 (|V|^2) 估计和显式 `sigma_visibility2`。这是长曝光相关器的充分统计量，不生成一整夜的 625 MS/s 原始波形。
    """),
    code(r"""
    def ab_photon_spectral_density(magnitude, wavelength_m=INST.wavelength_m):
        f_nu = 3631.0 * JY * 10**(-0.4*magnitude)
        return f_nu / (H_PLANCK * C / wavelength_m)

    def detected_star_rate_hz(magnitude, instrument=INST):
        return (instrument.collecting_area_m2 * instrument.throughput
                * ab_photon_spectral_density(magnitude, instrument.wavelength_m)
                * instrument.optical_bandwidth_hz)

    def unit_visibility_snr(magnitude, integration_s, instrument=INST,
                            spectral_channels=1, nsb_rate_hz=None):
        nsb = instrument.detected_nsb_rate_hz if nsb_rate_hz is None else nsb_rate_hz
        star = detected_star_rate_hz(magnitude, instrument)
        total = star + nsb
        one_channel = (star**2 / (total * instrument.optical_bandwidth_hz)
                       * np.sqrt(instrument.electronics_bandwidth_hz * integration_s / 2)
                       / instrument.excess_noise_factor)
        return one_channel * np.sqrt(spectral_channels)

    unit_snr = unit_visibility_snr(INST.source_ab_magnitude, INST.segment_s)
    sigma_stat = 1.0 / unit_snr
    sigma_total = np.sqrt(sigma_stat**2 + INST.calibration_floor_visibility2**2)
    simulated_uv = uv_coverage.copy()
    simulated_uv["sigma_visibility2_stat"] = sigma_stat
    simulated_uv["sigma_visibility2"] = sigma_total
    simulated_uv["visibility2_measured"] = (
        simulated_uv.visibility2_true + RNG.normal(0, sigma_total, len(simulated_uv))
    )
    simulated_uv["expected_snr"] = simulated_uv.visibility2_true / sigma_total
    simulated_uv.to_csv(OUTPUT_DIR / "simulated_uv_measurements.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.4))
    sc = axes[0].scatter(simulated_uv.u_lambda/1e6, simulated_uv.v_lambda/1e6,
                         c=simulated_uv.visibility2_measured, cmap="viridis",
                         vmin=0, vmax=1, s=6, alpha=0.70, linewidth=0)
    axes[0].set(title=rf"Simulated uv measurements, $m_{{AB}}={INST.source_ab_magnitude:g}$",
                xlabel=r"$u$ [M$\lambda$]", ylabel=r"$v$ [M$\lambda$]")
    axes[0].set_aspect("equal", adjustable="box")
    fig.colorbar(sc, ax=axes[0], label=r"Measured $|V|^2$ (not clipped)")
    sample = simulated_uv.sample(min(1800, len(simulated_uv)), random_state=20260823)
    axes[1].scatter(sample.visibility2_true, sample.visibility2_measured,
                    s=8, alpha=0.35, color=BLUE)
    axes[1].plot([0,1], [0,1], color=INK, linewidth=1)
    axes[1].axhline(0, color=GREY, linewidth=0.8)
    axes[1].set(title="Measurement closure", xlabel=r"True $|V|^2$", ylabel=r"Measured $|V|^2$")
    axes[1].text(0.03, 0.94, f"unit-|V|² SNR/segment = {unit_snr:.1f}\nσ = {sigma_total:.3f}",
                 transform=axes[1].transAxes, va="top")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "simulated_uv_measurements.png", bbox_inches="tight")
    plt.show()

    measurement_summary = pd.DataFrame([{
        "source_ab_magnitude": INST.source_ab_magnitude,
        "detected_star_rate_MHz": detected_star_rate_hz(INST.source_ab_magnitude)/1e6,
        "detected_nsb_rate_MHz": INST.detected_nsb_rate_hz/1e6,
        "unit_visibility_snr_per_segment": unit_snr,
        "sigma_visibility2_total": sigma_total,
        "negative_measurement_fraction": np.mean(simulated_uv.visibility2_measured < 0),
    }])
    display(measurement_summary.T.rename(columns={0:"value"}))
    """),
    md(r"""
    ## 8. 无相位、非参数重建

    使用本分支 `python/sii_reconstruction.py` 的 `positive_softmax_multistart_lbfgsb`：

    - 像素强度由 softmax 参数化，严格非负且总流量为 1；
    - 只在有限圆形支撑域内优化；
    - 损失是带每点误差权重的 Huber (|V|^2) 残差，加二次平滑正则；
    - 多起点 L-BFGS-B 缓解非凸局部极值；
    - 真值绝不进入优化，仅在完成后对允许的平移和 180°镜像进行验证。

    为控制矩阵规模，先在 (uv) 平面作逆方差加权分箱；分箱是数据压缩，不插值、不制造未观测 (uv) 点。零基线点作为独立的总流量标定加入。
    """),
    code(r"""
    from sii_reconstruction import UvData, reconstruct_uv_data

    def bin_uv_measurements(frame, cell_mlambda=120.0):
        data = frame.copy()
        data["ku"] = np.round(data.u_lambda/1e6/cell_mlambda).astype(int)
        data["kv"] = np.round(data.v_lambda/1e6/cell_mlambda).astype(int)
        data["inverse_variance"] = 1.0 / data.sigma_visibility2**2
        data["weighted_value"] = data.visibility2_measured * data.inverse_variance
        grouped = data.groupby(["ku","kv"], as_index=False).agg(
            u_lambda=("u_lambda", "mean"),
            v_lambda=("v_lambda", "mean"),
            weighted_value=("weighted_value", "sum"),
            inverse_variance=("inverse_variance", "sum"),
            multiplicity=("visibility2_measured", "size"),
        )
        grouped["visibility2"] = grouped.weighted_value / grouped.inverse_variance
        grouped["sigma"] = np.sqrt(1.0/grouped.inverse_variance)
        return grouped

    binned = bin_uv_measurements(simulated_uv)
    binned = pd.concat([binned, pd.DataFrame([{
        "ku":0, "kv":0, "u_lambda":0.0, "v_lambda":0.0,
        "weighted_value":10_000.0, "inverse_variance":10_000.0,
        "multiplicity":1, "visibility2":1.0, "sigma":0.01,
    }])], ignore_index=True)
    weights = 1.0 / binned.sigma.to_numpy()**2
    weights = np.clip(weights/np.mean(weights), 0.1, 10.0)
    weights /= weights.mean()
    uv_data = UvData(
        u_lambda=binned.u_lambda.to_numpy(float),
        v_lambda=binned.v_lambda.to_numpy(float),
        visibility_abs2=binned.visibility2.to_numpy(float),
        sigma=binned.sigma.to_numpy(float),
        weight=weights,
        multiplicity=binned.multiplicity.to_numpy(int),
        input_rows=len(simulated_uv),
        finite_rows=len(simulated_uv),
        physical_violations=int(((binned.visibility2<0)|(binned.visibility2>1)).sum()),
    )

    result = reconstruct_uv_data(
        uv_data,
        grid_size=28,
        fov_mas=0.70,
        support_radius_mas=0.32,
        starts=4,
        max_iter=900,
        smoothness=0.020,
        huber_delta=0.15,
        seed=20260823,
        peak_minimum_separation_mas=0.10,
    )
    print({k:v for k,v in result.metrics.items() if not isinstance(v, np.ndarray)})
    print(f"binned uv points: {len(binned):,} from {len(simulated_uv):,} measurements")
    """),
    code(r"""
    def grid_truth_on_reconstruction(theta, source=SOURCE):
        xx, yy = np.meshgrid(theta, theta)
        (x1,y1),(x2,y2) = binary_offsets_rad(source)
        x1,y1,x2,y2 = np.array([x1,y1,x2,y2])*RAD_TO_MAS
        pixel = abs(theta[1]-theta[0])
        s1 = max(source.primary_diameter_mas/2, 0.65*pixel)
        s2 = max(source.secondary_diameter_mas/2, 0.65*pixel)
        truth = np.exp(-((xx-x1)**2+(yy-y1)**2)/(2*s1**2))
        comp2 = np.exp(-((xx-x2)**2+(yy-y2)**2)/(2*s2**2))
        comp2 *= source.flux_ratio_secondary_to_primary * truth.sum()/comp2.sum()
        truth += comp2
        return truth/truth.sum()

    def ambiguity_align(candidate, truth):
        best = None
        for mirrored, image in [(False,candidate), (True,candidate[::-1,::-1])]:
            corr = np.fft.ifft2(np.fft.fft2(truth)*np.conj(np.fft.fft2(image))).real
            iy, ix = np.unravel_index(np.argmax(corr), corr.shape)
            sy = iy if iy <= image.shape[0]//2 else iy-image.shape[0]
            sx = ix if ix <= image.shape[1]//2 else ix-image.shape[1]
            aligned = np.roll(image, (sy,sx), axis=(0,1))
            coefficient = np.corrcoef(aligned.ravel(), truth.ravel())[0,1]
            if best is None or coefficient > best[0]:
                best = coefficient, aligned, mirrored, (sy,sx)
        return best

    truth_reco = grid_truth_on_reconstruction(result.theta_mas)
    truth_corr_raw, aligned_reco_raw, used_mirror, shift = ambiguity_align(result.image, truth_reco)
    # Compare at the actual nominal fringe scale of this array, lambda/Bmax.
    reconstruction_pixel_mas = abs(result.theta_mas[1] - result.theta_mas[0])
    comparison_beam_fwhm_mas = INST.wavelength_m / lact_baselines.baseline_m.max() * RAD_TO_MAS
    comparison_beam_sigma_pixels = comparison_beam_fwhm_mas / (2.355 * reconstruction_pixel_mas)
    aligned_reco = gaussian_filter(aligned_reco_raw, comparison_beam_sigma_pixels)
    truth_reco = gaussian_filter(truth_reco, comparison_beam_sigma_pixels)
    aligned_reco /= aligned_reco.sum()
    truth_reco /= truth_reco.sum()
    truth_corr = np.corrcoef(aligned_reco.ravel(), truth_reco.ravel())[0,1]
    truth_nrmse = np.sqrt(np.mean((aligned_reco-truth_reco)**2))/(np.sqrt(np.mean(truth_reco**2))+1e-15)
    extent_reco = [result.theta_mas[0], result.theta_mas[-1], result.theta_mas[0], result.theta_mas[-1]]

    fig, axes = plt.subplots(2, 3, figsize=(13.8, 8.6), constrained_layout=True)
    sc = axes[0,0].scatter(uv_data.u_lambda/1e6, uv_data.v_lambda/1e6,
                           c=uv_data.visibility_abs2, cmap="viridis", vmin=0, vmax=1, s=12)
    axes[0,0].scatter(-uv_data.u_lambda/1e6, -uv_data.v_lambda/1e6,
                      c=uv_data.visibility_abs2, cmap="viridis", vmin=0, vmax=1, s=6, alpha=.3)
    axes[0,0].set(title="Binned measured uv power", xlabel=r"$u$ [M$\lambda$]", ylabel=r"$v$ [M$\lambda$]")
    axes[0,0].set_aspect("equal", adjustable="datalim")
    fig.colorbar(sc, ax=axes[0,0], label=r"measured $|V|^2$")
    lim = max(abs(result.dirty_autocorrelation.min()), abs(result.dirty_autocorrelation.max()))
    im = axes[0,1].imshow(result.dirty_autocorrelation, origin="lower", extent=extent_reco,
                          cmap="coolwarm", vmin=-lim, vmax=lim)
    axes[0,1].set(title="Dirty autocorrelation", xlabel="ΔRA [mas]", ylabel="ΔDec [mas]")
    fig.colorbar(im, ax=axes[0,1])
    axes[0,2].scatter(uv_data.visibility_abs2, result.predicted_visibility_abs2,
                      s=10, alpha=.45, color=BLUE)
    axes[0,2].plot([-0.2,1.1],[-0.2,1.1], color=INK, linewidth=1)
    axes[0,2].set(title=f"Forward closure, weighted RMSE={result.metrics['weighted_fit_rmse']:.3f}",
                  xlabel=r"Measured $|V|^2$", ylabel=r"Reconstructed $|V|^2$", xlim=(-.2,1.1), ylim=(-.2,1.1))
    vmax = max(truth_reco.max(), aligned_reco.max())
    for axis, image, title in [
        (axes[1,0], truth_reco, f"Truth at common {comparison_beam_fwhm_mas:.3f} mas beam"),
        (axes[1,1], aligned_reco, "Reconstruction, ambiguity-aligned"),
    ]:
        art = axis.imshow(image, origin="lower", extent=extent_reco, cmap="magma", vmin=0, vmax=vmax)
        axis.set(title=title, xlabel="ΔRA [mas]", ylabel="ΔDec [mas]")
        fig.colorbar(art, ax=axis)
    diff = aligned_reco-truth_reco
    dl = max(abs(diff.min()),abs(diff.max()))
    art = axes[1,2].imshow(diff, origin="lower", extent=extent_reco, cmap="coolwarm", vmin=-dl, vmax=dl)
    axes[1,2].set(title=f"Residual; corr={truth_corr:.3f}, NRMSE={truth_nrmse:.3f}", xlabel="ΔRA [mas]", ylabel="ΔDec [mas]")
    fig.colorbar(art, ax=axes[1,2])
    fig.savefig(OUTPUT_DIR / "phaseless_reconstruction.png", bbox_inches="tight")
    plt.show()

    reconstruction_summary = pd.DataFrame([{
        "true_separation_mas": SOURCE.separation_mas,
        "recovered_separation_mas": result.metrics.get("two_peak_separation_mas"),
        "true_position_angle_deg": SOURCE.position_angle_deg,
        "recovered_position_angle_deg": result.metrics.get("two_peak_position_angle_deg"),
        "true_flux_ratio": SOURCE.flux_ratio_secondary_to_primary,
        "recovered_peak_ratio": result.metrics.get("two_peak_brightness_ratio"),
        "weighted_forward_rmse": result.metrics["weighted_fit_rmse"],
        "raw_truth_correlation": truth_corr_raw,
        "common_beam_truth_correlation": truth_corr,
        "common_beam_truth_nrmse": truth_nrmse,
        "used_180_degree_mirror_for_validation": used_mirror,
    }])
    display(reconstruction_summary.T.rename(columns={0:"value"}))
    reconstruction_summary.to_csv(OUTPUT_DIR / "reconstruction_summary.csv", index=False)

    np.save(OUTPUT_DIR / "reconstruction_image.npy", result.image)
    pd.DataFrame({
        "theta_x_mas": np.tile(result.theta_mas, len(result.theta_mas)),
        "theta_y_mas": np.repeat(result.theta_mas, len(result.theta_mas)),
        "brightness": result.image.ravel(),
    }).to_csv(OUTPUT_DIR / "reconstruction_image.csv", index=False)
    """),
    md(r"""
    ## 9. 几种情况：能看到多少星等，角分辨率是多少

    这里区分三个概念：

    1. `fringe_scale = λ/Bmax`：最长投影基线对应的名义角尺度；
    2. `uniform_disk_first_null = 1.22 λ/Bmax`：均匀圆盘第一零点可被最长基线采到时的直径尺度；
    3. 星等极限：在未解析源 (|V|^2=1) 上达到指定 SNR 的 AB 星等。它随源被分辨后的 (|V|^2)、NSB、系统效率和系统误差恶化。

    光谱复用假设各通道具有相同电子带宽且统计独立，因此总 SNR 按 (sqrt{N_{ch}}) 提升。这是硬件情景，不是免费增益。
    """),
    code(r"""
    def limiting_magnitude(target_snr, hours, spectral_channels=1, nsb_rate_hz=None):
        return brentq(
            lambda mag: unit_visibility_snr(mag, hours*3600, INST,
                                            spectral_channels, nsb_rate_hz)-target_snr,
            -8.0, 18.0,
        )

    sensitivity_rows = []
    for nsb_label, nsb_rate in [("ideal dark (0 MHz)",0.0), ("reference NSB (70.527 MHz)",70.527e6)]:
        for hours in [1,10,50]:
            for channels in [1,16]:
                sensitivity_rows.append({
                    "NSB_case": nsb_label, "hours":hours, "spectral_channels":channels,
                    "m_AB_at_SNR5":limiting_magnitude(5,hours,channels,nsb_rate),
                    "m_AB_at_SNR10":limiting_magnitude(10,hours,channels,nsb_rate),
                })
    sensitivity_table = pd.DataFrame(sensitivity_rows)

    resolution_rows = []
    for name, baselines in [("branch_7",branch_baselines),("layout_0803_reco32",lact_baselines)]:
        bmax = baselines.baseline_m.max()
        positive = baselines.loc[baselines.baseline_m>0,"baseline_m"]
        bmin = positive.min()
        resolution_rows.append({
            "layout":name, "Bmin_m":bmin, "Bmax_m":bmax,
            "fringe_scale_mas":INST.wavelength_m/bmax*RAD_TO_MAS,
            "uniform_disk_first_null_mas":1.22*INST.wavelength_m/bmax*RAD_TO_MAS,
            "largest_scale_lambda_over_Bmin_mas":INST.wavelength_m/bmin*RAD_TO_MAS,
        })
    resolution_table = pd.DataFrame(resolution_rows)
    display(Markdown("### 星等极限（单条未解析基线）"))
    display(sensitivity_table.round(3))
    display(Markdown("### 角尺度"))
    display(resolution_table.round(4))

    fig, axes = plt.subplots(1,2,figsize=(11.4,3.8))
    reference = sensitivity_table[sensitivity_table.NSB_case.str.startswith("reference")]
    for channels,color,marker in [(1,BLUE,"o"),(16,ORANGE,"s")]:
        subset=reference[reference.spectral_channels==channels]
        axes[0].plot(subset.hours,subset.m_AB_at_SNR5,marker=marker,color=color,label=f"{channels} spectral channel(s)")
    axes[0].set_xscale("log")
    axes[0].invert_yaxis()
    axes[0].set(title="Limiting AB magnitude at SNR=5",xlabel="Integration [h]",ylabel="Faint limit [mag; fainter downward]")
    axes[0].legend(frameon=False)
    axes[0].grid(alpha=.2)
    axes[1].bar(resolution_table.layout,resolution_table.fringe_scale_mas,color=[GOLD,BLUE])
    axes[1].set_yscale("log")
    axes[1].set(title="Nominal angular scale at 400 nm",ylabel=r"$\lambda/B_{max}$ [mas]")
    for i,val in enumerate(resolution_table.fringe_scale_mas): axes[1].text(i,val*1.08,f"{val:.3f}",ha="center")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR/"sensitivity_and_resolution.png",bbox_inches="tight")
    plt.show()

    sensitivity_table.to_csv(OUTPUT_DIR/"limiting_magnitude_cases.csv",index=False)
    resolution_table.to_csv(OUTPUT_DIR/"angular_resolution_cases.csv",index=False)
    """),
    md(r"""
    ## 10. Checks

    下列检查是本 notebook 的最低物理闭合门槛：

    - (V(0,0)=1) 且 (|V(u,v)|^2=|V(-u,-v)|^2)；
    - 32 镜编号连续、无重复，且 NWU→ENU 坐标转换逐行一致；
    - 每条基线旋转后满足 (u^2+v^2+w^2=|\mathbf B|^2)；
    - 5 星等的光子率严格为 0 星等的 1%；
    - SNR 对积分时间按 (sqrt T) 缩放；
    - 重建只用观测功率，前向闭合 RMSE 有界，并能在固有歧义下恢复双源结构。
    """),
    code(r"""
    checks = {}
    checks["layout_has_32_unique_telescopes"] = (
        len(LACT32) == 32 and LACT32.telescope_id.is_unique and LACT32.lactsim_index.is_unique
    )
    checks["layout_nwu_to_enu_transform"] = (
        np.allclose(LACT32.east_m, -LACT32.corsika_west_cm/100.0)
        and np.allclose(LACT32.north_m, LACT32.corsika_north_cm/100.0)
        and np.allclose(LACT32.up_m, LACT32.corsika_up_cm/100.0)
    )
    checks["zero_baseline_visibility"] = np.isclose(abs(binary_visibility(0,0))**2,1,atol=1e-12)
    probe_u = RNG.normal(0,1e9,100)
    probe_v = RNG.normal(0,1e9,100)
    checks["power_conjugate_symmetry"] = np.allclose(
        np.abs(binary_visibility(probe_u,probe_v))**2,
        np.abs(binary_visibility(-probe_u,-probe_v))**2,atol=1e-12)
    checks["uvw_norm_invariant"] = invariant_error < 1e-6
    checks["five_mag_rate_ratio"] = np.isclose(detected_star_rate_hz(5)/detected_star_rate_hz(0),0.01,rtol=1e-12)
    checks["snr_sqrt_time"] = np.isclose(
        unit_visibility_snr(2,3600)/unit_visibility_snr(2,900),2.0,rtol=1e-12)
    checks["short_waveform_has_expected_samples"] = len(short_record["sample_time_ns"]) in (106, 107)
    checks["recovery_loss_is_tiny_at_reference_rate"] = (
        recovery_table.analytic_mean_charge_fraction.min() > 0.9999
    )
    checks["reconstruction_forward_closure"] = result.metrics["weighted_fit_rmse"] < 0.12
    checks["reconstruction_recovers_two_peaks"] = len(result.metrics.get("peaks",[])) == 2
    checks["truth_validation_correlation"] = truth_corr > 0.55
    validation = pd.DataFrame({"check":checks.keys(),"passed":checks.values()})
    display(validation)
    assert all(checks.values()), checks
    validation.to_csv(OUTPUT_DIR/"validation_checks.csv",index=False)
    print("ALL CHECKS PASSED")
    """),
    md(r"""
    ## 11. Takeaways & required caveats

    1. 统一版以最新 `main` 为底座，保留旧 v1 的全阵列源/相干/重建框架，并吸收旧 v2 的连续电子学时间尺度；旧 v2 把 ENU 坐标直接代入缺少站点纬度的矩阵，统一版已经用 topocentric ENU 投影替换。
    2. 7 镜验证阵列适合标定和直径测量；`layout_0803` 的 32 镜实际模拟坐标提供公里级基线，进入约 (0.1\,\mathrm{mas}) 及以下的成像区间。表中有物理最长基线和本次时角/赤纬的投影覆盖；它仍是生产输入坐标，不应冒充最终现场测绘值。
    3. 单通道 6 m 望远镜的星等极限主要由口径×效率、电子带宽、NSB 与时间决定。对已被分辨的目标，实际阈值还要减去 (|V|^2) 带来的 SNR；“能检测星”不等于“能重建复杂表面”。
    4. (|V|^2) 成像天然缺相位。当前 MAP 重建适合闭合研究，但正式论文还应做多次噪声 realization、超参数证据/交叉验证、bootstrap 图像不确定度，并与 MiRA/SQUEEZE 或 HIO/ER 相位恢复交叉比较。
    5. 微单元恢复现已按指数模型进入 C++ 电子学链和短波形验证；10 ns 只是没有 S17351 公开数据时的暂定值，并已扫描 1/10/30 ns。当前仍未包括 LACT 实测镜面到达时间展宽、串扰、后脉冲、通带随角度变化、完整复传递函数、时钟漂移和标定星系统误差；这些是从“研究预测”升级到“仪器性能声明”的必要输入。
    """),
]

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
nbf.write(nb, OUTPUT)
print(OUTPUT)
