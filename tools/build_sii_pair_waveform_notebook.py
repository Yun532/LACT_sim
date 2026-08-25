#!/usr/bin/env python3
"""生成双望远镜 HBT 逐 p.e. 到完整波形的中文验证 notebook。"""

from pathlib import Path
import textwrap

import nbformat as nbf


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "notebooks" / "lact_sii_two_telescope_waveforms.ipynb"


def md(text):
    return nbf.v4.new_markdown_cell(textwrap.dedent(text).strip())


def code(text):
    return nbf.v4.new_code_cell(textwrap.dedent(text).strip())


nb = nbf.v4.new_notebook()
nb["metadata"]["kernelspec"] = {
    "display_name": "Python 3", "language": "python", "name": "python3"}
nb["cells"] = [
    md(r"""
    # 两台 LACT：从热光光子到两路真实波形和 HBT 互相关

    ## 结论先行

    本 notebook 不把相关性人为放大。它生成恒星、NSB 和稀疏 HBT 相关光子对，
    再依次通过 LACT 光线追迹得到的时间核、S17351 微单元恢复、实测单 p.e. 电荷、
    SPE 波形、加性噪声与 ADC。最后只用两路波形计算互相关。

    执行结果：mAB=2、$|V|^2=0.5$ 的 200 µs 记录只预期0.512个相关对，本次抽到1个；
    原始互相关仍看不到6 ns处的HBT峰。main当前实测SPE的自相关FWHM约47.5 ns，
    远宽于两镜DC相对展宽0.852 ns，因此DC只让实测波形相关峰高度再下降约0.5%。
    简单零滞后的等效带宽为14.46 MHz；用实测传递函数做频域匹配后，原生4 ns数据
    可达87.98 MHz。相比旧200 MHz矩形带宽模型，推荐匹配算法的SNR约为79.7%、
    所需时间约为1.58倍；简单零滞后则约需9.59倍时间。
    """),
    md("## 1. 设置和 main 参数"),
    code("""
    from dataclasses import replace
    from pathlib import Path
    import json, shutil, subprocess, sys
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from IPython.display import display
    from scipy.signal import correlate

    REPO_ROOT = Path.cwd().resolve()
    if not (REPO_ROOT / "python" / "sii_unified.py").exists():
        REPO_ROOT = REPO_ROOT.parent
    sys.path.insert(0, str(REPO_ROOT / "python"))
    from sii_unified import (
        BinarySource, Instrument, Observation, detected_star_rate_hz,
        generate_uvw, hbt_correlated_pair_rate_hz,
        load_measured_spe_template, load_optical_timing_mixture,
        optical_timing_transfer_efficiency, render_pe_waveform,
        sample_optical_delays_ns,
        simulate_hbt_primary_pe, simulate_uv_observation, unit_visibility_snr,
        waveform_cross_correlation, write_main_primary_pe_csv,
    )

    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'Noto Sans SC', 'SimHei']
    plt.rcParams['axes.unicode_minus'] = False
    OUTPUT_DIR = REPO_ROOT / "run_logs" / "sii_pair_waveform"
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    instrument = Instrument.from_repository(REPO_ROOT)
    timing = load_optical_timing_mixture(instrument.optical_timing_kernel_path)
    measured_template = load_measured_spe_template(instrument.spe_template_path)
    """),
    md(r"""
    ## 2. 正确的相关光子流

    当光学相干时间远小于电子学分辨率时，热光的二阶超额相关可以等价写成稀疏相关
    光子对：

    $$R_{\rm pair}=R_{\star,1}R_{\star,2}\tau_c|V|^2.$$

    其余恒星光子和 NSB 都是两镜独立 Poisson 流。相关光子不是“复制全部波形”，
    而只是总光子流里极少的一部分。每个光子随后独立抽取各自望远镜的光程延迟。
    """),
    code("""
    magnitude = 2.0
    visibility2 = 0.50
    geometric_delay_ns = 6.0
    duration_ns = 200_000.0       # 200 us：完整保存但仍远短于天文积分时间
    star_rate_hz = detected_star_rate_hz(magnitude, instrument)
    nsb_rate_hz = instrument.detected_nsb_rate_hz
    pair_rate_hz = hbt_correlated_pair_rate_hz(
        star_rate_hz, star_rate_hz, instrument.coherence_area_s, visibility2)

    hits, stream_info = simulate_hbt_primary_pe(
        np.random.default_rng(20260825), duration_ns,
        star_rate_hz, nsb_rate_hz, visibility2, instrument,
        optical_timing_mixture=timing,
        geometric_delay_ns=geometric_delay_ns, padding_ns=220.0)
    primary_csv = write_main_primary_pe_csv(
        hits, OUTPUT_DIR / "two_telescope_primary_pe.csv")
    display(pd.Series(stream_info, name="value").to_frame())
    print(f"两镜逐 p.e. 行数：{len(hits):,}")
    print(f"main 原生输入：{primary_csv}")
    """),
    md(r"""
    对当前 $m_{AB}=2$、$|V|^2=0.5$，200 µs 中预期的 HBT 对仍远少于随机恒星和
    NSB 光子。因此下面两条完整波形看起来几乎独立，这是天文强度干涉的正常尺度。
    """),
    md("## 3. 两台望远镜：main 当前实测 SPE 与原生 4 ns 波形"),
    code("""
    native_sample_ns = instrument.sample_width_ns
    waveforms_measured = []
    for telescope_id in (0, 1):
        times = hits.loc[hits.telescope_id == telescope_id, "time_ns"].to_numpy()
        waveforms_measured.append(render_pe_waveform(
            np.random.default_rng(3100+telescope_id), times, duration_ns,
            instrument, sample_width_ns=native_sample_ns,
            template=measured_template))

    # 完整记录是 200 us；这里明确画 2 us（2000 ns），不是只模拟 400 ns。
    view_start, view_end = 20_000.0, 22_000.0
    fig, axes = plt.subplots(2, 1, figsize=(12, 5.8), sharex=True, constrained_layout=True)
    for telescope_id, (axis, waveform) in enumerate(zip(axes, waveforms_measured)):
        keep = ((waveform["sample_time_ns"] >= view_start) &
                (waveform["sample_time_ns"] <= view_end))
        axis.plot(waveform["sample_time_ns"][keep]-view_start,
                  waveform["adc_mv"][keep], lw=.8)
        axis.set(ylabel="ADC [mV]", title=f"望远镜 {telescope_id+1}：实测SPE + 4 ns + NSB + 噪声 + ADC")
    axes[-1].set(xlabel="局部时间 [ns]（本图连续显示 2000 ns）")
    fig.savefig(OUTPUT_DIR/"two_measured_spe_waveforms.png", dpi=160)
    plt.show()
    print({"spe_file": instrument.spe_template_path,
           "spe_support_ns": [measured_template[0].min(), measured_template[0].max()],
           "spe_peak_mv": measured_template[1].max(),
           "spe_area_mv_ns": np.trapezoid(measured_template[1], measured_template[0]),
           "native_sample_ns": native_sample_ns})
    """),
    md(r"""
    ## 4. 只从两路波形计算互相关

    下面没有读取 `hbt_pair_id`，只对两条 ADC 波形去均值后做 FFT 互相关。
    由于第二路晚到 6 ns，按这里的顺序峰应位于 $+6$ ns；但单个 200 µs 记录的
    统计噪声远大于 HBT 超额相关，所以不能把最高的随机尖峰误认成天文信号。
    """),
    code("""
    lag_measured, corr_measured = waveform_cross_correlation(
        waveforms_measured[1]["adc_mv"], waveforms_measured[0]["adc_mv"],
        native_sample_ns, 200.0)
    raw_peak_lag = lag_measured[np.argmax(corr_measured)]
    raw_peak_value = corr_measured.max()
    raw_rms = np.std(corr_measured[
        np.abs(lag_measured-geometric_delay_ns) > 60.0])
    fig, ax = plt.subplots(figsize=(9, 3.7), constrained_layout=True)
    ax.plot(lag_measured, corr_measured, lw=1)
    ax.axvline(geometric_delay_ns, color="crimson", ls="--", label="已知几何时延")
    ax.set(title="200 µs 两路实测SPE/4 ns ADC波形互相关：仍由随机噪声主导",
           xlabel="B 相对 A 的滞后 [ns]", ylabel="相关系数")
    ax.legend()
    fig.savefig(OUTPUT_DIR/"measured_spe_raw_cross_correlation.png", dpi=160)
    plt.show()
    print({"raw_max_lag_ns": raw_peak_lag,
           "raw_max_correlation": raw_peak_value,
           "off_peak_rms": raw_rms})
    """),
    md(r"""
    ## 5. HBT 波形相关的期望形状：DC 前后

    为了看清被随机背景淹没的期望值，这里只计算“一个相关光子对经过两路完整脉冲
    响应后的平均互相关”，然后乘以实际的 $R_{\rm pair}$ 即得到物理协方差。
    这不是放大后的假观测，而是对微弱信号项做重要抽样。
    """),
    code("""
    response_dt_ns = 0.05
    response_lags = np.arange(-250.0, 250.0+response_dt_ns/2, response_dt_ns)
    rng_response = np.random.default_rng(52)
    response_pairs = 300_000
    relative_delay_dc = (geometric_delay_ns
        + sample_optical_delays_ns(rng_response, response_pairs, timing)
        - sample_optical_delays_ns(rng_response, response_pairs, timing))
    relative_delay_iso = np.full(response_pairs, geometric_delay_ns)

    def expected_pair_waveform_correlation(relative_delays):
        edges = np.r_[response_lags-response_dt_ns/2,
                      response_lags[-1]+response_dt_ns/2]
        delay_hist = np.histogram(relative_delays, bins=edges)[0].astype(float)
        delay_hist /= delay_hist.sum()
        pulse = np.interp(response_lags, measured_template[0],
                          measured_template[1], left=0.0, right=0.0)
        pulse_ac = correlate(pulse, pulse, mode="same", method="fft")
        # 以零滞后为中心的循环卷积；范围远大于脉冲和时间核支撑。
        return np.fft.fftshift(np.fft.ifft(
            np.fft.fft(np.fft.ifftshift(delay_hist))
            * np.fft.fft(np.fft.ifftshift(pulse_ac))).real)

    expected_iso = expected_pair_waveform_correlation(relative_delay_iso)
    expected_dc = expected_pair_waveform_correlation(relative_delay_dc)
    normalization = expected_iso.max()
    expected_iso /= normalization
    expected_dc /= normalization
    dc_waveform_peak_retention = expected_dc.max()
    pair_delay_rms_ns = np.std(relative_delay_dc-geometric_delay_ns)
    iso_half = response_lags[expected_iso >= 0.5]
    dc_half = response_lags[expected_dc >= 0.5*expected_dc.max()]
    measured_spe_ac_fwhm_ns = iso_half[-1]-iso_half[0]
    measured_spe_dc_fwhm_ns = dc_half[-1]-dc_half[0]

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 3.8), constrained_layout=True)
    axes[0].hist(relative_delay_dc-geometric_delay_ns, bins=120, density=True,
                 alpha=.75, label="两镜独立 DC 时间核之差")
    axes[0].axvline(0, color='.2', ls='--')
    axes[0].set(title="相关光子对的相对到达时间", xlabel="去除几何时延后的差 [ns]",
                ylabel="概率密度")
    axes[1].plot(response_lags, expected_iso, label="等时镜面")
    axes[1].plot(response_lags, expected_dc, label="LACT 时间核")
    axes[1].axvline(geometric_delay_ns, color='.2', ls='--')
    axes[1].set(xlim=(-80,92), title="main 实测 SPE 后的 HBT 期望相关峰",
                xlabel="B 相对 A 的滞后 [ns]", ylabel="相对协方差")
    axes[1].legend()
    fig.savefig(OUTPUT_DIR/"expected_hbt_peak_dc_comparison.png", dpi=160)
    plt.show()
    print({"single_telescope_timing_rms_ns": timing["rms_spread_ns"],
           "pair_delay_rms_ns": pair_delay_rms_ns,
           "measured_SPE_peak_retention_after_DC": dc_waveform_peak_retention,
           "SPE_autocorrelation_FWHM_ns": measured_spe_ac_fwhm_ns,
           "SPE_plus_DC_FWHM_ns": measured_spe_dc_fwhm_ns})
    """),
    md(r"""
    实测SPE的自相关本身约几十ns宽，已经远宽于0.852 ns两镜光学展宽。因此在这条
    实测波形上，DC只造成很小的额外峰高损失。这个结果不表示DC没有影响，而是说明
    当前电子脉冲响应比DC更慢，系统的有效相关带宽主要由SPE与相关算法决定。
    """),
    md("## 6. 实测 SPE 的有效带宽：简单零滞后与匹配相关"),
    code("""
    fine_dt_ns = 0.01
    fine_time_ns = np.arange(measured_template[0].min(),
                             measured_template[0].max()+fine_dt_ns/2,
                             fine_dt_ns)
    fine_pulse = np.interp(fine_time_ns, measured_template[0],
                           measured_template[1])
    fine_dt_s = fine_dt_ns*1e-9
    # 零填充只提高频率插值精度，不改变脉冲或带宽积分。
    fft_length = 1 << int(np.ceil(np.log2(len(fine_pulse)*16)))
    full_frequency_hz = np.fft.rfftfreq(fft_length, fine_dt_s)
    full_pulse_spectrum = np.fft.rfft(
        fine_pulse, n=fft_length)*fine_dt_s
    relevant_frequency = full_frequency_hz <= 2.5e9
    frequency_hz = full_frequency_hz[relevant_frequency]
    pulse_spectrum = full_pulse_spectrum[relevant_frequency]

    weights = np.asarray(timing["weights"])
    means = np.asarray(timing["mean_delay_ns"])
    sigmas = np.asarray(timing["std_delay_ns"])
    optical_transfer = (
        np.exp(-2j*np.pi*(frequency_hz*1e-9)[:,None]*means[None,:])
        * np.exp(-2*np.pi**2*(frequency_hz*1e-9)[:,None]**2
                 * sigmas[None,:]**2)) @ weights
    optical_power = np.abs(optical_transfer)**2

    pulse_power = np.abs(pulse_spectrum)**2
    direct_noise_integral = np.trapezoid(pulse_power**2, frequency_hz)
    direct_signal_no_dc = np.trapezoid(pulse_power, frequency_hz)
    direct_signal_dc = np.trapezoid(
        pulse_power*optical_power, frequency_hz)
    direct_bandwidth_no_dc_hz = (
        direct_signal_no_dc**2/direct_noise_integral)
    direct_bandwidth_dc_hz = direct_signal_dc**2/direct_noise_integral

    adc_step_mv = instrument.adc_full_scale_mv/(2**instrument.adc_bits)
    sample_noise_variance_mv2 = (
        instrument.electronic_noise_rms_mv**2 + adc_step_mv**2/12)
    total_rate_hz = star_rate_hz+nsb_rate_hz

    def matched_bandwidth(sample_width_ns):
        sample_rate_hz = 1e9/sample_width_ns
        # main每个样点保存采样窗内平均电压，故有sinc采样窗传递函数。
        sampled_pulse = pulse_spectrum*np.sinc(frequency_hz/sample_rate_hz)
        shot_psd = (2*total_rate_hz*instrument.excess_noise_factor**2
                    * np.abs(sampled_pulse)**2)
        white_noise_psd = 2*sample_noise_variance_mv2/sample_rate_hz
        shot_fraction = shot_psd/(shot_psd+white_noise_psd)
        keep = frequency_hz <= sample_rate_hz/2
        return np.trapezoid(
            optical_power[keep]**2*shot_fraction[keep]**2,
            frequency_hz[keep])

    matched_4ns_hz = matched_bandwidth(native_sample_ns)
    matched_1p6ns_hz = matched_bandwidth(1.6)
    matched_0p25ns_hz = matched_bandwidth(0.25)
    assumed_200mhz_equivalent_hz = (
        200e6*optical_timing_transfer_efficiency(timing, 200e6)**2)
    bandwidth_table = pd.DataFrame([
        {"estimator":"实测SPE，简单零滞后", "effective_bandwidth_MHz":direct_bandwidth_dc_hz/1e6},
        {"estimator":"实测SPE，4 ns匹配频谱", "effective_bandwidth_MHz":matched_4ns_hz/1e6},
        {"estimator":"实测SPE，1.6 ns匹配频谱", "effective_bandwidth_MHz":matched_1p6ns_hz/1e6},
        {"estimator":"实测SPE，0.25 ns匹配频谱", "effective_bandwidth_MHz":matched_0p25ns_hz/1e6},
        {"estimator":"旧模型：200 MHz矩形带宽+DC", "effective_bandwidth_MHz":assumed_200mhz_equivalent_hz/1e6},
    ])
    bandwidth_table["SNR_relative_to_old_model"] = np.sqrt(
        bandwidth_table.effective_bandwidth_MHz
        /(assumed_200mhz_equivalent_hz/1e6))
    bandwidth_table["time_factor_vs_old_model"] = (
        1/bandwidth_table.SNR_relative_to_old_model**2)
    display(bandwidth_table)

    half_power_index = np.flatnonzero(
        pulse_power/pulse_power[0] <= 0.5)[0]
    spe_half_power_mhz = frequency_hz[half_power_index]/1e6
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 3.9), constrained_layout=True)
    axes[0].plot(measured_template[0], measured_template[1])
    axes[0].set(title="main 当前实测单p.e.波形", xlabel="相对到达时间 [ns]",
                ylabel="幅度 [mV]")
    axes[1].semilogx(frequency_hz[1:]/1e6,
                     pulse_power[1:]/pulse_power[0], label="|P(f)|²")
    axes[1].axvline(spe_half_power_mhz, color='crimson', ls='--',
                    label=f"半功率 {spe_half_power_mhz:.1f} MHz")
    axes[1].axvline(125, color='.4', ls=':', label="4 ns Nyquist")
    axes[1].set(xlim=(1,500), ylim=(1e-5,1.2), yscale='log',
                title="实测SPE功率频谱", xlabel="频率 [MHz]", ylabel="归一功率")
    axes[1].legend()
    fig.savefig(OUTPUT_DIR/"measured_spe_and_bandwidth.png", dpi=160)
    plt.show()
    """),
    md(r"""
    这解释了为什么“170 ns长尾”和“4 ns采样”不能只看一个数字：直接取零滞后时，
    实测SPE只给出十几MHz量级的等效带宽；若在频域用实测传递函数对白噪声和光子
    shot noise做匹配加权，可以恢复到接近百MHz。后者是建议的真实相关器算法，但其
    数值仍依赖当前尚未实测的0.35 mV白噪声假设和ADC模型。
    """),
    md("## 7. 为什么短波形看不见：理论 SNR 对照"),
    code("""
    def baseline_snr(integration_s, effective_bandwidth_hz, magnitude_value=magnitude):
        trial = replace(
            instrument, electronics_bandwidth_hz=effective_bandwidth_hz)
        return (unit_visibility_snr(magnitude_value, integration_s, trial)
                * visibility2)

    # effective bandwidth已经包含DC、SPE、采样、电子噪声和ADC传递损失。
    def five_sigma_time(effective_bandwidth_hz, magnitude_value):
        return (5.0/baseline_snr(
            1.0, effective_bandwidth_hz, magnitude_value))**2

    snr_record = baseline_snr(
        duration_ns*1e-9, matched_4ns_hz)
    snr_table = pd.DataFrame([
        {"mAB":mag, "estimator":"简单零滞后",
         "effective_bandwidth_MHz":direct_bandwidth_dc_hz/1e6,
         "hours_for_5sigma":five_sigma_time(direct_bandwidth_dc_hz, mag)/3600}
        for mag in (2.0, 3.0)
    ] + [
        {"mAB":mag, "estimator":"4 ns匹配频谱",
         "effective_bandwidth_MHz":matched_4ns_hz/1e6,
         "hours_for_5sigma":five_sigma_time(matched_4ns_hz, mag)/3600}
        for mag in (2.0, 3.0)
    ] + [
        {"mAB":mag, "estimator":"旧200 MHz模型",
         "effective_bandwidth_MHz":assumed_200mhz_equivalent_hz/1e6,
         "hours_for_5sigma":five_sigma_time(
             assumed_200mhz_equivalent_hz, mag)/3600}
        for mag in (2.0, 3.0)
    ])
    display(snr_table)
    print({"200us_matched_SNR": snr_record,
           "4ns_matched_vs_old_SNR": np.sqrt(
               matched_4ns_hz/assumed_200mhz_equivalent_hz),
           "4ns_matched_time_factor_vs_old":
               assumed_200mhz_equivalent_hz/matched_4ns_hz})
    print("短记录 HBT SNR / 原始互相关随机峰 SNR =",
          snr_record/(raw_peak_value/raw_rms))
    """),
    md(r"""
    ## 8. `main` 完整光学与当前时间核：哪些相同，哪些不同

    当前默认时间核来自一次 `run_optical_sim` 的 590,239 条轴上 400 nm 焦面命中，
    因而**光程到达时间分布是真实光线追迹的压缩结果**。但是旧运行使用理想误差、
    白板焦面，没有完整保留 LACT2 实测逐镜曲率/偏转、结构遮挡、真实 1656 像素、
    集光器和 PDE。下面的新配置会用 `main` 的完整链重新生成响应；若已经在同一仓库
    运行并提取，provenance 会自动显示真实结果，否则明确显示“待运行”，不会拿旧核冒充。
    """),
    code("""
    full_cfg = REPO_ROOT/"configs/optics/lact2_measured_full_response_400nm.cfg"
    full_kernel = REPO_ROOT/"configs/optics/lact2_measured_full_response_400nm.csv"
    full_provenance = REPO_ROOT/"configs/optics/lact2_measured_full_response_400nm.provenance.json"
    old_provenance = json.loads((REPO_ROOT/"configs/optics/lact_1229_onaxis_timing_kernel.provenance.json").read_text(encoding="utf-8"))
    optical_audit = pd.DataFrame([
        {"项目":"逐镜几何与光程", "当前默认核":"有：main完整追迹后压缩", "新完整响应":"有"},
        {"项目":"LACT2实测镜面/仰角形变", "当前默认核":"无：旧理想误差", "新完整响应":"有：70°实测插值"},
        {"项目":"3D支架遮挡", "当前默认核":"无", "新完整响应":"有"},
        {"项目":"真实相机/集光器/PDE", "当前默认核":"无", "新完整响应":"有"},
        {"项目":"实测SPE/4 ns/ADC", "当前波形":"有", "新完整响应":"光学后由电子学独立处理"},
        {"项目":"HBT相关性的来源", "当前与新流程":"入射热光统计", "新完整响应":"光学本身不会产生HBT"},
    ])
    display(optical_audit)
    if full_provenance.exists() and full_kernel.exists():
        full_info = json.loads(full_provenance.read_text(encoding="utf-8"))
        display(pd.Series(full_info, name="完整main光学响应").to_frame())
    else:
        full_info = None
        print("待运行：build/run_optical_sim", full_cfg)
        print("然后运行：python tools/derive_full_optical_response.py <hits.csv> <kernel.csv> --input-photons 1000000 --provenance-json <json> --source-config", full_cfg)
    """),
    md(r"""
    ## 9. 完整响应对 UV 平面的真正影响

    光学和电子学**不会移动**基线的 $(u,v)$ 坐标，也不会改变源模型给出的真实
    $|V(u,v)|^2$。它们改变的是每个 UV 点的误差、权重、系统偏差和最终能否重建。
    下面固定同一套32镜坐标、同一随机种子和2小时观测，只替换相关器等效带宽。
    """),
    code("""
    layout_path = REPO_ROOT/"configs/arrays/layout_0803_reco32_coordinates.csv"
    uv_source = BinarySource(ab_magnitude=3.0)
    uv_observation = Observation(hours_per_night=2.0, nights=1, segment_s=1200.0)
    uvw_real = generate_uvw(layout_path, uv_observation, instrument)

    def effective_instrument(bandwidth_hz, label):
        # 带宽已经包含光学核、SPE、采样、加性噪声和ADC传递；关闭这些项的二次扣除。
        return replace(instrument, electronics_bandwidth_hz=bandwidth_hz,
                       optical_timing_kernel_path=None, spe_template_path=None,
                       parameter_source=label)

    uv_cases = [
        ("旧200 MHz假设", assumed_200mhz_equivalent_hz),
        ("实测SPE + 4 ns匹配", matched_4ns_hz),
        ("实测SPE + 简单零滞后", direct_bandwidth_dc_hz),
    ]
    uv_results = []
    for label, bandwidth in uv_cases:
        measured, meta = simulate_uv_observation(
            uvw_real, uv_source, uv_observation,
            effective_instrument(bandwidth, label), seed=9102)
        uv_results.append((label, measured, meta))

    common_limit = max(np.quantile(np.abs(item.visibility2_measured-item.visibility2_true), .99)
                       for _, item, _ in uv_results)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2), sharex=True, sharey=True,
                             constrained_layout=True)
    for axis, (label, measured, meta) in zip(axes, uv_results):
        residual = measured.visibility2_measured-measured.visibility2_true
        points = axis.scatter(measured.u_m, measured.v_m, c=residual, s=8,
                              cmap="coolwarm", vmin=-common_limit, vmax=common_limit)
        axis.set(title=f"{label}\\n中位σ={meta['sigma_visibility2_stat']:.3g}",
                 xlabel="u [m]", ylabel="v [m]")
        axis.set_aspect("equal", adjustable="box")
    fig.colorbar(points, ax=axes, label="模拟值 - 真实 |V|²")
    fig.savefig(OUTPUT_DIR/"uv_noise_full_response_comparison.png", dpi=160)
    plt.show()
    display(pd.DataFrame([{
        "情况": label, "等效带宽_MHz": bandwidth/1e6,
        "UV点数": meta["uv_measurements"],
        "中位统计误差": meta["sigma_visibility2_stat"],
        "总积分_h": meta["total_integration_hours"],
    } for (label, bandwidth), (_, _, meta) in zip(uv_cases, uv_results)]))
    """),
    md("## 10. 可选：把同一批光子直接交给 main 电子学"),
    code("""
    preview = hits[(hits.time_ns >= -220) & (hits.time_ns <= 2_220)].copy()
    preview_csv = write_main_primary_pe_csv(
        preview, OUTPUT_DIR/"two_telescope_primary_pe_main_preview.csv")
    executable_candidates = [
        REPO_ROOT/"build"/"run_camera_electronics",
        REPO_ROOT/"build"/"run_camera_electronics.exe"]
    executable = next((path for path in executable_candidates if path.exists()), None)
    main_output = OUTPUT_DIR/"main_electronics_preview"
    main_config = Path(instrument.parameter_source)
    command = ([str(executable)] if executable else ["build/run_camera_electronics"])
    command += [str(main_config), str(preview_csv), str(main_output),
                "-C", "electronics.n_pixels=1",
                "-C", "electronics.sampling.start_ns=0",
                "-C", "electronics.sampling.end_ns=2000"]
    print(" ".join(command))
    if executable:
        completed = subprocess.run(command, cwd=REPO_ROOT, check=True,
                                   text=True, capture_output=True)
        print(completed.stdout)
    else:
        print("当前工作区没有已编译的 main 可执行文件；上面的 CSV 和命令可在服务器 build 后直接运行。")
    """),
    md(r"""
    ## 11. 检查与边界

    - 两镜相关对数必须服从 $R_1R_2\tau_c|V|^2$，而不是共享全部 Poisson 光子。
    - 两个独立 LACT 时间核之差的 RMS 应接近 $\sqrt{2}\times0.602$ ns。
    - 200 µs 原始波形不应显著探测到天文 HBT；否则很可能把同源噪声或模拟器共享随机数
      误当成了恒星相关。
    - 本 notebook 确实读取 main 当前 `single_pe.template` 指向的实测CSV；频域匹配结果
      使用0.35 mV独立白噪声和8-bit均匀量化假设，这两项尚不是实测噪声功率谱。
    - 当前时间核仅为 400 nm、轴上、理想误差配置。论文级性能仍需要逐镜、随仰角/离轴角
      的实测时间核，以及实测快速前端传递函数与跨镜时钟稳定性。

    方法依据：[H.E.S.S. 两镜实测与相干时间定义](https://academic.oup.com/mnras/article/527/4/12243/7455900)、
    [时间分辨率、光子率和 IACT 非等时性的比较](https://academic.oup.com/mnras/article/512/2/1722/6534919)、
    [有限仪器响应下相关峰面积与 SNR](https://academic.oup.com/mnras/article/537/3/2527/8003771)。
    """),
    code("""
    checks = {
        "pair_rate_formula": np.isclose(pair_rate_hz,
            star_rate_hz**2*instrument.coherence_area_s*visibility2),
        "pair_delay_rms_is_sqrt2_single": np.isclose(
            pair_delay_rms_ns, np.sqrt(2)*timing["rms_spread_ns"], rtol=.02),
        "short_raw_record_is_below_one_sigma": snr_record < 1.0,
        "measured_spe_is_current_main_file": (
            Path(instrument.spe_template_path).name == "spe_template_measured.csv"),
        "measured_spe_dc_peak_is_physical": (
            0.0 < dc_waveform_peak_retention < 1.0),
        "matched_bandwidth_exceeds_raw_zero_lag": (
            matched_4ns_hz > direct_bandwidth_dc_hz),
        "old_model_is_more_optimistic_than_measured_4ns": (
            assumed_200mhz_equivalent_hz > matched_4ns_hz),
        "main_csv_exists": primary_csv.exists(),
        "waveforms_are_finite": all(np.isfinite(w["adc_mv"]).all()
                                     for w in waveforms_measured),
    }
    display(pd.Series(checks, name="passed").to_frame())
    assert all(checks.values())
    print("全部物理和程序检查通过。")
    """),
]

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
nbf.write(nb, OUTPUT)
print(OUTPUT)
