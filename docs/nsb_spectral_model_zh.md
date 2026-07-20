# NSB 光谱模型第一版记录

本文记录已经接入 LACT_sim 的夜空背景光 NSB 第一版计算方案。目标是把
`nsb.rate_pe_per_ns_per_pixel` 从手动测试常数，改成由 LoNS 光谱、固定有效面积、
像素视场和效率曲线计算得到的物理量。

## 输入来源

参考文件：

```text
/Users/yun/Library/Mobile Documents/com~apple~CloudDocs/LACT/nsb/LoNS.ipynb
```

notebook 中使用 ESO SkyCalc 生成 LoNS 光谱。SkyCalc 原始 flux 单位为：

```text
ph / (s um arcsec^2 m^2)
```

程序中需要转换为：

```text
ph / (s nm sr m^2)
```

转换公式：

```text
LoNS_flux_ph_s_nm_sr_m2 =
    LoNS_flux_skycalc * (206265)^2 / 1000
```

其中：

- `(206265)^2`：`1 sr` 中的 `arcsec^2` 数。
- `/1000`：把每 `um` 转成每 `nm`。

注意：当前目录里的 `lons_eso_sky.txt` 是 404 HTML，不是有效光谱文件。第一版应以
notebook 的 SkyCalc 生成逻辑或导出的有效 LoNS 光谱表为准。

## 第一版固定有效面积

先不在每次 CORSIKA 模拟里重新做平行光标定。固定使用 official 纯光学平行光
得到的白板有效收集面积：

```text
不考虑遮挡:
  A_eff_no_obstruction = 29.623570 m^2

考虑遮挡:
  A_eff_with_obstruction = 24.576860 m^2
```

来源：

```text
run_logs/official_tests/perfect_parallel/run.log
run_logs/official_tests/raytrace_structure_parallel/run.log
```

标定运行沿用上述 official 配置，仅把平行光光子数提高到 `10,000,000` 以降低蒙特卡洛统计误差。

计算方式：

```text
A_eff = source_sampling_area_m2 * hit_output_plane / total_photons
```

official 平行光测试中：

```text
source_sampling_area_m2 = 50.265482 m^2
total_photons           = 10,000,000
```

统计结果：

```text
无遮挡:
  hit_output_plane = 5893422
  A_eff = 29.623570 m^2

考虑遮挡:
  hit_output_plane = 4889411
  A_eff = 24.576860 m^2
```

遮挡透过比例：

```text
24.576860 / 29.623570 = 0.829639
```

遮挡损失比例：

```text
17.0361 %
```

这两个面积暂时只代表纯光学白板有效面积，不再额外包含镜面反射率、滤光片、
SiPM PDE、大气等波长相关效率。这些效率在 NSB 光谱积分时单独相乘。

## 像素视场

notebook 中示例使用：

```text
fov_x_deg = 0.18
fov_y_deg = 0.18
Omega = deg2rad(fov_x_deg) * deg2rad(fov_y_deg)
      = 9.8696044e-6 sr
```

LACT `new_camera` 当前像素边长为 `2.44 cm`，白板/焦平面距离为 `8 m`，可用小角近似：

```text
Omega_pixel = (pixel_size_m / focal_length_m)^2
            = (0.0244 / 8)^2
            = 9.3025e-6 sr
```

第一版推荐在程序里默认用相机配置自动计算 `Omega_pixel`。为了和 notebook 结果对照，
独立测试图中可以同时标出 notebook 使用的 `0.18 deg x 0.18 deg` 视场。

## 计算公式

LoNS 光谱积分得到每像素 p.e. rate：

```text
rate_pe_per_ns_per_pixel =
    1e-9 * A_eff * Omega_pixel *
    ∫ LoNS(lambda)
      * mirror_reflectivity(lambda)
      * filter_transmission(lambda)
      * sipm_pde(lambda)
      * atmosphere_factor(lambda)
      d_lambda
```

其中：

- `LoNS(lambda)` 单位为 `ph / (s nm sr m^2)`。
- `A_eff` 单位为 `m^2`。
- `Omega_pixel` 单位为 `sr`。
- 各效率曲线为无量纲。
- 积分结果先是 `pe / s / pixel`，乘 `1e-9` 后为 `pe / ns / pixel`。

notebook 备注中提到：如果 LoNS 来自 SkyCalc 的 sky radiance，通常不再额外乘 SkyCalc
输出的 `trans`，因为 sky radiance 已经是到达观测点的天空辐射。LACT_sim 中是否再乘
额外大气项，应由 `nsb.atmosphere_factor` 显式配置控制。

## 从 notebook 提取的参考数值

`background_summary.txt`，无月示例：

```text
Omega_sr                      = 9.8696044011e-06
Area_m2                       = 25.0
dt_s                          = 30e-9
Integral of black curve       = 1.0184330663e13 ph s^-1 sr^-1 m^-2
Integral of red curve         = 3.5366885168e11 pe s^-1 sr^-1 m^-2
Detected p.e. in 30 ns        = 2.6179287413
```

`background_summary_with_moon.txt`，有月示例：

```text
Omega_sr                      = 9.8696044011e-06
Area_m2                       = 25.0
dt_s                          = 30e-9
Integral of black curve       = 1.6261508662e13 ph s^-1 sr^-1 m^-2
Integral of red curve         = 1.2325963209e12 pe s^-1 sr^-1 m^-2
Detected p.e. in 30 ns        = 9.1239285550
```

把 notebook 的红线积分换到 LACT `new_camera` 的像素视场 `9.3025e-6 sr`，并使用上面的
两个固定面积，可得到第一版 NSB rate 量级：

```text
无月，无遮挡:
  rate = 0.097462 pe/ns/pixel
  mean = 2.437 pe / 25 ns / pixel
  mean = 2.924 pe / 30 ns / pixel

无月，考虑遮挡:
  rate = 0.080858 pe/ns/pixel
  mean = 2.021 pe / 25 ns / pixel
  mean = 2.426 pe / 30 ns / pixel

有月示例，无遮挡:
  rate = 0.339671 pe/ns/pixel
  mean = 8.492 pe / 25 ns / pixel
  mean = 10.190 pe / 30 ns / pixel

有月示例，考虑遮挡:
  rate = 0.281804 pe/ns/pixel
  mean = 7.045 pe / 25 ns / pixel
  mean = 8.454 pe / 30 ns / pixel
```

这些数值比当前 smoke test 的：

```ini
nsb.rate_pe_per_ns_per_pixel=0.05
```

更接近物理输入，但仍依赖 SkyCalc 条件、视场定义、固定面积选择和效率曲线。

## 配置接口

第一版保留当前 constant-rate 路径，同时新增 spectral-rate 计算路径：

```ini
nsb.enabled=true
nsb.model=spectral_flux
nsb.spectrum_csv=configs/nsb/lons_skycalc_dark.csv
nsb.spectrum_unit=ph_s_nm_sr_m2

# 二选一：无遮挡或考虑遮挡。
nsb.effective_area_m2=29.623570
# nsb.effective_area_m2=24.576860

nsb.pixel_solid_angle=auto
nsb.seed=12345
```

主程序启动时先把光谱积分成：

```text
computed_rate_pe_per_ns_per_pixel
```

后续 NSB 采样仍沿用当前逻辑：

```text
waveform 开启:
  NSB_pe(time_bin, pixel) ~ Poisson(rate * time_bin_width_ns)

waveform 关闭:
  NSB_pe(pixel) ~ Poisson(rate * nsb.window_ns)
```

独立检查命令：

```bash
build/compute_nsb_rate configs/nsb/spectral_rate_check_with_obstruction.cfg
```

当前 official NSB+trigger 测试已经使用 spectral-rate 配置：

```text
configs/official_tests/corsika_nsb_trigger_camera.cfg
configs/official_tests/corsika_obstruction_nsb_trigger_camera.cfg
```

## 独立 NSB 测试建议

独立测试程序和画图脚本：

```bash
build/compute_nsb_rate configs/nsb/spectral_rate_check_with_obstruction.cfg

MPLBACKEND=Agg python3 python/plot_nsb_spectral_rate.py \
  --effective-area-m2 24.576860 \
  --output run_logs/manual_checks/nsb_spectrum_from_notebook/nsb_spectral_response_with_poisson.png \
  --diagnostic-csv run_logs/manual_checks/nsb_spectrum_from_notebook/diagnostic_from_script.csv \
  --summary run_logs/manual_checks/nsb_spectrum_from_notebook/summary_from_script.txt
```

建议输出：

```text
run_logs/official_tests/nsb_spectral/diagnostic.csv
run_logs/official_tests/nsb_spectral/summary.txt
run_logs/official_tests/nsb_spectral/nsb_spectral_response.png
```

`diagnostic.csv` 建议包含：

```csv
wavelength_nm,
lons_flux_ph_s_nm_sr_m2,
mirror_reflectivity,
filter_transmission,
sipm_pde,
atmosphere_factor,
total_efficiency,
contribution_pe_ns_pixel_nm
```

图像建议画在一张论文风格图中：

1. 左轴 log：LoNS 光谱和 `LoNS × total_efficiency`。
2. 右轴 linear：mirror/filter/PDE/total efficiency。
3. 图内标注：
   - `A_eff`
   - `Omega_pixel`
   - `rate_pe_per_ns_per_pixel`
   - `mean_pe_per_pixel` for 25 ns 或 30 ns。

理论与实际模拟对比：

1. 用计算出的 `rate_pe_per_ns_per_pixel` 作为 NSB 输入。
2. 生成大量 pixel/time-bin Poisson 样本。
3. 对比：
   - 理论均值 `rate * window_ns`。
   - 模拟均值。
   - 理论 RMS `sqrt(rate * window_ns)`。
   - 模拟 RMS。
4. 画 histogram：模拟 NSB p.e. 分布与 Poisson 理论曲线叠加。

这样可以独立验证 NSB rate 计算和 Poisson 注入是否一致，再接入正式 CORSIKA
`run_corsika_trace`。
