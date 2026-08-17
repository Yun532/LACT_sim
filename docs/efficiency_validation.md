# Efficiency 曲线验证方法

本项目里的效率项都是乘法权重，不改变光线传播路径。对一个光子或 photon bunch，
最终进入相机图像的权重大致是：

```text
最终 p.e. 权重 =
  photon weight
  * mirror_reflectivity(wavelength)
  * facet_reflectivity_scale
  * filter_transmission(wavelength)
  * atmosphere_transmission(wavelength)
  * light_collector_acceptance(angle, material)
  * sipm.pde(wavelength)
  * efficiency.constant_scale
```

其中 `mirror_reflectivity`、`filter_transmission`、`atmosphere_transmission`、
`sipm.pde` 都可以关闭、设成常数，或设成两列表格：

```text
wavelength_nm,efficiency
400,0.86
401,0.861
```

## 当前程序的曲线规则

- 表格支持逗号或空格分隔。
- 以 `#` 开头或 `#` 后面的内容会被忽略。
- 表头等非数字行会被忽略。
- 波长会按从小到大排序。
- 同一波长出现多次时，会先取算术平均值，再参与插值。
- 表格范围内使用线性插值。
- 表格范围外效率为 `0`。
- 某个效率项不配置或设为 `none/off/false` 时，该项等于 `1`。

重复波长平均这一条很重要：你现在的 `sipm_pde.csv` 是从曲线数据整理来的，
里面存在同一波长多个 PDE 点。平均后曲线是确定的，不会受编译器排序细节影响。

## 单项测试

单项测试的目标是回答：“每一个输入表格被程序读入后，随波长变化的效率是否符合表格？”

运行：

```bash
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_efficiency_curves.py \
  --output-dir run_logs/official_tests/efficiency_curves \
  2>&1 | tee run_logs/official_tests/efficiency_curves/run.log
```

默认读取：

```text
configs/efficiency/mirror_reflectivity_dm0113_13point_mean.csv
configs/efficiency/filter_transmission.csv
configs/efficiency/sipm_pde.csv
```

DM0113 反射率曲线取自镜面 13 个测点的逐波长算术平均。实测覆盖
252--750 nm；表内缺失的整数波长采用相邻实测点线性插值，250--251 nm
按 252--253 nm 的局部斜率线性外推，从而保持模拟默认的 250--750 nm
有效范围。该曲线只描述平均反射率；统一或逐镜确定性乘数分别由
`efficiency.mirror_reflectivity_scale` 和
`efficiency.mirror_reflectivity_scale_csv` 控制，镜面间随机缩放误差仍由 error cfg
中的 `reflectivity_scale_sigma` 独立控制。具体优先级见
`docs/facet_parameter_precedence_zh.md`。

输出：

```text
run_logs/official_tests/efficiency_curves/efficiency_summary.png
run_logs/official_tests/efficiency_curves/efficiency_summary.pdf
run_logs/official_tests/efficiency_curves/mirror_reflectivity.png
run_logs/official_tests/efficiency_curves/filter_transmission.png
run_logs/official_tests/efficiency_curves/sipm_pde.png
run_logs/official_tests/efficiency_curves/atmosphere_transmission.png
run_logs/official_tests/efficiency_curves/efficiency_curve_report.txt
```

`efficiency_summary.png/pdf` 是推荐检查图：所有效率项画在同一张图上，同色虚线表示
理论/输入文件线，同色实线表示程序实际使用曲线，黑色粗线表示总效率。单项图里灰点是
原始输入表格点，橙色虚线是把输入表格按波长顺序连起来的理论/文件线，蓝线是程序实际使用的
“重复波长合并 + 线性插值”曲线。如果虚线和实线明显不一致，通常说明输入表格有重复波长或
异常点；如果报告里出现大量 `rows outside [0, 1]`，就需要回头检查输入文件。

## 总效率测试

同一个脚本还会输出总效率：

```text
run_logs/official_tests/efficiency_curves/total_efficiency.png
run_logs/official_tests/efficiency_curves/efficiency_curve_samples.csv
```

总效率按逐波长乘积计算：

```text
total =
  constant_scale
  * mirror_reflectivity
  * filter_transmission
  * sipm_pde
  * atmosphere_transmission
```

如果没有给 `--atmosphere-csv`，大气项默认是 `1`。如果某个表格的波长范围不覆盖当前
采样波长，该项在范围外为 `0`，所以总效率也会变成 `0`。
总效率图里的灰点表示所有输入表格唯一波长位置上的总效率，黑线表示当前采样网格上的总效率。

示例：加入一个大气透过率表格一起检查：

```bash
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_efficiency_curves.py \
  --atmosphere-csv configs/efficiency/atmosphere_transmission.csv \
  --output-dir run_logs/official_tests/efficiency_curves_with_atmosphere \
  2>&1 | tee run_logs/official_tests/efficiency_curves_with_atmosphere/run.log
```

## C++ 回归测试

`ctest` 里新增了 `test_efficiency_curves`，用于锁住以下行为：

- 表头和注释不会影响读取。
- 表格范围内线性插值正确。
- 表格范围外返回 `0`。
- 重复波长会被平均。
- 总效率是各单项效率的乘积。
- 所有效率项关闭时，总效率为 `1`。

运行：

```bash
cmake --build build -j4
ctest --test-dir build --output-on-failure -R efficiency
```

## 后续更严格的物理验证

当前测试验证的是“表格读取、插值、乘积”是否正确。更接近物理链路的验证可以再加两类：

1. 单色光源扫描：固定波长分别跑 300、350、400、450、500 nm，检查输出 p.e. 比例是否等于
   对应总效率比例。
2. 已知光谱加权：给一个人工光谱，比较模拟输出总 p.e. 与
   `sum(spectrum * total_efficiency)` 的解析计算。

这两类测试会更慢一些，适合作为 full-response benchmark，而不是每次快速回归都跑。
