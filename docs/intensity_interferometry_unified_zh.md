# LACT 恒星强度干涉统一模拟

## 结论先行

本实现以 `main` 的光学、S17351 微单元、NSB 与实测 SPE 标定为底座，合并：

- 旧 v1 的多源复可见度、全阵列所有基线、地球自转覆盖和无相位重建；
- 旧 v2 的飞秒光学相干、纳秒电子波形、ADC/噪声/延迟与小时级相关器分层。

统一分支不把一整夜生成为 625 MS/s 波形。短波形用于校准电子学偏差；长曝光使用
带显式误差的 `|V|²` 充分统计量。这两层由同一源模型、通带和仪器参数连接。

旧分支最终提交为：

- `intensity-interferometry`: `a6e2817b42669be4c8d6f46e34d99e1495098e2e`
- `intensity-interferometry-v2`: `d64f78681e30facd4e33f1defb336ce2c767dd6a`

## 一键生成和执行 notebook

在仓库根目录运行：

```bash
python tools/build_sii_paper_notebook.py
python tools/execute_notebook.py notebooks/lact_sii_paper_simulation.ipynb --timeout 1200
```

notebook 包括：

1. 参数来源与“已有/假设/缺失”表；
2. `layout_0803` 的 32 台望远镜实际坐标和 496 条基线；
3. 理论双星天图和完整理论 `(u,v)` 功率图；
4. 任意两镜的 ENU 基线、正确 topocentric `(u,v,w)`、`w/c` 几何时延；
5. 全阵列 6 小时地球自转覆盖；
6. 590239条轴上光线导出的镜面到达时间核，以及170 ns实测 SPE、恢复、电子噪声和 ADC 波形诊断；
7. 逐望远镜 Poisson 光子/NSB计数、共享增益/时钟/透明度误差和带协方差来源的长曝光 `|V|²` 模拟；
8. 正值、有限支撑、平滑正则、多起点的无相位重建；
9. 星等极限和角尺度情景；
10. 2/6/10 h 单夜、6/30/60 h 多夜、0/1/2 倍 NSB、理想/现实电子学的逐情景重建；
11. 100次望远镜级参数 Monte Carlo；
12. 16个“星等×积分时间”组合的960张盲图：`m_AB=2–6`、2/6/30/60小时、双星注入与单星空白对照，以及精确95%检出区间；
13. 自动物理闭合检查。

输出表和图片位于 `run_logs/sii_paper_notebook/`。

## 更暗双星的盲检出边界

这里的“极限”专指 notebook 默认的 0.20 mas、次/主星流量比 0.55 双星，不可直接外推到
更小分离、更极端流量比或一般表面结构。所有星等使用同一重建网格、支撑与正则参数。
每个“星等×积分时间”组合包含30次双星注入和30次单星空白；只有真阳性率的95%置信
下限不低于80%，同时假阳性率的95%置信上限不高于20%，才标记为可靠。

| 总积分时间 | 本次离散网格中最暗的可靠测试点 | 下一档失败原因 |
|---:|---:|---|
| 2 h | `m_AB=4` | 5等：27/30双星检出，但16/30单星误报 |
| 6 h | `m_AB=4` | 5等：29/30双星检出，但4/30单星误报，置信上限未过门槛 |
| 30 h（5夜） | `m_AB=5` | 6等：27/30双星检出，但21/30单星误报 |
| 60 h（10夜） | `m_AB=5` | 6等：27/30双星检出，但15/30单星误报 |

这是已测试的整数星等网格，不表示连续极限精确等于4.00或5.00等。参考暗天 NSB 为
70.527 MHz/pixel；2–6等的已探测星光率依次约为213.45、84.98、33.83、13.47和
5.36 MHz/pixel。越过边界后，问题不是优化器完全画不出双峰，而是单星噪声也频繁产生
满足阈值的伪双峰。

## 两台望远镜实际测量什么

望远镜不直接测复电场，也不直接得到傅里叶相位。每台镜得到随时间变化的光电流或
ADC 电压 `I_i(t)`。先用阵列坐标和源方向计算几何时延

```text
tau_g = dot(B_ij, source_direction) / c
```

校正两路时间轴，再累计归一化互相关：

```text
g2_ij(tau) - 1 = <delta I_i(t) delta I_j(t+tau)> / (<I_i><I_j>)
```

校正后零时延峰的面积与 `|V(u,v)|²` 成正比。该时刻的 `(u,v)` 是基线在天空东、
天空北方向的投影并除以波长。地球转动后，同一对镜子会写出一段 `(u,v)` 轨迹。

## S17351 恢复模型

厂家公开网站没有检索到 S17351 的恢复时间数据表。Hamamatsu 的 MPPC 技术说明指出，
恢复常数由淬灭电阻与结电容的 RC 决定，并给出典型像元约 15 ns 的量级。因此：

- 默认 `recovery_time_ns=10`，明确标记为暂定值；
- notebook 固定扫描 1、10、30 ns；
- C++ 链路使用 `1-exp(-delta_t/tau)` 作为重复 avalanche 的电荷恢复比例；
- 获得实测数据后只替换配置，不改变算法。

可用配置为：

```ini
microcell.saturation_enabled=true
microcell.model=explicit_exponential_recovery
microcell.recovery_enabled=true
microcell.recovery_time_ns=10.0
```

示例文件：`configs/electronics/measured_spe_recovery_10ns_1p6ns.cfg`。

在正常星光和 NSB 下，光电子随机分配到 270336 个微单元，平均单元占用
`rate * tau / N_cell` 很小。恢复时间通常不会主导强度干涉灵敏度；ADC 带宽、NSB、
通光效率、时延稳定性和系统标定更重要。强激光、局部聚焦或高亮脉冲则必须逐单元模拟。

## 仍不能宣称为实测性能的部分

当前电子带宽、ADC 规格和 0.35 mV 噪声沿用工程假设；短波形逐 PE 抽取 537 个
实测电荷因子，长曝光用 SPE shot-noise 方差近似传播电子噪声和 ADC 量化损失。正式
1229分片镜面轴上理想误差配置的590239条400 nm光线给出约0.602 ns RMS、1.732 ns
峰峰值到达时间展宽，已经进入短波形和长曝光传递效率；原始128 MB逐光子表压缩为
54个逐镜面混合分量并保存来源哈希。
2 nm 窄带滤光片的角度响应、
实际过压/温度下暗计数、串扰、后脉冲、随仰角/离轴角变化的光学时间核、实测跨镜时钟漂移与零基线标定
仍缺实测输入。因此 notebook 给出的是可复现的研究预测和参数敏感性，不是 LACT 已验收
的极限星等。

Hamamatsu 参考：
[Physics and operation of the MPPC silicon photomultiplier](https://hub.hamamatsu.com/us/en/technical-notes/mppc-sipms/physics-and-operation-of-the-MPPC-silicon-photomultiplier.html)。
