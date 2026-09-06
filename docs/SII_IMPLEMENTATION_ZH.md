# SII 代码如何把一个天体算成 LACT 观测与推断

本文按数据流阅读程序：先给定天体与仪器，把光变成探测事件和电压，处理两镜时差，从整段相关曲线估计平方可见度，再构造整夜阵列观测，最后用同一测量模型反推天体。每一节说明实际调用、输入、计算和输出；较长的物理推导链接到[物理说明](SII_PHYSICS_ZH.md)。

主入口是[完整 notebook](../notebooks/sii_complete_waveform_report.ipynb)，其源码由[生成脚本](../tools/build_sii_science_notebook.py)维护。Notebook先看理论和中间图，再做模型与性能比较；本文则把分散在不同单元和模块中的同一条计算链连起来。复现命令在末节。

需要先分清两种产物。微秒记录真的逐事件生成、逐采样求电压、逐滞后算相关。六小时数据使用短波形标定的响应与误差，抽样长曝光的充分统计量，并没有生成六小时连续 ADC。两者通过同一 SPE、光学时间核、光子率、采样相位和 GLS 定义连接。以下所有“观测”均指这套已注明条件的模拟，实测输入不等于实测观测结果。

## 1. 程序启动时，先固定一套不会中途改变的场景

Notebook设置随机种子 `SEED=20260905`，记录 Python、NumPy、SciPy 版本，将 BLAS 线程设为1，然后调用 `verify_main_parameters(ROOT)` 和 `Instrument.from_repository(ROOT)`。前者核对输入内容，后者真正解析配置。它们读取本地检出的文件，不在每次运行时联网获取 main。

基准计算的输入如下。表中数值是本例，不是所有库函数的内建默认值；例如直接调用 `Instrument()` 不会自动获得 main 的4 ns采样和实测模板。

| 量 | 本例取值 | 在程序中的位置 |
|---|---|---|
| 天体 | AB 2等，直径0.16 mas均匀盘 | `BinarySource(primary_diameter_mas=.16)`，`source_case='single_disk'` |
| 仪器 | 400 nm中心，2 nm全宽，双偏振 | `Instrument.from_repository`，`polarization_factor=.5` |
| 阵列 | 用户提供的36台input，630条镜对基线 | `configs/arrays/lact36_20260906.input` |
| 观测 | 赤纬22°，纬度29.36°，过中天前后合计6小时 | `Observation`，曝光单位为SI秒 |
| 输出分段 | 每段1200 s，一夜18段，共11340条测量 | `segment_s`、`nights` |
| 可见度积分 | 每段9个时间中点，基准通带5个波长点 | 每条测量45个子采样，均不增加独立观测数 |
| 电压采样 | 250 MS/s，间隔4 ns | 实测SPE在采样中心取值 |
| 教学短记录 | Tel.1与Tel.2，时角0.5 rad，24 μs，注入倍率1 | `sii_walkthrough`；这一时刻不是0.5小时 |
| 图像计算 | 32×32像素，0.70 mas视场，0.32 mas支撑半径 | `reconstruct_uv_data`的显式参数 |

这里把亮度形态与总通量分开。改变圆盘直径只改变归一化可见度；保持AB星等不变就保持总星光率。双星、椭圆和遮挡盘在后文作为明确更换的源模型，不会偷偷改变基准单镜电压的参数。

输入核对清单是[main_parameter_manifest.json](../configs/sii/main_parameter_manifest.json)，对应 main 提交 `2926031b14c8aa2164a4d9233f5c2d0a75324127`。清单以统一LF换行后的SHA-256检查10个关键文件。额外光追缓存、布局和数值代码也写入本次结果摘要。缓存的原始C++光追并未在本次notebook中重跑，来源与复现边界见物理说明第1、3节。

## 2. 先把配置变成三个真正控制信号的量

接下来的事件生成首先需要星光率、背景率和相干面积，而不是一个笼统的“探测效率”。`Instrument.from_repository`从 `corsika_lact_pylast_root_only_measured_waveform.cfg` 展开组件，读镜面反射率、滤光片透过率、PDE、天空谱、像素尺寸与焦距，再读电子学文件中的SPE、电荷样本和微单元结构。

记波长为 $\lambda$，三条效率曲线乘积再乘可选SII滤光片得到 $R(\lambda)$。$R$ 无量纲，曲线之间用线性插值。恒星采用平坦 $F_\nu$ 的AB谱，因而真正积分的是

$$
J_1=\int \frac{R(\lambda)}{\lambda}\,d\lambda,
\qquad J_2=\int R^2(\lambda)\,d\lambda.
$$

同一式内波长单位必须一致。$J_1$ 无量纲，$J_2$ 的单位随波长网格为nm；程序中的比例式消去该单位。`throughput` 中的通带校正为 $J_1/[R(\lambda_0)\Delta\lambda/\lambda_0]$，这里 $\lambda_0$ 是中心，$\Delta\lambda$ 是全宽。它使外层保留中心谱密度接口，同时得到通带积分后的星光率。

中心有效探测面积优先取400 nm光追缓存 `lact2_measured_single_pixel_400nm.provenance.json` 中的5.92276165 m²。这一面积已经含镜面、遮挡、PSF、像素接收、集光器、滤光片与PDE，代码用它除以面积字段24.57686 m²得到中心效率，再乘源透射假设0.7836336和上述通带校正。不能再乘一遍中心PDE。完整核缺失时库有旧光学核回退，所以正式结果还核对实际输入路径和哈希。

`detected_star_rate_hz`使用AB零点和光子能量把结果转成Hz。基准得到 $R_\star=1.99889298\times10^8\ \mathrm{s^{-1}}$。NSB另算 $\int L_\lambda R\,d\lambda$，乘面积、集光器和像素立体角近似 $(\mathrm{pixel\ size}/\mathrm{focal\ length})^2$，得到0.573214 MHz；不再乘源的上游透射系数。天空谱是既有SkyCalc情景，不能称为LACT当晚实测。启用暗计数时仅加一次 $R_b=R_{\rm NSB}+R_{\rm dark}$。

相关信号使用的是相干面积 $\tau_{\rm eff}$，不是ADC宽度。代码设置

$$
\tau_{\rm eff}=\frac{p}{\Delta\nu}\,
\frac{\Delta\lambda J_2}{\lambda_0^2 J_1^2},
\qquad \Delta\nu=\frac{c\Delta\lambda}{\lambda_0^2}.
$$

$p=.5$ 是偏振系数，$c$ 为光速，$\tau_{\rm eff}$ 单位为s。它来自探测光谱平方的积分，推导见[物理说明附录C](SII_PHYSICS_ZH.md#appendix-c)。本例为 $1.334291922\times10^{-13}$ s，已经包含偏振，后面不能再乘一次0.5。

这一步同时生成 `visibility_wavelength_nm` 和 `visibility_spectral_weights`。后者按 $R^2$ 加权并归一化，而不是按光子率的一次方加权。默认平滑2 nm通带使用密集梯形积分计算率和面积、5点Gauss–Legendre计算可见度。若给 `sii_bandpass_path`，程序改在所有输入曲线的相邻物理节点之间做5点Gauss–Legendre积分，并用同一组节点积分率、面积和可见度权重。因此399.2 nm附近很窄的透射峰也不会落在全部节点之间。全零或非有限响应明确报错。窄通带改变积分点数，但不代表增加探测通道。

这一节的输出仍只有仪器常数和求积网格，没有随机数、没有天体图像，也没有观测误差。

## 3. 用同一套坐标决定“在哪里看”和“何时到达”

`read_corsika_layout`读取TELESCOPE卡片，保留原始坐标、标签、行号和半径。原始单位为cm，三个方向为北、西、上；转成ENU米遵循原项目的转换：

$$
(E,N,U)=(-y_{\rm input},x_{\rm input},z_{\rm input})/100.
$$

程序没有重新居中或额外旋转。Tel.1为 $(-434.92063,417.46032,0)$ m，Tel.36为 $(266.66667,136.50794,0)$ m。400 cm是4 m半径，不能当作镜直径。派生坐标见[CSV](../configs/arrays/lact36_20260906_coordinates.csv)；检查脚本会从input重新转换后比较。

`generate_uvw`遍历所有 $i<j$，定义基线 $\boldsymbol B_{ij}=\boldsymbol r_j-\boldsymbol r_i$。每段取时间中心，调用 `celestial_tangent_axes_enu` 得天球的三个单位轴，再由 `uvw_from_enu`点乘得到米单位的 $(u_m,v_m,w_m)$。前两轴固定于赤道坐标切平面，第三轴指向源，不用当地天顶与源的叉乘引入额外视场转动。

时角定义为 $H=\mathrm{LST}-\mathrm{RA}$，西向为正。SI秒转时角用平均恒星日86164.0905 s，因此6小时窗口的角宽并非直接把SI小时当恒星时小时。程序要求总时长是完整段数，并拒绝源在该观测窗口的检查点落到地平线以下。

几何表保留 `u_m/v_m/w_m`、基线ENU、镜号、段号和时角。用中心波长除前两项得 `u_lambda/v_lambda`，单位是“波长数”。表内 `geometric_delay_ns=w_m/c` 是几何投影约定；事件到达的右镜减左镜时间差为它的负值。`geometric_arrival_delays_ns`与`tracking_geometry`使用后一约定：

$$
d_i=-\frac{(\boldsymbol r_i-\boldsymbol r_{\rm ref})\cdot\boldsymbol s(H)}{c}.
$$

$\boldsymbol s$ 是源方向单位向量，$d_i$ 是相对参考镜的到达时差，换成ns交给波形。几何投影符号与到达时差符号必须区分；否则补偿方向会反。源方向与切向轴的具体分量见物理说明第4节。

`tracking_geometry`还对时角求导，输出每秒变化多少ns的时延率及二阶曲率界。短记录用 $d_i(t)=d_i(0)+\dot d_i t$，局部时间从零起算，避免把整夜绝对纳秒数放进插值丢精度。其近似误差用曲率界乘记录时长平方的一半检查。

## 4. 在这些位置算出理想可见度，然后才做曝光平均

`source_visibility`返回复振幅 $V$；外层取模平方得到 $P=|V|^2$。可见度是归一化亮度的傅里叶变换，零基线为1，公式及推导见[物理说明附录A](SII_PHYSICS_ZH.md#appendix-a)。各分支的真实运算如下。

均匀盘以 $q=\sqrt{u^2+v^2}$ 及弧度直径 $d$ 计算 $V=2J_1(\pi dq)/(\pi dq)$。`uniform_disk_visibility`单独处理零自变量为1。mas先乘 `MAS_TO_RAD`；贝塞尔函数变负时保留振幅符号，到外层才平方。

双星先分别算两个圆盘的 $V_1,V_2$，再计算 $(V_1+rV_2e^{-2\pi i(u\Delta E+v\Delta N)})/(1+r)$。$r$ 为次星/主星流量比，偏移 $(\Delta E,\Delta N)$ 由分离和从北向东的位置角决定，单位为rad。因此有干涉交叉项，不能平均两个盘的功率。椭圆先旋转UV坐标，再按长、短轴缩放自变量，使用同一圆盘贝塞尔形式。

静态遮挡盘用恒星圆盘振幅减去带位移相位的行星圆盘振幅，按剩余通量归一化。挡光比例取面积比；程序要求行星完整位于恒星盘内，不把部分入凌几何交集当成已实现。它是单个时刻的遮挡形态，不是随轨道演化的凌星光变模型。

`segment_uv_samples`只使用几何与仪器，不读源。它在每段9个时间中点重新投影，并按各波长重新除以波长。`segment_sampling_weights`给出时间等权乘光谱权重。第 $a$ 条测量的期望是

$$
\overline P_a=\sum_{q=1}^{n_a}\omega_{aq}|V(u_{aq},v_{aq})|^2,
\qquad \sum_q\omega_{aq}=1.
$$

$q$ 标识该条测量的时间/波长节点，$n_a=45$ 是基准节点数。必须逐节点平方后平均；先平均复振幅或只在代表UV坐标计算均不等价。`simulate_uv_observation`保存中心真值、平均真值以及 `uv_samples_u/v/weight`。中心真值用于量化曝光模糊；反演使用的是完整节点，真值列不传给反演器。

到这里可以画连续理论平面，也可以画11340个有限曝光理论点。两者都没有加入噪声；图上的像素密度不等于观测数。

## 5. 产生光电子：二阶信号如何进入连续事件时刻

理论 $P$ 进入短记录时，`simulate_hbt_primary_pe`处理两镜，`simulate_array_photon_times`处理共享的多镜事件流。两者都用弱光下的稀疏共享对近似，而不逐飞秒生成电磁场。

两镜共享对率为

$$
R_{ij}^{\rm pair}=R_{\star,i}R_{\star,j}\tau_{\rm eff}\overline P_{ij}.
$$

其单位是s⁻¹。生成器先抽泊松个数，再在连续记录区间均匀抽参考时刻，把同一个参考时刻加入两个镜的事件流。各镜还生成独立星光和独立背景。为了保持给定单镜星光率，多镜情况下必须扣除该镜参与的全部共享对率：$R_i^{\rm single}=R_{\star,i}-\alpha\sum_{j\ne i}R_{ij}^{\rm pair}$。$\alpha$ 是已知注入倍率；若独立率变负，立即拒绝该倍率。

正常物理记录取 $\alpha=1$。增大倍率只用于短样本响应标定，输出同时记录物理率和注入率，并在估计响应时除以倍率。它不能使记录变成等效真实长曝光，也不能把放大注入的高阶协方差当作天体的协方差。

多镜入口先从 `source_coherence_spectrum`取得每个波长的复相干矩阵 $\Gamma^{(r)}$，检查Hermitian、单位对角和半正定，再算 $\sum_r w_r|\Gamma_{ij}^{(r)}|^2$。每台望远镜只有一个事件流，参与它的所有基线都使用这同一个流，不能为每一镜对重新渲染一台虚拟望远镜。

星光事件按 $t\mapsto t+d_i+\dot d_i t$ 到达探测端，再抽独立光学延迟。时间伸缩的Jacobian会使接收星光率成为 $R_{\star,i}/(1+\dot d_i10^{-9})$，元数据分别记录参考率和接收率；背景按接收时间独立生成。生成区间两端扩展到包含几何时差、SPE支持区及光学核保护区，避免记录开头无故少光子。

光学延迟来自 `load_optical_timing_mixture`读出的混合分布：按权重选分量，再按该分量均值、标准差抽高斯延迟。相关对的参考时刻共享，两端光学延迟独立。因此保留超额对数，却把飞秒相关面积展宽到纳秒响应里。

基准退化参数 $R_\star\tau_{\rm eff}\approx2.67\times10^{-5}$，弱对模型适合验证这一条件下的二阶响应；它没有生成完整热光自聚束、三阶闭合相位和全部四阶累积量。独立复高斯热光模检验及弱光误差预算放在第13节，不把缺失项藏在代码标签里。兼容main事件CSV时 `origin='cherenkov'` 只是原接口信号标签，此处源仍为恒星。

## 6. 同一批光电子怎样变成单镜电压

`render_pe_waveform`接收连续ns时刻，先加入已启用的本征高斯抖动，再为每个事件抽电荷因子。`load_empirical_charge_factors`把实测电荷样本归一到均值1，因此保留幅度涨落而不重复改变平均SPE增益。二阶矩为 $F_q^2=\langle q^2\rangle=1.03254526$，`excess_noise_factor`保存其平方根。

若微单元恢复时间非零，为事件抽微单元号，逐单元按事件时序更新恢复状态；恢复因子为 $1-e^{-\Delta t/\tau_r}$。$\Delta t$ 是该单元前后事件间隔，$\tau_r$ 是独立给定的恢复时间。当前 $\tau_r=0$ 时全部因子为1，也不抽无意义的微单元随机号。SPE电压尾部不能直接解释为这一恢复时间。

随后在采样中心 $t_n=(n+1/2)\Delta t_{\rm ADC}$ 求

$$
x_i[n]=\sum_k q_{ik}f_{ik}h(t_n-t_{ik})+\epsilon_i[n].
$$

$h$ 为实测SPE模板，单位mV；$t_{ik}$ 是第 $k$ 个事件到达时刻，$q_{ik}$ 是电荷因子，$f_{ik}$ 是恢复比例，$\epsilon$ 是可选白噪声。代码只计算每个SPE支持范围内的采样点，分块构造索引，用线性插值读连续平移模板，再用 `bincount` 累加。没有先把光子时间取整到4 ns格；同一脉冲的不同采样相位确实能产生不同数组。

噪声RMS为0时直接复制模拟电压，不额外抽噪声。ADC位数为0时保留浮点电压；启用时 `digitize_adc`按满量程剪裁和量化。输出包含 `sample_time_ns`、`analog_mv`、`noisy_mv`、`adc_mv`、电荷与恢复因子。24 μs的每镜6000点是真正计算的ADC数组。

教学图调用 `single_record`一次，把事件和ADC放入同一个 `lesson_scene`；随后 `pair_correlation`直接读取这一对象。保存的 [walkthrough_record.npz](../validation/sii_science/walkthrough_record.npz) 包含两镜事件、ADC、采样中心、时差和种子。检查脚本既按种子重放，也从保存ADC重算相关，防止前后图暗中换了一次随机记录。

## 7. 先补偿几何时延，再计算相关曲线

主相位处理使用 `integer_align`。对每镜到达时差 $d_i$ 取最接近的整数样本 $n_i=\mathrm{round}(d_i/\Delta t_{\rm ADC})$，读 $x_i[n+n_i]$，并裁剪所有镜共同有真实输入的区间。没有周期卷绕，也不对ADC做分数插值。两镜剩余时差为

$$
\delta_{ij}=(d_j-n_j\Delta t_{\rm ADC})-(d_i-n_i\Delta t_{\rm ADC}).
$$

它可落在 $[-\Delta t_{\rm ADC},\Delta t_{\rm ADC}]$，会改变采样相关峰，所以必须交给下一节的连续模板。教学记录右镜平移58点即232 ns，剩余 $-0.967249$ ns；6000点裁成5942点，实际曝光23.768 μs。估计误差用裁后的曝光。

`waveform_cross_correlation`对两条数组分别减块均值，使用FFT求完整线性互相关，以每个滞后的实际重叠点数 $N-|\ell|$ 除回，并以两条波形的样本标准差归一：

$$
\widehat C[\ell]=\frac{\sum_{n\in\mathrm{overlap}}(x[n+\ell]-\bar x)(y[n]-\bar y)}
{(N-|\ell|)s_xs_y}.
$$

$\ell$ 是整数滞后，乘采样宽度得到ns；$s_x,s_y$单位mV，所以 $\widehat C$ 无量纲。这是电压涨落的相关系数，不是光学 $g^{(2)}$；分母含散粒噪声、SPE电荷涨落和已启用的电噪声。默认只保留±200 ns滞后，用整段相关曲线进行拟合，不能把某个随机最大值当作HBT峰。

`correlate_blocks`可对阵列使用同一组完整块，返回“块×基线×滞后”数组和等曝光平均，不足一个块的尾部单独记录后丢弃。单个常量通道不能做有意义的归一化，阵列入口会拒绝。

另一个 `align_waveforms` 入口是分数Lanczos–sinc插值及块内时延率追踪，用于处理对照。它按ADC半样本中心计算实际读取坐标，显式裁掉无完整抽头支持的边缘；半宽默认16，并做16至512的对照。插值会同时改变响应和噪声，必须标定插值后的处理。它没有与整数移位的标定静默混用。

## 8. 从相关曲线估计一个平方可见度，而非读取零滞后

`WaveformGLSCalibration`的核心是：滞后坐标、零信号均值 $\boldsymbol C_0$、单位可见度模板 $\boldsymbol k$、标定块协方差 $\boldsymbol\Sigma_b$、块时长 $T_b$，以及共享响应增益误差。物理模型为 $E[\widehat{\boldsymbol C}]=\boldsymbol C_0+P\boldsymbol k$。

相邻滞后由同一波形构成，噪声不独立。`waveform_gls_weights`求线性方程得到

$$
\boldsymbol w=\frac{\boldsymbol\Sigma_b^{-1}\boldsymbol k}
{\boldsymbol k^T\boldsymbol\Sigma_b^{-1}\boldsymbol k},
\quad \widehat P=\boldsymbol w^T(\widehat{\boldsymbol C}-\boldsymbol C_0),
\quad \sigma_P(T)=\sqrt{\frac{T_b/T}{\boldsymbol k^T\boldsymbol\Sigma_b^{-1}\boldsymbol k}}.
$$

$T$为实际曝光秒，转置记为上标 $T$；逆矩阵在代码中通过线性求解实现。$\boldsymbol w^T\boldsymbol k=1$保证模板幅度正确恢复，最小方差推导见[物理说明附录H](SII_PHYSICS_ZH.md#appendix-h)。负权重可用于抑制共同噪声，负 $\widehat P$ 则是有限噪声估计，均不应截断。

### 8.1 连续SPE决定任意分数时延下的峰形

`analytic_waveform_calibration`先在0.05 ns细网格插值SPE，计算自相关 $A_h(\tau)=\int h(t)h(t+\tau)dt$。星光共享对的相关响应是它与两端光学延迟差分布的卷积，频域通过乘光学传递函数的模平方完成。本征时间抖动也乘入同一传递函数。单镜平稳泊松自相关保留SPE散粒噪声，光学到达抖动不额外降低单镜泊松噪声。

这个连续相关核在 $\ell\Delta t_{\rm ADC}+\delta$ 处采样，再按光子率与波形方差归一，并保留块内减均值的一阶有限时长修正，得到 $\boldsymbol k(\delta)$。零信号协方差由单镜归一自相关的Bartlett卷积得到。该解析模型支持线性SPE和独立白噪声；非线性ADC、恢复、串扰和后脉冲不在此闭式模型内。

主标定现在明确使用这一连续峰形，再让波形Monte Carlo校准增益与协方差。它不靠平移已经粗采样的零相位峰去猜2 ns处的形状。`phase_template_bank`默认在−4到4 ns放65个节点，为每个节点求模板、权重和块误差。输入一个已有MC标定时，要求其 `phase_model='linear_spe'`，且零相位模板确实是连续模型的单一增益缩放，块时长也必须相同。只有任意离散峰、没有连续相位模型的旧标定不能进入长曝光主流程。

### 8.2 波形Monte Carlo怎样确定剩余未知量

`simulate_waveform_gls_calibration`生成2048条零信号与512条已知注入两镜记录，默认每条20 μs，信号倍率10000。每条都经过事件、SPE、ADC和互相关，结果是两张“记录×滞后”矩阵。

零信号前一半估计样本协方差；同一滞后差的矩阵对角线取平均，施加平稳Toeplitz结构以减少有限样本噪声。协方差再按比例 $\rho$ 向其对角阵收缩，候选 $\rho=10^{-4},10^{-3},10^{-2},.05$；特征值下限设为最大值的 $10^{-10}$，以便稳定求解。

第二个四分之一零信号与信号校正集共同选择“噪声标准差/正响应增益”最小的候选。信号校正集是512条中的第三个四分之一；用其平均投影除已知注入幅度，得到共同增益，乘回连续模板。均值标准误除均值给出 `response_relative_uncertainty`。前一半信号保留原始MC响应诊断，在线性主路径不再把其粗采样峰作为连续形状来源。

最后四分之一零信号没有参加收缩选择，用它的实际投影标准差与预测标准差之比平方缩放协方差。其有限尺度误差约为 $1/\sqrt{2(n-1)}$，$n$为这批记录数；保存为 `sigma_relative_uncertainty`。最后四分之一信号仅做响应检验，不再修改增益。

图中同时画连续理论、乘MC增益的模板、原始MC样本平均，明确后两者是否独立。使用了连续形状的主标定与连续理论不能仅凭峰形吻合宣称独立验证；独立新种子的零信号、不同倍率及块长记录才继续检验响应与误差。

仪器签名同时记录数值设置与模板/电荷/光学核内容。长曝光入口核对签名、源率和总背景率；改星等或SPE之后不能继续用旧标定。允许自定义短记录模板，不意味着主相位模型已经标定了该模板。

## 9. 把短块标定沿整个观测窗口累积起来

`run_sii_pipeline`先调用 `generate_uvw`，然后调用 `simulate_uv_observation(estimator='waveform_gls', waveform_calibration=...)`。内部先计算第4节的 $\overline P_a$，再调用同一个 `phase_template_bank` 和 `tracked_segment_precision`，而不是给所有基线复制零时延误差。

相位积分比可见度积分更密：每段1200个时间中点，重新计算源方向、基线到达时差和采样周期内的残差，再插值相位表的块标准差。基线间到达时差折回一个采样周期；相差整数滞后只平移峰，±200 ns窗口足够容纳±4 ns残差。该近似的截窗和节点收敛由相位留出、1200/2400点对照检查。

每条长曝光测量对应等时长块平均，所以方差必须平均块方差：

$$
\sigma_a^2=\frac{T_b}{T_a}\,
\frac{1}{N_t}\sum_{r=1}^{N_t}\sigma_b^2(\delta_a(t_r)),
\qquad T_a=\mathrm{segment\_s}\times\mathrm{nights}.
$$

$N_t$为相位积分点数，$\sigma_b$是标定块误差。不能改为平均信息量的倒数，否则测量会变成依赖相位的加权时间平均，而第4节期望仍是等时间平均。持续观测假定读出缓冲能覆盖几何移位；不会把孤立24 μs测试记录的边缘损失每隔24 μs再扣一遍。

然后只为整批数据抽一个标定增益 $g\sim\mathcal N(1,s_g^2)$，为每行抽独立标准正态 $z_a$，生成

$$
y_a=g\overline P_a+\sigma_a z_a.
$$

这正是写入 `visibility2_measured` 的量。`sigma_visibility2_stat` 和主路径的 `sigma_visibility2`均为独立统计误差；`calibration_relative_sigma`携带共同增益先验。逐点显示用的 `sigma_visibility2_calibration` 不作为独立方差加进合并权重。摘要中的单个sigma是各行的中位数，不能代替表内逐基线误差。

`sample_waveform_gls_visibility2`仍保留为固定模板、固定相位的标量统计工具，主整夜入口已经不再用它代替相位追踪。这里使用高斯与独立统计项是弱光、平稳线性条件下的长曝光近似；后文36镜性能另加共同噪声预算与噪声标定尺度对照。

## 10. 整理、保存再读取，必须还是同一个测量模型

`prepare_reconstruction_uv`只读取测量值、误差和几何，不读取真图或 `visibility2_true`。默认以120百万波长为网格宽度，把代表UV坐标舍入成组号。分组用于控制计算量，不是把整个格当成一个精确UV位置。

同组第 $a$ 条测量的独立方差为 $s_a^2$，令 $f_{ga}=s_a^{-2}/\sum_{b\in g}s_b^{-2}$。输出测量为 $y_g=\sum_{a\in g}f_{ga}y_a$，独立方差为 $1/\sum_{a\in g}s_a^{-2}$。与此同时，每个原始积分节点的权重变为 $f_{ga}\omega_{aq}$，坐标本身不变。这样分组后的预测和分组数据经历同一个线性算子。

`UvData.sampling`保存四个一维数组：节点u、节点v、所属输出组、节点权重。`u_lambda/v_lambda`只是代表坐标，可用于画图；存在sampling时，拟合不拿代表坐标代替节点。权重必须有限非负，每个输出组总和为1。

### 10.1 共享零点误差如何传递

旧解析观测入口允许 `baseline_zero_point_rms` 非零。此时同一基线跨时段共享一个零点随机数，而不同基线独立。输出增加 `baseline_zero_point_sigma`；合并使用 `sigma_visibility2_stat`作权重，把共享部分单独传播。

记 $A$为上述分组线性算子，$D$为原始独立方差对角阵，$B$把每个单位正态基线零点映射到对应测量，其元素就是该行的零点RMS。则分组协方差为

$$
C_{\rm group}=ADA^T+(AB)(AB)^T.
$$

代码直接累积小矩阵 $AB$，不先分配11340×11340矩阵。全零接口时完全不分配这项。共享零点0.05即使合并很多时段也保留约0.05的误差底限；同一基线进入不同组还留下非零非对角项。`UvData.sigma`改为该完整协方差的对角平方根，`covariance`保存完整矩阵。

当前全协方差高斯拟合可使用固定平滑强度；自动80/20留出选择平滑尚不支持跨训练/验证集的协方差，因而会明确拒绝 `smoothness='cv'`。其余每镜增益、每夜增益、透明度、NSB及未知残余时延漂移的非零长曝光似然尚未标定，两种长曝光入口都拒绝，默认零值照常运行。

### 10.2 用完整NPZ交付可复现的推断输入

`write_uv_data`保存版本号1、所有测量数组、上述sampling、可选协方差、共同增益先验、合并数和输入质量计数。`read_uv_data`以 `allow_pickle=False`读取，检查形状、有限性、正误差、组号、权重和协方差正定性，不以默认零值补掉缺失物理字段。

主notebook实际把四个源的 `UvData`各保存为 `*_uv.npz`，重新读回后才构造图像核并拟合。因此发布结果验证的是磁盘输入流程。最小调用为：

```python
uv = sii.prepare_reconstruction_uv(result.measurements, cell_mlambda=120.)
reco.write_uv_data('observation_uv.npz', uv)
uv = reco.read_uv_data('observation_uv.npz')
fit = reco.reconstruct_uv_data(uv, grid_size=32, fov_mas=.70,
    support_radius_mas=.32, starts=3, max_iter=8000)
```

`read_uv_measurements`遇到 `.npz`也转入完整读取。传统CSV入口只表示单点UV测量；一旦发现积分节点或共享误差字段，就拒绝并指向完整NPZ，而不静默换模型。发布的 `walkthrough_disk_measurements.csv`和`binary_measurements.csv`是浏览简表，完整推断分别用 `single_disk_uv.npz`和`binary_uv.npz`。简表不足以独立恢复积分观测，不能因列名相似就拿去替代。

## 11. 反演如何逐次计算一幅候选图的预测

`reconstruct_uv_data`先建立固定视场、像素坐标和圆形支撑。像素存的是该格归一化流量 $I_p$，不是未经面积归一的亮度密度；支撑外固定零，支撑内 $I_p\ge0$，总和为1。

对一次优化，UV节点不变，只是图像反复变化。若每次都对几十万个节点计算傅里叶会很慢。`power_sampling_kernel`用图像自相关的等价表达，把所有积分几何预计算为

$$
K_{g,\Delta}=\sum_{a,q\in g}f_{ga}\omega_{aq}
\cos\{2\pi[u_{aq}\Delta\theta_E+v_{aq}\Delta\theta_N]\},
\qquad m_g(I)=\sum_\Delta K_{g,\Delta}A_I(\Delta).
$$

$\Delta$标识像素间角位移，单位rad；$A_I$是图像与自身翻转的全线性卷积。`power_from_image`用 `fftconvolve`算自相关再矩阵乘核。该式与逐节点算 $|V|^2$再平均代数等价，不是插值补UV，也不是对观测测量再次自相关。32×32图对应63×63位移格，核大小为“合并测量数×3969”；预计超过1 GB则明确报错。

`_power_gradient`把测量损失的导数经核转置传回位移格，消除浮点非对称后，与图像作卷积并乘2，得到对像素流量的导数。这是自相关两端都依赖图像产生的因子2。预计算只改变速度，不改变曝光积分或误差定义。

## 12. 怎样比较预测与数据、选参数并判断优化结果

记测量向量为 $\boldsymbol y$，图像预测为 $\boldsymbol m$，统计协方差为 $C$。`statistical_loss`计算绝对高斯目标的一半：

$$
L=\tfrac12(g\boldsymbol m-\boldsymbol y)^TC^{-1}(g\boldsymbol m-\boldsymbol y)
+\frac{(g-1)^2}{2s_g^2}.
$$

独立情况按sigma逐点白化；有全协方差时用Cholesky和三角求解。权重不截断、不归一为均值1。共同增益在 `profile_calibration_gain`中解析消去：

$$
g_* =\frac{\boldsymbol m^TC^{-1}\boldsymbol y+s_g^{-2}}
{\boldsymbol m^TC^{-1}\boldsymbol m+s_g^{-2}}.
$$

$s_g=0$直接固定 $g=1$。对图像求导时用包络定理，不再对最优增益求额外导数。完整推导见[物理说明附录I](SII_PHYSICS_ZH.md#appendix-i)。噪声测量可为负或大于1，原值保留；只约束候选图的物理性质。

图像平滑项惩罚相邻像素亮度密度梯度。代码把用户强度乘 $(N-1)^4$，补偿网格改变时像素流量的尺度，$N$是每边像素数。默认 `smoothness='cv'`用固定随机划分的80%数据训练、20%验证，在0、0.0001、0.01中选择，不按真图相似度调参。每个候选只在训练集拟合增益；验证评分接收训练增益后验的秩一不确定度及其行列式项，不能在验证集重新校准以降低残差。

选择平滑后，在全数据上运行多个初值：集中、弥散、随机平滑图。默认流量参数化将有界非负变量除其总和；也保留softmax对照。优化用带解析梯度的L-BFGS-B，取目标最低的起点。本例最多8000次迭代；失败状态、每起点目标和停止消息都保存，不只保留好看的图。

程序另外计算单纯形驻点间隙 $G=\sum_p I_p\partial L/\partial I_p-\min_p\partial L/\partial I_p$。它衡量向梯度最小像素挪一点流量还能降低多少目标；通过阈值只是局部必要条件，不证明非凸问题全局最优。

优化结束后只允许不损失流量、不越出支撑的整数平移作居中，再重新计算保存图像对应的预测。真图直到拟合后才用于比较；配准只允许平移和180°中心反演，不周期卷绕。原因是这些变换无法由两镜平方可见度区分。峰数、图像误差、配准保留通量及重复模拟的稳定性用于描述解，不能把算法收敛等同于数据唯一重建了遮挡细节。

圆盘参数拟合与图像优化共享同一个平均观测模型，但候选只是一维直径。`profile_model_grid`对每个直径剖面共同增益，`profile_grid_interval`取最小点与 $\Delta\chi^2\le3.8414588$ 的区间，并标记截边和不连通。95%是常规单参数似然近似，最终通过重复模拟检验覆盖，不能从一次拟合推断准确率。

## 13. 最终36镜性能怎样与短波形连接，哪些验证回答哪些问题

前面的流程已经能从源生成并反演UV数据。`tools/validate_sii_performance.py`进一步使用全部36台实际共享ADC来学习响应与噪声，检查这一外推在固定范围内是否可用。它与主流程复用相位模板和段内方差积分，标定记录与图像演示各自独立。

三个时角各生成96条点源训练、192条零信号、96条圆盘留出，合计1152条24 μs阵列记录。点源和圆盘倍率300，零信号单独生成；点源真功率1用于训练，待测圆盘真值不用于校准。630条基线共享36通道，所以响应与噪声学习误差按记录级阵列统计，不把一条记录的630基线当作630条独立标定记录。原始共同时间0.027648 s，裁后0.026164224 s，不能乘36或630。

各相位留出组检查恢复偏差及误差尺度。得到的共同增益与噪声尺度随后用于整夜模型。`weak_light_covariance_budget`计算线性弱光Bartlett跨基线谱界，并用未分辨多模热光计数的精确二至四阶矩作另一对照；计数模型的界不冒充SPE滤波后全部四阶项的严格界。

性能表预定27个场景：星等2、4、6；直径0.08、0.16、0.32 mas；曝光1、3、6小时，每格500次，共13500次。直径候选0至0.64 mas，步长0.0005 mas；似然细化到0.00005 mas，并用非网格直径的直接模型检查插值。误差同时报告偏差、RMSE、区间覆盖、覆盖率二项区间、边界和多峰标记。

`compressed_profiles`为了批量模拟，以候选模型白化后的向量建立QR子空间，真值不参与选基。正交补只给所有候选同一个卡方常数，故似然比可在低维子空间计算；逐候选重构误差超过阈值会报错。每次抽一个共享增益、一个正值噪声尺度和子空间高斯噪声；剖面增益与固定增益两种分析使用同一批观测。拟合sigma还做噪声学习误差及弱光共同协方差的保守修正。这些是已说明的条件模型，不是模拟了未知真实天气或电子学漂移。

有限瞳面另由 `pupil_difference_quadrature`计算两独立均匀圆瞳入射点之差的分布，再由 `aperture_disk_power`对基线加该差向量求平均，量化8 m均匀瞳面对镜心近似的影响。真实遮挡、不同入瞳响应与滤光片角响应的联合分布未给全，故该对照没有冒称自动重现真实瞳面。

其他实际执行的检验各有明确作用，完整数表分别保留：

| 入口/结果目录 | 它检验什么 | 不能据此推出什么 |
|---|---|---|
| `thermal_mode_counts`、`joint_thermal_mode_counts` | 独立复高斯场→光强→泊松计数的二阶、三阶及边缘矩 | 没有生成整夜飞秒电场或纳秒SPE波形 |
| `validate_sii_observation.py` / `sii_observation` | 三镜共享波形、补偿、插值核收敛、24/96 μs曝光缩放 | 短块通过不证明整夜平稳性 |
| `validate_sii_tracking.py` / `sii_tracking` | 七时刻固定相干矩阵下更新几何优于固定初始时延 | 固定相干不是恒星随UV变化的结果 |
| 同脚本 `--source-case single_disk` / `sii_source_tracking` | 七时刻真实源相干，点源训练与圆盘留出分开 | 4410个36镜功率预测不等于4410组独立ADC |
| `test_sii_review_regressions.py` | 分数相位、完整数据往返、共享零点、非零噪声/抖动、窄通带 | 六类回归通过不取代实际仪器标定 |
| notebook的图像重复与网格对照 | 多种数值设置、随机观测下的解及驻点稳定性 | 局部解稳定不证明相位唯一或行星检出 |

## 14. 保留的零值接口、替代路径与规格书边界

当前实测SPE、电荷涨落、光追时间展宽保持启用。下面的零值只是缺少实测输入时关闭效应的约定，不是把这些效应测成了零。

| 接口 | 当前值 | 已实现范围与下一步需要 |
|---|---:|---|
| `dark_count_rate_hz` | 0 | 独立泊松暗计数可用；按实际汇总通道、温度、过压填入并重标定 |
| `electronic_noise_rms_mv` | 0 | 独立白噪声可用；有色谱或跨通道噪声需要对应模型 |
| `intrinsic_time_jitter_ns` | 0 | 独立高斯抖动已同时进入波形、连续模板和匹配带宽 |
| `microcell_recovery_time_ns` | 0 | 短波形有指数恢复；非线性长曝光相位标定尚未支持 |
| `adc_bits/adc_full_scale_mv` | 0/0 | 短波形有量化剪裁；线性解析相位模型不能直接用于非零ADC |
| 串扰、后脉冲概率 | 0 | 非零显式拒绝，单个总概率不足以决定时间/电荷分布 |
| 镜增益、夜增益、时延、透明度和NSB漂移 | 0 | 长曝光非零似然尚未标定，显式拒绝 |
| `baseline_zero_point_rms` | 0 | 旧解析路径支持共享基线零点及完整协方差传播；主波形路径仍需其专门标定 |

S17351规格书的结构化摘录为 [s17351_datasheet.json](../configs/sii/s17351_datasheet.json)，保留原PDF哈希。8通道、每通道33792微单元与main一致。25℃、过压8.5 V下每通道暗计数典型1.2 MHz、最大3 MHz，通过 `datasheet_dark_scenario`单独比较1或8通道汇总；默认仍关闭。典型串扰3%不会自动再叠加到已有实测电荷涨落里。PDF的PDE点用于核对完整曲线，不能凭典型点补齐真实噪声谱、恢复时间、延迟分布或整机ADC。详细条件与缺项见物理说明第13节。

旧 `estimator='analytic'`按星光/背景计数及解析SNR生成观测，属于替代比较路径。`with_matched_effective_bandwidth`先积分光学及本征抖动传递函数的四次模、散粒噪声占总谱比例的平方，把电子/量化噪声计入 $B_{\rm eff}$。标记 `optical_timing_in_effective_bandwidth=True` 后，解析观测不再乘一次电子学效率或光学衰减。旧带宽模型含采样箱sinc近似，主波形则在采样中心求SPE；不能把匹配带宽后的仪器再送给主波形当作另一段前端滤波。

`likelihood='legacy'`还保留历史幅度/Huber及权重缩放对照，会加旧零基线约束，不是本文的统计推断路径。无噪声数据须显式选择 `likelihood='noiseless'`；不能缺sigma却默认进行有统计意义的性能推断。

库中历史 `simulate_short_pair_waveforms`先在ADC时间箱内抽共享计数，再分别在两镜箱内展开事件时刻；它保存的是箱计数的相关，连续对事件的两端时刻不完全相同。因此它与主流程的 `simulate_hbt_primary_pe`具有不同的亚采样响应，不能用来替代本文标定。`make_fast_spe_template`提供测试用的双指数形状，正式例子读实测模板。`mean_recovery_fraction`计算均匀泊松照明下的平均恢复比例，仅作解析对照；实际启用恢复的电压仍逐微单元处理。这些历史/辅助入口存在于同一个模块，不表示notebook同时混用它们。

## 15. 从干净环境复现，并定位每一步的输出

仓库的 [requirements-sii-validated.txt](../requirements-sii-validated.txt)固定已验证的直接依赖版本，包含ipykernel；它不是操作系统、BLAS及全部传递依赖的完整环境镜像。在仓库根目录运行：

```powershell
python -m pip install -r requirements-sii-validated.txt
python -m pytest tests -k sii -q
python tools/build_sii_science_notebook.py
python tools/execute_notebook.py notebooks/sii_complete_waveform_report.ipynb --cwd . --timeout 7200
python tools/check_sii_science_artifacts.py
```

生成脚本会清空旧notebook输出；执行脚本按顺序运行全部38个单元中的20个代码单元，失败时保留已执行输出并返回错误。完整波形验证和图像重复需要较长时间。修改数值代码或输入后应从头重算，不能只运行末尾单元覆盖摘要哈希。

若只运行一条基准圆盘链，下面把前述调用接成可直接执行的脚本。它生成完整观测并保存后读回，不执行27场景性能扫描及附加检验：

```python
from pathlib import Path
import sys

root = Path.cwd()  # 仓库根目录
sys.path.insert(0, str(root / 'python'))
import sii_unified as sii
import sii_reconstruction as reco

instrument = sii.Instrument.from_repository(root)
observation = sii.Observation(hours_per_night=6., segment_s=1200,
                              visibility_subsamples_per_segment=9)
source = sii.BinarySource(ab_magnitude=2., primary_diameter_mas=.16)
calibration, diagnostics = sii.simulate_waveform_gls_calibration(
    instrument, source_ab_magnitude=2., null_records=2048,
    signal_records=512, covariance_shrinkage='auto', seed=20260905)
result = sii.run_sii_pipeline(
    root / 'configs/arrays/lact36_20260906.input', source, observation,
    instrument, source_case='single_disk', estimator='waveform_gls',
    waveform_calibration=calibration, do_reconstruction=False, seed=20261306)
uv = sii.prepare_reconstruction_uv(result.measurements, cell_mlambda=120.)
reco.write_uv_data('observation_uv.npz', uv)
uv = reco.read_uv_data('observation_uv.npz')
fit = reco.reconstruct_uv_data(uv, grid_size=32, fov_mas=.70,
    support_radius_mas=.32, starts=3, max_iter=8000,
    smoothness='cv', smoothness_candidates=(0.,1e-4,.01), seed=20261406)
print(fit.metrics['optimizer_success'], fit.metrics['stationarity_passed'])
```

主要结果位于 [validation/sii_science](../validation/sii_science)：

| 顺着本文的数据流 | 文件 |
|---|---|
| 实际输入与来源 | `parameters.csv`、`walkthrough_parameters.json`、`summary.json` |
| 理论光学相关 | `walkthrough_ideal_hbt.csv` |
| 同一事件到ADC到相关 | `walkthrough_record.npz`、`walkthrough_pair_correlation.csv`、`walkthrough_pair_result.json` |
| 标定峰与协方差 | `response.csv`、`covariance.csv`、`waveform_calibration.npz` |
| 独立波形检验 | `independent_null.csv`、`injection_scale.csv`、`block_duration.csv` |
| 实际阵列与浏览测量 | `array_baselines.csv`、`walkthrough_disk_measurements.csv`、`binary_measurements.csv` |
| 完整反演输入 | `single_disk_uv.npz`、`binary_uv.npz`、`ellipse_uv.npz`、`transit_uv.npz` |
| 参数与图像推断 | `diameter_profile.csv`、`diameter_coverage.csv`、`reconstruction.csv`、各源`*_image.npy` |
| 数值及重复稳定性 | `time_quadrature.csv`、`convergence.csv`、`image_alignment.csv`、`image_repeats.csv` |

最终36镜条件性能另位于 [validation/sii_performance](../validation/sii_performance)，包含逐记录630基线投影 `raw_projections.npz`、相位留出 `raw_phase_holdout.csv`、性能表 `performance.csv`、弱光与瞳面对照预算和规格书情景。主摘要绑定这些子验证摘要的哈希。性能CSV按实际字节记哈希，`.gitattributes`固定其CRLF，避免检出换行造成假不一致。

检查脚本核对公式分隔符、文档链接、notebook顺序执行、输入/代码哈希，从保存ADC重算相关和GLS，并重算完整阵列真功率；完整NPZ往返由回归测试和notebook实际读回覆盖。代码主体在 [sii_unified.py](../python/sii_unified.py)，共享阵列事件在 [sii_observation.py](../python/sii_observation.py)，图像推断在 [sii_reconstruction.py](../python/sii_reconstruction.py)，相位及性能在 [sii_performance.py](../python/sii_performance.py)，独立参考在 [sii_validation.py](../python/sii_validation.py)。这些文件边界服从前述同一数据流。

新的实测输入到来后，先核对来源并更新清单，再重算事件/响应标定、相位表、观测、推断和覆盖率。只改变注释或PDF摘录不会自动进入数值计算。原始备份标签保留为 `backup/sii-gls-before-fixes-20260905`；研究计划与临时审查文件留在本地，不属于本说明的用户复现入口。
