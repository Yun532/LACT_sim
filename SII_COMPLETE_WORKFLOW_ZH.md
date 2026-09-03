# LACT 单像素恒星强度干涉：从源模型到图像重建的完整计算流程

本文说明当前统一版本实际计算的物理量、参数来源和数值结果。全文按因果顺序组织：一个参数只在定义后使用，同一个物理量始终使用同一符号。

```text
源的归一化空间亮度 I(l,m) 与总星光通量
  → 连续理论复可见度 V(u,v) 和功率可见度 P(u,v)=|V|²
  → 32 台望远镜随时间变化的 (u,v,w) 采样坐标
  → 每台望远镜的恒星、NSB 和暗计数率
  → 光谱相干尺度 τ_eff
  → 光学时间核、SPE、采样和电子学共同决定的有效相关带宽 B_eff
  → 微秒级双镜完整波形，验证相关器和时间响应
  → 小时级 P(u,v) 充分统计量及其不确定度
  → 全阵列稀疏 UV 数据
  → 无相位、带约束的图像重建
```

对应实现如下：

- 完整可执行 notebook：[`notebooks/lact_sii_paper_simulation.ipynb`](notebooks/lact_sii_paper_simulation.ipynb)
- 统一模拟模块：[`python/sii_unified.py`](python/sii_unified.py)
- 独立重建模块：[`python/sii_reconstruction.py`](python/sii_reconstruction.py)
- 自动测试：[`tests/test_sii_unified.py`](tests/test_sii_unified.py)

本文数值和图片来自 `intensity-interferometry-unified` 分支的完整 notebook 重算，随机实验使用固定种子。图片不是手工示意图。

为避免缩写造成歧义：HBT 指 Hanbury Brown–Twiss 二阶强度相关；PE 或 p.e. 指一个光电子；SPE 指单光电子脉冲；NSB 指夜天光背景；ADC 指模数转换器；DC 指 Davies–Cotton 反射镜结构；PDE 指光子探测效率；PSF 指点扩散函数；FWHM 指半高全宽；mas 指毫角秒，$1\,\mathrm{mas}=10^{-3}$ 角秒。SNR 指信噪比，即信号期望值与其标准差之比。

---

## 0. 先明确模拟的三个层次

完整模拟包含三个不同层次，不能把它们混成一个步骤。

1. **天空真值层**：给定恒星在天空中的强度分布，计算任意空间频率处应有的理论 $P(u,v)=|V(u,v)|^2$。
2. **几何采样层**：望远镜坐标、源方向和观测时间决定理论 UV 平面上哪些点被采样。这个坐标本身与光学结构和电子学无关。
3. **测量层**：星光、夜天光背景、光学到达时间展宽、SiPM、SPE 波形、采样和电子学决定每个 UV 点测得多准，以及未校正时是否被衰减。

因此，“UV 点在哪里”和“这个点的数值及误差是多少”是两件事：前者由阵列几何决定，后者才由源和仪器决定。

### 0.1 全文统一使用的核心符号

下表先定义会贯穿全文的量；局部参数仍会在首次出现的步骤中解释。

| 符号 | 定义 | 单位 |
|---|---|---|
| $I(l,m)$ | 天空切平面上的归一化空间亮度 | $\mathrm{rad}^{-2}$ |
| $V(u,v)$ | 归一化复可见度 | 无量纲 |
| $P(u,v)$ | 功率可见度，定义为 $|V(u,v)|^2$ | 无量纲 |
| $\lambda_0$ | 窄带滤光片中心波长 | m |
| $\Delta\lambda$ | 滤光片的等效波长宽度；当前用 FWHM 近似 | m |
| $\Delta\nu$ | 与 $\Delta\lambda$ 对应的光学频率带宽 | Hz |
| $\tau_{\rm eff}$ | 未被电子学分辨的 HBT 峰的有效相干面积 | s |
| $K_{\rm opt}(t)$ | 单镜归一化光学到达时间核 | $\mathrm{s}^{-1}$ |
| $h_{\rm SPE}(t)$ | 一个光电子产生的实测单 PE 电压脉冲 | V |
| $B_{\rm eff}$ | 实际时间响应在最优相关器中的有效相关带宽 | Hz |
| $r_{\star,i}$ | 第 $i$ 台望远镜探测到的恒星光电子率 | $\mathrm{s}^{-1}$ |
| $r_{{\rm tot},i}$ | 第 $i$ 台望远镜探测到的总光电子率 | $\mathrm{s}^{-1}$ |
| $T_{\rm rec}$ | 一次显式短波形的有效记录时长 | s |
| $T_{\rm seg}$ | 一个长曝光统计时间段的积分时长 | s |
| $F_{\rm EN}$ | SiPM 单 PE 电荷涨落造成的 SNR 惩罚因子 | 无量纲 |

其中 $\tau_{\rm eff}$ 描述飞秒量级的光场相干，$K_{\rm opt}$ 描述纳秒量级的光路传播展宽，$h_{\rm SPE}$ 描述约百纳秒的电子脉冲。三者的来源和用途完全不同，后文会分别推导。

---

## 1. 生成源的空间强度模型

### 1.1 形态和总亮度分开描述

令 $l$ 和 $m$ 分别表示相对相位中心在东西、南北方向的小角偏移，单位为弧度。在小视场近似下，可分别理解为 $\Delta\alpha\cos\delta$ 和 $\Delta\delta$，其中 $\alpha$ 是赤经，$\delta$ 是赤纬。

用 $I(l,m)$ 表示源的**归一化形态**，规定

$$
\iint I(l,m)\,dl\,dm=1.
$$

这样，$I$ 只决定源长什么样；总 AB 星等 $m_{\rm AB}$ 单独决定到达望远镜的总星光通量。这个分离避免把“形态归一化”和“光子数”重复计算。

当前默认源由两颗均匀圆盘组成：

$$
I(l,m)=\frac{I_1(l-l_1,m-m_1)+r_F I_2(l-l_2,m-m_2)}{1+r_F},
$$

其中 $(l_k,m_k)$ 是第 $k$ 颗星的中心位置，$I_k$ 是积分归一化为 1 的圆盘亮度，$r_F=F_2/F_1$ 是次星和主星的总流量比。默认案例使用：

- 总亮度 $m_{\rm AB}=2$；
- 两星中心分离角 $s=0.20\,\mathrm{mas}$；
- 位置角 $\mathrm{PA}=35^\circ$，从北向东计算；
- 流量比 $r_F=0.55$；
- 主星角直径 $0.060\,\mathrm{mas}$，次星角直径 $0.040\,\mathrm{mas}$。

离散天图只用于显示真值和检验重建；理论可见度使用解析圆盘公式，因此不会受天图像素大小限制。

![理论双星源天图](docs/sii_workflow_figures/01_theoretical_source_sky.png)

**图 1。** 默认双星的归一化空间亮度。横轴为 $\Delta\alpha\cos\delta$，纵轴为 $\Delta\delta$，单位均为 mas。

---

## 2. 从源模型得到连续理论 UV 平面

### 2.1 归一化复可见度

van Cittert–Zernike 定理把空间亮度和复可见度联系起来：

$$
V(u,v)=\iint I(l,m)
\exp[-2\pi i(ul+vm)]\,dl\,dm.
$$

由于第 1 节已经令 $\iint I\,dl\,dm=1$，此处不再需要分母；若使用未归一化天图，则应再除以总亮度积分。$u$ 和 $v$ 是以波长为单位的投影基线，因此是无量纲空间频率。

强度干涉不直接给出 $V$ 的相位，而是测量

$$
P(u,v)\equiv |V(u,v)|^2.
$$

全文后续统一用 $P$ 表示功率可见度。

### 2.2 均匀圆盘和双星解析式

角直径为 $\theta$ 的均匀圆盘在径向空间频率 $q_{\rm uv}=\sqrt{u^2+v^2}$ 处的可见度为

$$
V_{\rm disk}(q_{\rm uv};\theta)=
\frac{2J_1(\pi\theta q_{\rm uv})}{\pi\theta q_{\rm uv}},
$$

其中 $J_1$ 是第一类一阶贝塞尔函数，$\theta$ 必须用弧度表示。

两颗圆盘组成的双星可见度为

$$
V(u,v)=
\frac{
V_1(q_{\rm uv})e^{-2\pi i(ul_1+vm_1)}+
r_F V_2(q_{\rm uv})e^{-2\pi i(ul_2+vm_2)}
}{1+r_F}.
$$

分离角主要决定 UV 平面中条纹的间距，位置角决定条纹方向，恒星直径决定长基线处的圆盘包络，流量比决定条纹调制度。

### 2.3 为什么要同时画正、负基线

真实亮度 $I(l,m)$ 满足 Hermitian 对称性：

$$
V(-u,-v)=V^*(u,v),
$$

所以

$$
P(-u,-v)=P(u,v).
$$

这是关于原点的中心对称。程序把每个 $(u,v)$ 和 $(-u,-v)$ 用相同显著度绘制，因为两者表示同一条无序基线的共轭点，不是两个独立测量。图上有时看起来更像左右对称，是阵列布局投影、源赤纬、时角范围和双星位置角共同形成的视觉效果；地球自转让各基线沿轨迹移动，但不会把数学上的中心对称改成单纯镜面对称。

![理论 UV 功率平面](docs/sii_workflow_figures/02_theoretical_uv_power.png)

**图 2。** 连续理论 $P(u,v)$。此图只由源模型决定，还没有加入阵列采样或仪器噪声。

---

## 3. 用 32 台望远镜和观测时间生成 UVW

### 3.1 物理基线

令第 $i$ 台望远镜在站点局部东、北、天顶坐标系中的位置为

$$
\mathbf r_i=(E_i,N_i,U_i),
$$

单位为 m。望远镜 $i$ 到 $j$ 的有向物理基线是

$$
\mathbf B_{ij}=\mathbf r_j-\mathbf r_i=(B_E,B_N,B_U).
$$

32 台望远镜共有

$$
N_{\rm bl}=\frac{32\times31}{2}=496
$$

条不重复基线；$N_{\rm bl}$ 表示基线数。

### 3.2 源方向和地球自转

令 $H$ 为源时角，定义为

$$
H(t)=\mathrm{LST}(t)-\alpha,
$$

其中 $\mathrm{LST}$ 是当地恒星时，$\alpha$ 是源赤经，$t$ 是观测时刻。站点纬度记作 $\phi$，源赤纬记作 $\delta$。

在东、北、天顶坐标系中，指向源的单位向量为

$$
\hat{\mathbf s}=
\begin{bmatrix}
\cos\delta\sin H\\
\sin\delta\cos\phi-\cos\delta\cos H\sin\phi\\
\sin\delta\sin\phi+\cos\delta\cos H\cos\phi
\end{bmatrix}.
$$

天空切平面上的两个正交单位向量定义为

$$
\hat{\mathbf u}=
\begin{bmatrix}
\cos H\\
\sin\phi\sin H\\
-\cos\phi\sin H
\end{bmatrix},
$$

$$
\hat{\mathbf v}=
\begin{bmatrix}
-\sin\delta\sin H\\
\cos\delta\cos\phi+\sin\delta\cos H\sin\phi\\
\cos\delta\sin\phi-\sin\delta\cos H\cos\phi
\end{bmatrix}.
$$

于是以 m 为单位的投影坐标为

$$
u_m=\mathbf B_{ij}\cdot\hat{\mathbf u},\qquad
v_m=\mathbf B_{ij}\cdot\hat{\mathbf v},\qquad
w_m=\mathbf B_{ij}\cdot\hat{\mathbf s}.
$$

选择滤光片中心波长 $\lambda_0$ 后，理论可见度所需的无量纲空间频率为

$$
u=\frac{u_m}{\lambda_0},\qquad
v=\frac{v_m}{\lambda_0}.
$$

$w_m$ 是沿视线方向的基线分量，不进入二维天图的 $P(u,v)$，但决定两台望远镜的几何到达时延：

$$
\tau_g=\frac{w_m}{c},
$$

其中 $c=299\,792\,458\ \mathrm{m\,s^{-1}}$ 是真空光速。

### 3.3 观测时间增加的作用

随着 $H(t)$ 改变，同一物理基线在 UV 平面上画出地球自转轨迹。更长的单夜观测同时带来两种收益：

1. 轨迹更长，方向覆盖通常更完整；
2. 对同一统计量累积更久，随机误差按 $T^{-1/2}$ 下降。

但低高度角处的大气消光、背景、滤光片角响应和光学 PSF 可能变差，因此程序设置最低高度角门槛，而不是无条件使用地平线附近数据。跨夜重复相同时角主要提高统计量，不会创造全新的 UV 轨迹。

![32 台望远镜的累计 UV 覆盖](docs/sii_workflow_figures/03_uv_coverage_cumulative_time.png)

**图 3。** 不同观测时长下 32 台望远镜的累计几何 UV 覆盖。该图只表示采样位置；正、负基线显著度相同。

---

## 4. 从 AB 星等得到单镜光电子率

这一节先确定每台望远镜每秒收到多少恒星和背景光电子。只有先有光电子率，后面才能推导相关对数和 SNR。

### 4.1 光学频带

当前模拟场景选择的中心波长为

$$
\lambda_0=400\ \mathrm{nm},
$$

并假设一个等效宽度为

$$
\Delta\lambda=2\ \mathrm{nm}.
$$

这里必须区分“LACT 当前已有硬件”和“SII 性能研究的场景参数”：main 提供了镜面、滤光片透过率和 PDE 曲线在 400 nm 附近的响应，但没有提供一片已经确定用于恒星强度干涉、等效带宽恰好为 2 nm 的实物窄带滤光片。因此 $\lambda_0=400\ \mathrm{nm}$、$\Delta\lambda=2\ \mathrm{nm}$ 目前是可替换的 SII 基准假设，不是 LACT 已定型滤光片的实测结论。

SII 常使用 nm 量级或更窄的光学通道，原因包括延长光学相干时间、抑制 NSB、减轻宽带色散和允许多个独立光谱通道分别相关。但“单通道越窄，理想 SNR 就越高”并不成立：在平坦连续谱和散粒噪声极限下，窄带使 $r_\star\propto\Delta\nu$ 下降，同时使 $\tau_{\rm eff}\propto1/\Delta\nu$ 上升，因此 $r_\star\tau_{\rm eff}$ 近似不变。真实系统中 NSB、谱线、PDE、电子噪声、计数率饱和和滤光片角响应会打破这种理想抵消，所以最终带宽需要扫描优化，不能把 2 nm 当成普适最优值。

窄带近似下，波长宽度对应的光学频率带宽为

$$
\Delta\nu\simeq
\left|\frac{d\nu}{d\lambda}\right|_{\lambda_0}\Delta\lambda
=\frac{c\,\Delta\lambda}{\lambda_0^2}.
$$

代入当前参数得到

$$
\Delta\nu=3.7474\times10^{12}\ \mathrm{Hz}
=3.7474\ \mathrm{THz}.
$$

$\Delta\nu$ 是**光学频带宽度**，不是 ADC 或相关器带宽。它同时进入第 4.2 节的光子率和第 5.2 节的相干面积；附录 A 给出为什么 $\int|g^{(1)}|^2d\tau\propto1/\Delta\nu$ 的完整推导。

### 4.2 AB 星等和恒星光电子率

AB 星等 $m_{\rm AB}$ 对应的谱能流密度为

$$
F_\nu=3631\ \mathrm{Jy}\times10^{-0.4m_{\rm AB}},
$$

其中 $1\ \mathrm{Jy}=10^{-26}\ \mathrm{W\,m^{-2}\,Hz^{-1}}$。中心波长处单个光子的能量为

$$
E_\gamma=h\nu_0=\frac{hc}{\lambda_0},
$$

其中 $h$ 是普朗克常数，$\nu_0=c/\lambda_0$ 是中心光频率。因此入射光子谱密度为

$$
n_\nu=\frac{F_\nu}{E_\gamma},
$$

单位是 $\mathrm{s^{-1}\,m^{-2}\,Hz^{-1}}$。

令 $A_{\rm col}$ 表示遮挡后用于归一化的收集面积，$\eta_{\rm det}$ 表示从该面积到中心像素产生光电子的总探测效率。窄带内近似不变时，第 $i$ 台望远镜的恒星光电子率为

$$
r_{\star,i}=A_{{\rm col},i}\,\eta_{{\rm det},i}\,n_\nu\,\Delta\nu.
$$

当前 400 nm 单像素光追给出中心像素有效探测面积 $5.92276\ \mathrm{m^2}$；再乘当前上游大气/未单列传输因子 $0.7836336$，得到净面积效率乘积

$$
A_{\rm col}\eta_{\rm det}=4.64128\ \mathrm{m^2}.
$$

代码为兼容 main 配置，将它写成

$$
A_{\rm col}=24.57686\ \mathrm{m^2},\qquad
\eta_{\rm det}=0.188847,
$$

两种写法的乘积完全相同。默认 $m_{\rm AB}=2$ 时得到

$$
r_\star=2.01548\times10^8\ \mathrm{s^{-1}}
=201.548\ \mathrm{MHz}.
$$

这里的 MHz 表示“每秒百万个光电子”，不是电子带宽。

### 4.3 NSB、暗计数和总计数率

令 $r_{{\rm NSB},i}$ 为单像素夜天光背景光电子率，$r_{{\rm dark},i}$ 为 SiPM 暗计数率。总光电子率为

$$
r_{{\rm tot},i}=r_{\star,i}+r_{{\rm NSB},i}+r_{{\rm dark},i}.
$$

当前 main 的宽带相机响应和暗夜天空模型积分给出

$$
r_{\rm NSB}=70.5275\ \mathrm{MHz},
$$

这个数并不是经过假设 2 nm SII 窄带片后的 NSB。当前 notebook 将它保留为保守背景压力值，所以时间响应研究仍有效，但光谱灵敏度不是最终自洽硬件预测。最终版本需要定义同一个总光谱响应 $T(\lambda)$，并用它分别积分恒星和 NSB；在此之前，不应把 70.5275 MHz 称为“2 nm 通道的背景率”。

暗计数因尚无本型号实测输入而暂设为

$$
r_{\rm dark}=0.
$$

所以默认单镜总率为

$$
r_{\rm tot}=272.075\ \mathrm{MHz}.
$$

恒星光子、NSB 和暗计数都按 Poisson 点过程产生，因此每台望远镜自身的光子数涨落已包含。NSB 和暗计数增加分母中的总散粒噪声，但不产生跨望远镜的天体 HBT 信号。

---

## 5. 光学频带如何决定 HBT 信号面积

### 5.1 从二阶相关到功率可见度

令 $J_i(t)$ 是第 $i$ 台望远镜的瞬时探测强度，$\delta J_i(t)=J_i(t)-\langle J_i\rangle$ 是强度涨落。使用 $J$ 是为了不与第 1 节的天空亮度 $I(l,m)$ 混淆。归一化二阶相关定义为

$$
g_{ij}^{(2)}(\tau)-1=
\frac{\langle\delta J_i(t)\,\delta J_j(t+\tau)\rangle}
{\langle J_i\rangle\langle J_j\rangle},
$$

其中 $\tau$ 是两路信号的相对时间延迟。对混乱热光，Siegert 关系给出

$$
g_{ij}^{(2)}(\tau)-1
=P(u,v)\,|g^{(1)}(\tau)|^2
$$

的单偏振理想形式，其中 $g^{(1)}(\tau)$ 是归一化时间一阶相干函数。实际系统把多个不可区分偏振模式和非矩形光谱合并后，最稳定的量不是无法被 ns 电子学直接分辨的峰高，而是峰面积。

### 5.2 有效相干面积 $\tau_{\rm eff}$

本文定义

$$
\tau_{\rm eff}
\equiv p_{\rm pol}
\int_{-\infty}^{\infty}|g^{(1)}(\tau)|^2\,d\tau
=\frac{p_{\rm pol}s_{\rm spec}}{\Delta\nu}.
$$

这里严格出现的是 $1/\Delta\nu$，不是直接的 $1/\Delta\lambda$。二者只通过窄带近似

$$
\frac{1}{\Delta\nu}\simeq
\frac{\lambda_0^2}{c\,\Delta\lambda}
$$

相互转换。直观上，频谱越宽，包含的不同光频越多，它们在延迟轴上越快失相干，所以相关峰越窄；频谱越窄，相关峰越宽。附录 A 从归一化矩形频谱的傅里叶变换和 Parseval 定理完整推导这个关系。

这里：

- $p_{\rm pol}$ 是偏振稀释因子；当前探测两种互不相干偏振，总计数率不分偏振，因此 $p_{\rm pol}=0.5$；
- $s_{\rm spec}$ 是实际或假设的归一化谱形相对理想矩形频带的峰面积修正；当前基准模型假设 $s_{\rm spec}=0.842$，它尚不是 LACT 实物 SII 滤光片的标定值；
- $\Delta\nu$ 是第 4.1 节已经定义的光学频率带宽。

代入当前值：

$$
\tau_{\rm eff}
=\frac{0.5\times0.842}{3.7474\times10^{12}}
=1.12344\times10^{-13}\ \mathrm{s}
=112.344\ \mathrm{fs}.
$$

代码属性 `coherence_area_s` 就是这里唯一的 $\tau_{\rm eff}$。它有时间单位，是“理想 HBT 峰的有效积分面积”，不是 DC 光学到达时间 RMS，也不是 SPE 脉冲宽度。

### 5.3 物理相关对率

先说明“超额相关”是什么意思。若两台望远镜接收到完全独立的光子流，它们仍会偶然在同一时间窗内各记录到光电子，这形成一个很大的无相关本底。热光的光子聚束使两路强度涨落在正确几何延迟处多出一个很小的正协方差：

$$
\langle\delta I_i(t)\,\delta I_j(t+\tau_g)\rangle>0,
$$

其中 $\delta I_i=I_i-\langle I_i\rangle$ 是望远镜 $i$ 的去均值强度涨落，$\tau_g$ 是第 3.2 节定义的几何时延。“HBT 超额相关”就是这部分**高于两路独立光子偶然符合本底的相关量**；它不是凭空增加新的光子，也不是指某个光子被两台望远镜同时探测。

校正几何时延后，这个超额协方差可以等价写成恒星热光产生的额外相关光电子对率：

$$
R_{{\rm pair},ij}
=r_{\star,i}r_{\star,j}\,\tau_{\rm eff}\,P(u_{ij},v_{ij}).
$$

$R_{{\rm pair},ij}$ 的单位是 $\mathrm{s^{-1}}$。默认案例在一条过中天基线上有 $P\simeq0.149476$，从而

$$
R_{\rm pair}=682.15\ \mathrm{s^{-1}}.
$$

在短波形 Monte Carlo 中，程序用稀疏“共享相关对”实现这个二阶矩；`hbt_pair_id` 只是模拟标签，真实仪器不能逐个判断哪个光电子属于 HBT 超额项。真实测量只能从大量样本的平均协方差中估计这部分微小超额。

这条公式是短波形和长曝光模拟共同的物理起点。后面两条路径必须从同一个 $R_{\rm pair}$ 或等价的 $\tau_{\rm eff}$ 出发，才能保证闭合。

---

## 6. 仪器时间响应和有效相关带宽

第 5 节给出了 HBT 信号的总面积。本节回答：这个面积经过 DC 光学时间弥散、SiPM 脉冲和 ADC 后在时间轴上如何展开，以及它的统计信息还能保留多少。

### 6.1 单镜光学到达时间核

令 $K_{\rm opt}(t)$ 是一个瞬时入射光子到达单像素后产生光电子的归一化延迟概率密度，满足

$$
\int K_{\rm opt}(t)\,dt=1.
$$

当前核来自 main 的 400 nm、纯平行光、中心 24.4 mm 单像素、100 万光子完整光追。光追包含镜面分面、遮挡、PSF、相机几何、集光器、滤光片和 PDE。公共传播时间已减去，只保留相对展宽。当前得到

$$
\sigma_{\rm opt}=0.548827\ \mathrm{ns},
$$

其中 $\sigma_{\rm opt}$ 是 $K_{\rm opt}(t)$ 的标准差。两台相同望远镜的相对延迟核是两个单镜核的互相关；若只用 RMS 近似，其宽度为

$$
\sigma_{\rm pair}\simeq\sqrt{2}\,\sigma_{\rm opt}=0.776\ \mathrm{ns}.
$$

因此 DC 结构确实影响 HBT 峰形：它把相关峰展宽并降低未校正的零延迟峰高。但归一化光学核不改变第 5 节的总相关面积。

当前只使用轴上响应。不同天顶角、方位角下的光学 PSF、到达时间核和滤光片角响应尚无 main 标定文件，因此接口保留、修正量当前设为 0；这不能被表述为“已证明方向无关”。

### 6.2 实测单 PE 波形和 SiPM 电荷

令 $h_{\rm SPE}(t)$ 表示一个光电子对应的实测电压脉冲模板。当前模板来自 537 个干净单 PE 脉冲，时间支撑约从 $-40\ \mathrm{ns}$ 延伸到 $+170$ 至 $183\ \mathrm{ns}$。支撑长度表示脉冲尾部持续多久，不等于可用相关带宽的简单倒数。

每个光电子还抽取归一化电荷因子 $a_k$。规定 $\langle a_k\rangle=1$，并定义

$$
F_{\rm EN}\equiv\sqrt{\langle a_k^2\rangle}.
$$

$F_{\rm EN}$ 是单 PE 电荷涨落的 SNR 惩罚因子；SNR 已在文首定义为信号期望值与标准差之比。当前 537 个实测电荷样本给出

$$
F_{\rm EN}=1.016142.
$$

微单元恢复使用指数模型。若同一微单元距上次击中的间隔为 $\Delta t_{\rm cell}$，恢复电荷比例为

$$
f_{\rm rec}=1-\exp\left(-\frac{\Delta t_{\rm cell}}{\tau_{\rm rec}}\right),
$$

其中 $\tau_{\rm rec}$ 是微单元恢复时间，当前暂取 $10\ \mathrm{ns}$。一个像素有 270336 个微单元；默认总率下平均占用极低，所以该效应很小，但短波形逐击中模拟仍保留它。

### 6.3 ADC 时间网格

令 $\Delta t_{\rm samp}$ 为 ADC 采样间隔。main 当前给出

$$
\Delta t_{\rm samp}=4\ \mathrm{ns}.
$$

因此采样率和 Nyquist 频率分别为

$$
f_s=\frac{1}{\Delta t_{\rm samp}}=250\ \mathrm{MS/s},
$$

$$
f_{\rm Nyq}=\frac{f_s}{2}=125\ \mathrm{MHz}.
$$

$f_{\rm Nyq}$ 是任何数字相关器可利用的最高频率上限。main 尚未给出独立标定的前端模拟带宽、ADC 位数、满量程和电子噪声，因此当前计算中 ADC 量化与附加电子噪声关闭；接口保留，数值均为 0。

短记录不能只等于 170 ns SPE 长度，否则边界脉冲会被截断，也不足以统计相关波动。当前显式短波形使用有效记录时长

$$
T_{\rm rec}=2\ \mu\mathrm{s},
$$

对应 500 个 4 ns 采样点，并在有效区两侧额外生成光电子以覆盖完整 SPE 尾部。若目的是直接看物理强度干涉信号，2 µs 仍然太短；它的作用是检查链路和统计闭合。

### 6.4 从时间核推导 $B_{\rm eff}$

先定义光学时间核的傅里叶变换

$$
H_{\rm opt}(f)=\int K_{\rm opt}(t)e^{-2\pi ift}\,dt.
$$

对两台相同望远镜，HBT 交叉谱的光学时间响应为 $|H_{\rm opt}(f)|^2$。令 $S_{\rm shot}(f)$ 表示由恒星、NSB 和暗计数产生的散粒噪声功率谱，$S_{\rm add}(f)$ 表示与光子流无关的附加电子噪声功率谱，并定义散粒噪声占比

$$
\eta_{\rm shot}(f)=
\frac{S_{\rm shot}(f)}{S_{\rm shot}(f)+S_{\rm add}(f)}.
$$

$\eta_{\rm shot}(f)$ 介于 0 和 1；当前因 $S_{\rm add}=0$，所以 $\eta_{\rm shot}(f)=1$。

对于使用已知峰形的匹配相关器，本文把有效相关带宽定义为

$$
B_{\rm eff}\equiv
\int_0^{f_{\rm Nyq}}
|H_{\rm opt}(f)|^4\eta_{\rm shot}^2(f)\,df.
$$

这个定义的含义是：一个带宽为 $B_{\rm eff}$ 的理想矩形、散粒噪声受限相关器，与当前真实响应具有相同的 $\mathrm{SNR}^2$。四次方来自两步：两镜交叉谱先得到 $|H_{\rm opt}|^2$，匹配滤波的 $\mathrm{SNR}^2$ 再对信号谱的平方积分。

当前 $K_{\rm opt}$ 是光追得到的离散混合分布。notebook 先计算 $H_{\rm opt}(f)$，再从 0 到 $125\ \mathrm{MHz}$ 数值积分上式，得到

$$
\boxed{B_{\rm eff}=110.910\ \mathrm{MHz}}.
$$

所以 110.91 MHz 不是手工选择或拟合参数，而是以下输入共同决定的计算结果：

1. 4 ns 采样给出的积分上限 $f_{\rm Nyq}=125\ \mathrm{MHz}$；
2. 100 万光子光追给出的 $H_{\rm opt}(f)$；
3. 当前附加电子噪声为 0，因此 $\eta_{\rm shot}(f)=1$。

等效地说，真实光学时间响应在这段频带内保留了约

$$
\frac{B_{\rm eff}}{f_{\rm Nyq}}=0.8873
$$

的匹配统计带宽。notebook 另外给出的 `rectangular_optical_retention=0.9407` 是 $|H_{\rm opt}|^2$ 的频带平均值，描述未匹配矩形零延迟峰高；它与积分 $|H_{\rm opt}|^4$ 定义的 $B_{\rm eff}/f_{\rm Nyq}$ 不是同一个量。

第 8.2 节将从计数统计推导出 $\mathrm{SNR}\propto\sqrt{B_{\rm eff}}$。因此已知且被匹配处理的 DC 时间展宽，主要通过降低 $B_{\rm eff}$ 影响长曝光结果。与无时间展宽的 125 MHz 上限相比，当前光学核使 SNR 乘以

$$
\sqrt{\frac{110.910}{125}}=0.942,
$$

即降低约 $5.8\%$；达到相同 SNR 所需的曝光时间增加为

$$
\frac{125}{110.910}=1.127,
$$

即约增加 $12.7\%$。这说明当前 $0.549\ \mathrm{ns}$ 单镜光学展宽有影响，但对匹配相关器不是灾难性的。若直接使用原始零延迟协方差，则实测长 SPE 尾部使带宽进一步降到下一节的 $B_{\rm raw}=14.02\ \mathrm{MHz}$，损失会大得多。

### 6.5 为什么 170 ns SPE 没把最优带宽降到几 MHz

令

$$
H_{\rm SPE}(f)=\int h_{\rm SPE}(t)e^{-2\pi ift}\,dt
$$

为实测 SPE 模板的傅里叶变换。在线性、稳定且散粒噪声占主导时，SPE 会同时滤波 HBT 信号和光子散粒噪声。已知 $H_{\rm SPE}$ 后，匹配相关器可以在有响应的频率上做最优加权，因此它不会仅因“波形尾巴长 170 ns”就把信息带宽设成 $1/170\ \mathrm{ns}$。

但是，若直接取两条原始电压波形的零延迟协方差而不做匹配加权，长 SPE 尾部会显著增加相邻样点相关性。本文为这个原始估计器定义

$$
B_{\rm raw}\equiv
\frac{
\left[\int_0^{f_{\rm Nyq}}
|H_{\rm SPE}(f)|^2|H_{\rm opt}(f)|^2\,df\right]^2
}{
\int_0^{f_{\rm Nyq}}|H_{\rm SPE}(f)|^4\,df
}.
$$

当前实测 SPE、4 ns 采样和光学核给出

$$
B_{\rm raw}=14.0217\ \mathrm{MHz}.
$$

因此同样积分时间下，原始零延迟估计器与匹配估计器的 SNR 比约为

$$
\sqrt{\frac{B_{\rm raw}}{B_{\rm eff}}}=0.356,
$$

即误差约大 2.81 倍，需要约 $B_{\rm eff}/B_{\rm raw}=7.91$ 倍曝光才能达到相同精度。

如果将来加入实测电子噪声，$H_{\rm SPE}$ 会通过 $\eta_{\rm shot}(f)$ 影响最优带宽，因为高频处信号和散粒噪声虽同时变小，附加电子噪声却不一定同步变小。那时必须用实测噪声功率谱重新计算 $B_{\rm eff}$，不能继续沿用 110.91 MHz。

---

## 7. 两台望远镜的完整短波形怎样生成和相关

### 7.1 逐光电子生成，而不是逐电磁波振荡生成

光学载波频率约为 $10^{15}\ \mathrm{Hz}$，没有必要也无法用 4 ns ADC 直接采样电场振荡。探测器输出由光电子点过程决定，因此短波形模拟按以下顺序进行：

1. 分别按 $r_\star$、$r_{\rm NSB}$ 和 $r_{\rm dark}$ 抽取 Poisson 光电子到达时间；
2. 按第 5.3 节的 $R_{{\rm pair},ij}$ 抽取额外的 HBT 相关对；
3. 对两镜各自的光电子抽取 $K_{\rm opt}(t)$ 延迟；
4. 随机分配 SiPM 微单元，应用 $\tau_{\rm rec}=10\ \mathrm{ns}$ 恢复；
5. 为每个光电子抽取实测电荷因子 $a_k$；
6. 将每个击中与实测 $h_{\rm SPE}(t)$ 卷积；
7. 加入已配置的电子噪声和 ADC 量化；当前两项均为 0；
8. 在 $\Delta t_{\rm samp}=4\ \mathrm{ns}$ 网格上得到两路电压波形 $x_i[n]$ 和 $x_j[n]$。

相关对构造是热光二阶统计量的 Monte Carlo 表示，不表示“同一颗光子被两台望远镜同时探测”。它保证额外协方差的期望值等于第 5.3 节的物理相关面积。

### 7.2 先校正几何时延，再相关

第 3.2 节已经定义几何时延 $\tau_g=w_m/c$。程序将其中一路移位到共同时间坐标后，计算去均值离散相关：

$$
\widehat C_{ij}[k]=
\frac{1}{N_k}
\sum_n
\bigl(x_i[n]-\bar x_i\bigr)
\bigl(x_j[n+k]-\bar x_j\bigr),
$$

其中 $k$ 是离散延迟索引，$N_k$ 是该延迟下重叠样点数，$\bar x_i$ 是记录内平均电压。归一化后为

$$
\widehat\rho_{ij}[k]=
\frac{\widehat C_{ij}[k]}
{\widehat\sigma_i\widehat\sigma_j},
$$

其中 $\widehat\sigma_i$ 和 $\widehat\sigma_j$ 是两路波形的样本标准差。

不能先把整段波形积分成单个总电荷

$$
Q_i=\int x_i(t)\,dt
$$

再只相关 $Q_i$ 和 $Q_j$，因为那会丢失纳秒尺度的时间结构。可以等价地相关细时间格 PE 序列或经过已知线性响应的电压波形，但必须保留相对于 $B_{\rm eff}$ 足够细的时间信息。

![短波形与 PE 恢复](docs/sii_workflow_figures/04_short_waveform_and_recovery.png)

**图 4。** 左：实测 SPE 模板。中：同一次 2 µs Monte Carlo 得到的 A、B 两台望远镜电压波形。右：直接由中图这两条波形计算的归一化延迟相关。右图确实与中图一一对应，但其起伏由有限记录噪声主导；零延迟附近没有可识别的物理 HBT 峰。

![两台望远镜的相关峰](docs/sii_workflow_figures/05_two_telescope_g2_peak.png)

**图 5。** 由 $\tau_{\rm eff}$、该基线的 $P$、几何时延、光学时间核和采样响应计算的**期望电子相关峰形**。它不是图 4 那一次随机波形的测量相关；图中没有叠加 2 µs 记录的随机噪声，所以能看见均值只有几 $\times10^{-6}$ 的理论峰。

### 7.3 图 4 和图 5 的关系

两图使用相同的源、基线、光子率和仪器响应，但表示不同统计对象：

- 图 4 是一次有限样本 realization；中图两条电压波形直接产生右图的噪声主导相关曲线。
- 图 5 是把同一随机实验无限重复后相关曲线的期望值，也可由已知传递函数直接算出。

因此“图 5 有理论峰”和“图 4 的 2 µs 实测相关看不见峰”可以同时成立。类比来说，单次抛硬币可能严重偏离 50%，但概率模型的期望仍是 50%。只有累积足够多独立时间样本，图 4 类型的测量平均值才逐渐趋近图 5。

更准确地说，每个 2 µs 片段有超额相关对的概率虽只有约 $0.1364\%$，但 20 min 包含

$$
N_{\rm rec}=\frac{1200\ \mathrm{s}}{2\ \mu\mathrm{s}}=6.0\times10^8
$$

个这样的片段。20 min 内的超额相关对期望总数为

$$
\langle N_{\rm HBT}\rangle
=R_{\rm pair}T
=682.15\times1200
\simeq8.19\times10^5.
$$

没有超额相关对的短片段占绝大多数，但它们的去均值随机相关有正有负，平均时主要互相抵消；少数 HBT 超额项总在校正后的同一延迟处贡献正均值，于是累计平均逐渐留下图 5 的峰形。峰的**期望存在**不等于一次有限曝光中峰已被显著测出，显著度仍由第 8.2 节的 SNR 决定。

### 7.4 为什么 2 µs 看不到真实 HBT 峰

一次记录中期望物理相关对数为

$$
\mu_{\rm pair}=R_{\rm pair}T_{\rm rec}.
$$

默认过中天示例代入 $R_{\rm pair}=682.15\ \mathrm{s^{-1}}$ 和 $T_{\rm rec}=2\ \mu\mathrm{s}$，得到

$$
\mu_{\rm pair}=1.364\times10^{-3}.
$$

若相关对计数服从 Poisson 分布，完全没有相关对的概率为

$$
\Pr(N_{\rm pair}=0)=e^{-\mu_{\rm pair}}=99.864\%.
$$

因此一条真实强度的 2 µs 波形几乎必然没有额外相关对，不能据此判断天体是否可见。程序使用两个互补检查：

- 零增强：测量纯噪声底和估计器方差；
- 仅增强 HBT 相关对：不改变恒星独立光子和 NSB，用已知倍数放大相关项，验证测量峰与预测值一致。

增强 10000 倍时，实测原始零延迟相关为 $0.02394\pm0.00629$，预测为 $0.02420$，说明从相关对、光学核、SPE 到相关器的实现闭合。

---

## 8. 从短记录推广到长曝光统计量

### 8.1 为什么可以不保存整夜波形

若系统在一个统计时间段内近似平稳、响应线性、时延已校正，而且 $K_{\rm opt}$、$h_{\rm SPE}$ 和噪声谱已经由短波形链路标定，那么关于 $P$ 的主要信息可压缩为相关估计值、有效带宽和方差。250 MS/s 连续记录 6 小时，每台望远镜有 $5.4\times10^{12}$ 个样点；32 台超过 $10^{14}$ 个样点。逐点保存和卷积不会增加新的天体统计信息，只会极大增加计算量。

因此当前设计是：短波形负责验证仪器和估计器，长曝光使用由同一物理模型推导的充分统计量。两者不是两套互相替代的物理模型。

### 8.2 单基线长曝光 SNR

令 $T_{\rm seg}$ 为一条基线在一个时间段内的有效积分时长。先用计数窗口的图像推导 SNR，再把真实响应等效为第 6 节的 $B_{\rm eff}$。

在本文的单边带宽约定下，定义等效独立时间窗口

$$
\Delta t_{\rm eff}\equiv\frac{2}{B_{\rm eff}}.
$$

$\Delta t_{\rm eff}$ 不是 4 ns ADC 采样间隔，也不是 SPE 长度；它是把匹配相关器的信息量换算成理想矩形相关器后，一个统计自由度所占的等效时间。令 $N_i$ 和 $N_j$ 是一个等效窗口内两台望远镜记录的总光电子数，其 Poisson 均值分别为

$$
\mu_i=r_{{\rm tot},i}\Delta t_{\rm eff},
\qquad
\mu_j=r_{{\rm tot},j}\Delta t_{\rm eff}.
$$

相关器实际使用去均值计数涨落

$$
\delta N_i=N_i-\mu_i,
\qquad
\delta N_j=N_j-\mu_j,
$$

而不是直接把原始 $N_iN_j$ 当成天体信号。一个窗口内，来自恒星的 HBT 超额相关使 $\delta N_i\delta N_j$ 出现非零期望值；由第 5.3 节的相关对率得到：

$$
\mu_{{\rm HBT},ij}
=r_{\star,i}r_{\star,j}\tau_{\rm eff}P_{ij}\Delta t_{\rm eff},
$$

其中 $P_{ij}=P(u_{ij},v_{ij})$ 是该基线和时刻的理论功率可见度。若先忽略很弱的 HBT 协方差，两路无相关 Poisson 涨落满足

$$
\mathrm{E}[\delta N_i\delta N_j]=0.
$$

但是它在一次有限测量中不会恰好等于 0。由于独立 Poisson 变量有

$$
\mathrm{Var}(N_i)=\mu_i,
\qquad
\mathrm{Var}(N_j)=\mu_j,
$$

乘积 $X_{ij}=\delta N_i\delta N_j$ 的方差为

$$
\begin{aligned}
\mathrm{Var}(X_{ij})
&=\mathrm{E}[\delta N_i^2]\,\mathrm{E}[\delta N_j^2]\\
&=\mu_i\mu_j\\
&=r_{{\rm tot},i}r_{{\rm tot},j}\Delta t_{\rm eff}^2.
\end{aligned}
$$

所以无相关随机符合留下的标准差为

$$
\sigma_{{\rm acc},ij}
=\sqrt{\mathrm{Var}(X_{ij})}
=\sqrt{r_{{\rm tot},i}r_{{\rm tot},j}}\,\Delta t_{\rm eff}.
$$

这个表达来自**去均值的两路 Poisson 涨落乘积的方差**。如果直接研究未去均值的原始乘积 $N_iN_j$，还会出现额外高阶项；当前 HBT 相关器先减均值，因此不使用那个表达。在 HBT 信号远小于总散粒噪声时，忽略信号本身对方差的微小修正是合理近似。

因此单个等效窗口的信噪比是

$$
\frac{\mu_{{\rm HBT},ij}}{\sigma_{{\rm acc},ij}}
=\frac{r_{\star,i}r_{\star,j}\tau_{\rm eff}}
{\sqrt{r_{{\rm tot},i}r_{{\rm tot},j}}}P_{ij}.
$$

$\Delta t_{\rm eff}$ 在单窗口信号和噪声中相消。一个 $T_{\rm seg}$ 时间段包含

$$
M_{\rm eff}=\frac{T_{\rm seg}}{\Delta t_{\rm eff}}
=\frac{B_{\rm eff}T_{\rm seg}}{2}
$$

个等效独立时间自由度，其中 $M_{\rm eff}$ 是有效独立样本数。独立样本的 SNR 按 $\sqrt{M_{\rm eff}}$ 合并；再除以单 PE 电荷涨落惩罚 $F_{\rm EN}$，便得到

$$
\boxed{
\mathrm{SNR}_{ij}=
\frac{r_{\star,i}r_{\star,j}\tau_{\rm eff}}
{\sqrt{r_{{\rm tot},i}r_{{\rm tot},j}}}
\,P_{ij}
\sqrt{\frac{B_{\rm eff}T_{\rm seg}}{2}}
\frac{1}{F_{\rm EN}}
}.
$$

这不是突然引入的经验公式，而是“每个窗口的 HBT 超额相关对 ÷ 随机符合涨落”，再乘“独立窗口数的平方根”。推导使用弱相关、Poisson 散粒噪声主导、时间段内平稳、匹配响应已知的近似；附录 B 再用量纲和极限情况检查该公式。

公式中的每个量都已在前文定义：

- $r_{\star}$ 和 $r_{\rm tot}$ 来自第 4 节；NSB 只增加 $r_{\rm tot}$，所以会降低 SNR；
- $\tau_{\rm eff}$ 来自第 5 节，包含光学频带、偏振和谱形；
- $P_{ij}$ 来自第 2 节并在第 3 节的 UV 点上取值；
- $B_{\rm eff}$ 来自第 6.4 节，包含光学时间展宽、采样和已知电子噪声权重；
- $T_{\rm seg}$ 是本次统计积分时间；
- $F_{\rm EN}$ 来自第 6.2 节的实测单 PE 电荷涨落；
- 因子 $1/2$ 来自本文采用的单边带宽和实值相关器方差约定。

这里没有重复计算：$\tau_{\rm eff}$ 决定超窄光学相干峰的总面积，$B_{\rm eff}$ 决定纳秒电子系统能以多少独立时间模态估计这块面积。

为方便表达误差，定义单位功率可见度的 SNR：

$$
S_{ij}^{(1)}\equiv
\frac{r_{\star,i}r_{\star,j}\tau_{\rm eff}}
{\sqrt{r_{{\rm tot},i}r_{{\rm tot},j}}}
\sqrt{\frac{B_{\rm eff}T_{\rm seg}}{2}}
\frac{1}{F_{\rm EN}}.
$$

于是

$$
\mathrm{SNR}_{ij}=P_{ij}S_{ij}^{(1)},
\qquad
\sigma_{P,ij}=\frac{1}{S_{ij}^{(1)}},
$$

其中 $\sigma_{P,ij}$ 是一次 $P_{ij}$ 测量的标准差。注意 $\sigma_P$ 与真实 $P$ 无关，而实际信号显著度随 $P$ 线性下降。

默认相同望远镜、$m_{\rm AB}=2$、参考 NSB、$T_{\rm seg}=1200\ \mathrm{s}$ 时：

$$
S^{(1)}=4.25818,
\qquad
\sigma_P=0.234842.
$$

前述示例基线 $P=0.149476$ 的实际单段显著度只有

$$
\mathrm{SNR}=0.149476\times4.25818=0.636.
$$

这解释了为什么单个 20 分钟 UV 点很噪，而组合 496 条基线和多个时间段后仍可能恢复整体形态。

$\mathrm{SNR}=1$ 只表示信号期望值等于一次测量的 $1\sigma$ 标准差，不是高显著性检出。若误差近似高斯，单个预先指定位置的双侧 $1\sigma$ 偏离在纯噪声下并不罕见；单点通常至少需要约 $3\sigma$ 才能称为有证据，严格发现往往要求更高并考虑多重检验。低于 $1\sigma$ 的点仍可作为带大误差的无偏数据参与联合拟合，但不能单独声称看见 HBT 峰。

因此，如果真的连续模拟 20 min 完整波形并使用与长曝光公式相同的匹配相关器，统计分布应与这里生成的 $\widehat P$ 和 $\sigma_P$ 一致，而不是每个 2 µs 片段各自先判断“有峰/无峰”再投票。对当前 $P=0.149476$ 基线，20 min 的匹配 SNR 仍只有 0.636，单次相关曲线通常不会出现可靠可辨认的峰；对 $P\simeq1$ 的短基线，同样 20 min 的匹配 SNR 为 4.258，峰在统计上应明显得多。逐样点生成 20 min 波形不会改变结果，只会把解析完成的每台望远镜约 $3\times10^{11}$ 个 ADC 样点显式算一遍。

### 8.3 短波形和长曝光的定量闭合

同一公式用于 $T_{\rm rec}=2\ \mu\mathrm{s}$ 时，单位 $P$ 的匹配 SNR 为

$$
S^{(1)}(2\ \mu\mathrm{s})=1.7384\times10^{-4},
$$

远小于 1，与第 7.4 节“几乎总是零相关对”的结论一致。将记录时长增加时，独立模拟测得的相关估计器标准差满足

$$
\sigma\propto T^{-1/2},
$$

而长曝光公式正是同一缩放。图 6 同时检查了相关对均值、原始波形带宽惩罚和 $T^{-1/2}$ 噪声缩放。

![短波形与长曝光一致性](docs/sii_workflow_figures/05a_short_long_consistency.png)

**图 6。** 短波形 Monte Carlo 和长曝光解析统计量的闭合。匹配相关器使用 $B_{\rm eff}=110.91\ \mathrm{MHz}$；直接原始零延迟相关使用 $B_{\rm raw}=14.02\ \mathrm{MHz}$。

### 8.4 未校准时间误差怎样进入长曝光

已知且稳定的 DC 光学时间核已经包含在 $B_{\rm eff}$ 中，不应再额外乘一个“光学效率”。只有未被模型校正的相对时间误差才需要单独衰减。

令 $\Delta\tau_{ij}$ 是基线 $ij$ 在校正 $\tau_g$ 后残留的相对延迟，定义归一化残余延迟响应

$$
A_{ij}(\Delta\tau_{ij})=
\frac{
\int_0^{f_{\rm Nyq}}
|H_{\rm opt}(f)|^2
\cos(2\pi f\Delta\tau_{ij})\,df
}{
\int_0^{f_{\rm Nyq}}|H_{\rm opt}(f)|^2\,df
}.
$$

若残余延迟未校正，则期望测量值变为 $A_{ij}P_{ij}$。当前 main 未提供时钟漂移或残余时延 RMS，故 $\Delta\tau_{ij}=0$、$A_{ij}=1$。这只是默认关闭，不代表真实系统一定没有该误差。

---

## 9. 组合全部基线和时间段，生成带噪稀疏 UV 数据

### 9.1 每个测量点

对每个有效时间段和每条基线，程序依次：

1. 用第 3 节公式计算 $(u_{ij},v_{ij},w_{ij})$ 和 $\tau_g$；
2. 用第 2 节解析模型计算真值 $P_{ij}$；
3. 用第 4 节光子率和第 8.2 节公式得到 $\sigma_{P,ij}$；
4. 抽取统计噪声，并加入启用的标定误差；
5. 保存 $(u,v)$、测得的 $\widehat P$、$\sigma_P$、望远镜编号和时刻。

统计模型可写成

$$
\widehat P_{ij,t}=A_{ij,t}G_{ij,t}P(u_{ij,t},v_{ij,t})
+z_{ij,t}+\epsilon_{ij,t},
$$

其中：

- $\widehat P_{ij,t}$ 是基线 $ij$ 在时间段 $t$ 的测量值；
- $A_{ij,t}$ 是第 8.4 节的残余时间响应；
- $G_{ij,t}$ 是望远镜增益、透明度等乘性标定因子的组合；
- $z_{ij,t}$ 是基线相关零点偏差；
- $\epsilon_{ij,t}$ 是均值为 0、标准差为 $\sigma_{P,ij,t}$ 的统计噪声。

当前缺少实测值的 $G$ 波动、$z$、时钟漂移、额外电子噪声、串扰、后脉冲和暗计数都设为 0，但接口保留。每台望远镜自己的 Poisson 光子涨落不是这些“可选系统误差”，它始终通过 $r_{\rm tot}$ 和显式短波形包含。

程序不把负的 $\widehat P$ 强行裁剪为 0。低 SNR 无偏相关估计出现负值是正常现象；提前裁剪会引入正偏差。重建阶段再通过非负天空图约束处理物理性。

### 9.2 时间分段和点数

默认长曝光单段使用

$$
T_{\rm seg}=1200\ \mathrm{s}=20\ \mathrm{min}.
$$

6 小时共有 18 个时间段，理论最多产生

$$
N_{\rm meas}=496\times18=8928
$$

个基线时间测量；$N_{\rm meas}$ 表示测量数。高度角筛选会删除不可用时段。每个时段不是只模拟一个“平均 UV 点”，而是每条基线都有自己的坐标、$P$ 真值和噪声。

![带噪声的稀疏 UV 测量](docs/sii_workflow_figures/06_simulated_uv_measurements.png)

**图 7。** 默认参考场景的 8928 个稀疏 $\widehat P(u,v)$。点的位置由几何决定，颜色由源真值、仪器响应和随机噪声共同决定。绘图色标限制在 0 到 1，超出范围的噪声值只在显示时颜色饱和，CSV 中的数据没有裁剪。右图显示测量值相对真值仍有明显散布。正负基线使用相同测量值和显著度。

---

## 10. 从无相位功率可见度重建天空图

### 10.1 离散正向模型

令 $x_p$ 是第 $p$ 个天空像素的非负亮度，像素角坐标为 $(l_p,m_p)$。程序施加

$$
x_p\ge0,\qquad \sum_p x_p=1,
$$

即非负性和总流量归一化。对第 $k$ 个 UV 点 $(u_k,v_k)$，预测复可见度为

$$
V_k(\mathbf x)=
\sum_p x_p\exp[-2\pi i(u_kl_p+v_km_p)],
$$

预测功率可见度为

$$
P_k(\mathbf x)=|V_k(\mathbf x)|^2.
$$

$\mathbf x$ 表示由所有像素值组成的图像向量。

### 10.2 UV 合并和逆方差加权

直接优化前，临近 UV 测量按网格合并。若某一格内包含测量 $k$，其权重定义为

$$
w_k=\frac{1}{\sigma_{P,k}^2}.
$$

合并值和误差为

$$
\overline P=
\frac{\sum_k w_k\widehat P_k}{\sum_k w_k},
\qquad
\sigma_{\overline P}=\left(\sum_k w_k\right)^{-1/2}.
$$

这一步减少重复点和优化规模，同时保留不同测量精度。

### 10.3 优化目标

重建最小化加权数据失配和弱正则项：

$$
\mathcal L(\mathbf x)=
\frac{1}{2}\sum_k
\left[
\frac{P_k(\mathbf x)-\widehat P_k}{\sigma_{P,k}}
\right]^2
+\lambda_{\rm TV}\,\mathrm{TV}(\mathbf x)
+\lambda_2\|\mathbf x\|_2^2.
$$

其中 $\mathrm{TV}(\mathbf x)$ 是总变分，抑制噪声产生的像素级振荡；$\lambda_{\rm TV}$ 是总变分权重；$\lambda_2$ 是二范数正则权重。每次迭代后投影回 $x_p\ge0$ 且 $\sum x_p=1$ 的集合。

因为目标函数是非凸的，程序使用多个初值并比较正向残差。图 8 对应的默认 6 小时 realization 中，加权正向 RMSE 为 0.0831；共同分辨率卷积后与真值的相关系数为 0.9647。恢复分离角为 0.2025 mas，真值为 0.2000 mas；恢复位置角为 $39.8^\circ$，真值为 $35.0^\circ$；恢复峰值比为 0.799，真值流量比为 0.55。后两项的偏差说明“能看出是双星”不等于位置角和测光都已准确。

### 10.4 为什么单点 SNR 为 0.636 仍能重建

先区分模拟真值和重建输入。前向模拟必须用源模型给出的

$$
P_k^{\rm true}=P(u_k,v_k)
$$

规定随机实验的均值，然后生成一次带噪测量

$$
\widehat P_k=P_k^{\rm true}+\epsilon_k,
\qquad
\epsilon_k\sim\mathcal N(0,\sigma_{P,k}^2),
$$

其中 $k$ 编号不同基线和时间段。重建器接收的是 $(u_k,v_k,\widehat P_k,\sigma_{P,k})$，不是 `visibility2_true`；它也不读取双星分离角、位置角、流量比或真傅里叶相位。程序额外加入 $P(0,0)=1$ 作为总流量归一化，并使用非负性、有限支撑和平滑先验。真值只在重建结束后用于计算误差指标，以及消除强度干涉无法确定的平移和 180° 镜像后画对比图；这种事后配准不会改变重建优化得到的图像内容。

第 8.2 节的 0.636 只对应一个 $P=0.149476$ 的特定基线、一个 20 min 时间段，不代表所有 UV 点。当前数据中 $P$ 从接近 0 延伸到 0.985；当 $P\simeq1$ 时，同一时间段的 SNR 是 4.258。8928 个测量的中位单点 SNR 为 1.381，约 $60.1\%$ 的点大于 1，约 $13.8\%$ 的点大于 3。

重建也不是从一个 0.636 的点作图，而是同时拟合所有测量。以 $120\ \mathrm{M}\lambda$ 网格合并后，8928 个测量形成 1021 个实际 UV 格；再加入由总流量归一化给出的零基线 $P(0,0)=1$，优化器共使用 1022 个约束点。一个 UV 格的测量重复数中位数为 6，因此其误差中位数从原始单点的 0.2348 降到约 0.0959。源的信号还具有跨 UV 平面连续、成体系的条纹结构，非负性、有限支撑、总流量归一化和平滑正则会排除大量纯噪声图像。

所以“单点不显著”和“整体模型显著”并不矛盾，类似很多低 SNR 像素共同检出一个弱而连续的整体结构。但每次 $\widehat P_k$ 都会偏离真值，因此每次重建也会变化；误差变大时，双星可能被重建成错误分离度、错误流量比、单峰或噪声伪双峰。Notebook 不是只展示一个幸运 realization：它还对每种星等和曝光生成多次独立带噪数据，并用单星空白样本估计假双峰率。只有双星真阳性率高、单星假阳性率低时，才把该情景列为可靠检出。

不过当前重建仍偏乐观，原因包括：源很亮（$m_{\rm AB}=2$）；缺失的增益漂移、时钟误差和附加电子噪声均设为 0；匹配时间响应假设完全已知；长曝光统计误差主要按基线独立高斯量抽取，尚未使用完整的跨基线协方差矩阵；重建还使用非负、有限视场和弱平滑先验。因此图 8 是在当前假设下的可行性结果，不是对真实 LACT 六小时必然成像质量的保证。

### 10.5 强度干涉固有歧义

仅有 $P=|V|^2$ 时，至少存在以下信息损失：

- 整体平移不改变 $P$；
- 180° 反演可给出相同 $P$；
- 稀疏 UV 覆盖和噪声可产生多个局部极小值；
- 不对称源的相位信息不由单个二阶相关数据集直接提供。

因此评价重建时，程序允许先做平移和 180° 镜像配准，再在共同分辨率下比较；这不是人为改善观测，而是消除不可辨识自由度。

![无相位图像重建](docs/sii_workflow_figures/07_phaseless_reconstruction.png)

**图 8。** 真值、带噪稀疏测量和无相位重建。图像清楚并不等于每个单独 UV 点都有高 SNR；重建使用全部基线、全部时段和先验约束共同求解。

---

## 11. 观测时间、星等极限和角分辨率

### 11.1 时间缩放

从第 8.2 节直接得到

$$
\mathrm{SNR}\propto\sqrt{T_{\rm seg}},
\qquad
\sigma_P\propto T_{\rm seg}^{-1/2}.
$$

在 UV 几何已接近饱和后，继续重复观测仍会降低误差，但收益是平方根而非线性。例如相同 UV 网格在 5 夜共 30 h 和 10 夜共 60 h 后，重复测量逆方差合并的等效误差分别约为 0.105 和 0.0743；6 h 单夜中一个 20 min 原始测量的误差仍为 0.2348。这里降低的是重复采样合并后的误差，不是把每个原始 20 min 点改写成更高精度。

### 11.2 角分辨率

令 $B_{\max}$ 为最长投影基线。常用条纹尺度为

$$
\theta_{\rm fringe}\simeq\frac{\lambda_0}{B_{\max}}.
$$

对于均匀圆盘第一零点，尺度为

$$
\theta_{\rm null}\simeq1.22\frac{\lambda_0}{B_{\max}}.
$$

当前 32 台 `layout_0803_reco32` 的物理最长基线约为 $1219.21\ \mathrm{m}$。在 $400\ \mathrm{nm}$ 下，名义条纹尺度约为 $0.0677\ \mathrm{mas}$，均匀圆盘第一零点尺度约为 $0.0826\ \mathrm{mas}$。这些是几何上限，不等于在有限 SNR、有限方向覆盖和无相位重建下必然可分辨的最小双星间距。

最短基线 $B_{\min}\simeq121.01\ \mathrm{m}$ 对应约 $\lambda_0/B_{\min}=0.682\ \mathrm{mas}$ 的最大敏感结构尺度；缺少更短基线会削弱大尺度形态恢复。

### 11.3 星等极限

星等极限由第 8.2 节公式反解，而不是只根据源图是否“看起来像双星”判断。对于未分辨源的 $P=1$、参考 NSB、单光谱通道，当前解析估计为：

| 总积分时间 | $\mathrm{SNR}=5$ 的 AB 星等 | $\mathrm{SNR}=10$ 的 AB 星等 |
|---:|---:|---:|
| 1 h | 2.327 | 1.733 |
| 10 h | 3.210 | 2.693 |
| 50 h | 3.763 | 3.289 |

这些数字是“单基线、单位 $P$、假设 2 nm 恒星通道、main 宽带 NSB 压力值和当前电子带宽”条件下的检测门槛，不是最终硬件星等极限，也不是复杂双星可靠成像极限。真实基线上 $P<1$ 时，显著度应再乘 $P$。多光谱通道若彼此统计独立且每个通道保持相同窄带性能，总 SNR 理想情况下按 $\sqrt{N_{\rm ch}}$ 增长，其中 $N_{\rm ch}$ 是独立光谱通道数；工程损耗和通道间相关会降低该收益。最终灵敏度必须在恒星和 NSB 共用同一个实测 $T(\lambda)$ 后重新计算。

---

## 12. 当前参数来源、已包含项和待标定接口

| 物理项 | 当前状态 | 当前处理 |
|---|---|---|
| 32 台望远镜坐标 | 已有 | 直接读取 `layout_0803` run-card 导出 |
| 遮挡后面积 $A_{\rm col}$ | 已有 | $24.57686\ \mathrm{m^2}$ |
| 单像素净面积效率乘积 | 光追标定 | $4.64128\ \mathrm{m^2}$，含中心像素光追和上游传输因子 |
| 镜面、宽带相机窗口、PDE、集光器、PSF | 已有 | 进入 main 的 400 nm 单像素光追响应 |
| SII 光谱通道 $\lambda_0,\Delta\lambda$ | 场景假设 | 当前取 400 nm、2 nm；不是已有 LACT 窄带片标定 |
| 谱形因子 $s_{\rm spec}$ | 场景假设 | 当前取 0.842；需由最终 SII 总光谱响应重新积分 |
| 光学到达时间核 $K_{\rm opt}$ | 光追标定 | 100 万平行光子，$\sigma_{\rm opt}=0.548827\ \mathrm{ns}$ |
| 实测 SPE $h_{\rm SPE}$ | 已有 | 537 个干净脉冲模板，用于完整短波形 |
| SPE 电荷涨落 $F_{\rm EN}$ | 已有 | $F_{\rm EN}=1.016142$，进入短波形和长曝光 SNR |
| 暗夜 NSB | main 宽带模型估计 | $70.5275\ \mathrm{MHz/pixel}$；当前作为保守压力值，尚未与假设的 2 nm SII 通道统一 |
| 微单元数 | 已有 | 270336 个/像素 |
| 微单元恢复时间 $\tau_{\rm rec}$ | 暂定 | 10 ns，可替换 |
| ADC 采样 | 已有 | 4 ns，即 250 MS/s、Nyquist 125 MHz |
| 有效相关带宽 $B_{\rm eff}$ | 由响应推导 | 当前 $110.910\ \mathrm{MHz}$，不是自由输入常数 |
| ADC 位数和满量程 | 缺失 | 接口保留，0 表示关闭量化 |
| 附加电子噪声谱 | 缺失 | 接口保留，当前为 0；加入后需重算 $\eta_{\rm shot}(f)$ 和 $B_{\rm eff}$ |
| 串扰、后脉冲、暗计数 | 缺失 | 接口保留，当前均为 0 |
| 时钟漂移和残余时延 | 缺失 | 接口保留，当前为 0 |
| 望远镜增益、透明度、NSB 稳定性 | 缺失 | 接口保留，当前随机 RMS 均为 0 |
| 随天顶角、方位角变化的光学核 | 缺失 | 当前只用轴上核；未来可按方向替换 $K_{\rm opt}$ |
| 窄带滤光片角响应 | 缺失 | 路径接口保留，当前只使用轴上频带 |

因此，当前时间响应链和短—长曝光统计闭合是物理一致的，但光谱部分仍是“2 nm 恒星通道 + main 宽带 NSB 压力值”的保守混合场景，不是最终自洽通带。缺失系统误差设 0 可避免虚构参数，但重建质量应理解为这些误差得到控制时的可行性结果。形成硬件灵敏度结论前，必须用同一个总响应 $T(\lambda)$ 同时计算 $r_\star$、$r_{\rm NSB}$、$s_{\rm spec}$ 和 $\tau_{\rm eff}$。

---

## 13. 如何运行自己的参数

### 13.1 一次完整前向模拟

`python/sii_unified.py` 将仪器、源、阵列和观测计划分开。典型调用关系是：

```python
base_instrument = Instrument.from_repository(repo_root)
source = BinarySource(...)
observation = Observation(...)

# from_repository 先给出 ADC Nyquist 上限；再用实际响应计算 B_eff。
instrument = with_matched_effective_bandwidth(
    base_instrument,
    source_ab_magnitude=source.ab_magnitude,
)

result = run_sii_pipeline(
    layout=positions,
    source=source,
    instrument=instrument,
    observation=observation,
    seed=20260903,
)
```

`Instrument.from_repository(...)` 会重新读取 main 当前使用的镜面/宽带窗口/PDE 配置、单像素光追响应、SPE、SiPM 和采样参数，但其中的 `electronics_bandwidth_hz` 首先只是 $f_{\rm Nyq}$。公共函数 `with_matched_effective_bandwidth(...)` 再执行第 6.4 节的积分，同时标记光学时间核已经进入 $B_{\rm eff}$，从而避免长曝光里重复乘一次光学时间效率。不能跳过该步直接把 125 MHz 当成最终 $B_{\rm eff}$。没有标定值的接口仍保持 0，除非用户显式赋值。

需要额外注意：main 参数变化可以自动传播，但当前假设的 2 nm SII 光谱通道并不来自 main。若更换 SII 滤光片，必须同时更新恒星光谱积分、NSB 光谱积分和 $\tau_{\rm eff}$，不能只改 `optical_width_nm` 一个数字。

### 13.2 重建独立运行

前向模拟输出带有 $(u,v,\widehat P,\sigma_P)$ 的表格。`python/sii_reconstruction.py` 只依赖这张表和成像网格，因此可以：

- 对同一批模拟数据比较不同重建超参数；
- 对不同源、观测时间、NSB 或电子学场景重复运行统一前向流程；
- 将来把模拟表替换成格式相同的真实相关测量，而不重写重建器。

notebook 把“源 → 理论 UV → 阵列 UVW → 仪器标定 → 短波形闭合 → 长曝光数据 → 重建 → 星等/分辨率扫描”完整串联，并输出所有表格和图。因此它既是论文结果的可复现入口，也是修改参数后的工作模板。

---

## 14. 一句话总结

当前程序先用源模型计算连续理论 $P(u,v)$，再用 32 台望远镜和地球自转决定采样位置；单镜光子率与光学相干面积给出 HBT 信号强度，光追时间核、实测 SPE、采样和噪声谱决定估计器的有效带宽；微秒完整波形验证这套响应，小时级充分统计量以同一公式生成每个带误差的 UV 点，最后用非负、归一化、正则化的无相位优化重建天空图。短波形回答“真实仪器链是否实现了预期相关统计”，长曝光回答“在实际观测时间内这些相关统计能测到多准”；两者必须闭合，但用途不同。

---

## 附录 A：为什么相干面积与 $1/\Delta\nu$ 成正比

### A.1 从归一化频谱到一阶相干函数

令 $\psi(\nu)$ 表示探测到的归一化光子频谱形状，$\nu$ 是光频率，并规定

$$
\int_{-\infty}^{\infty}\psi(\nu)\,d\nu=1.
$$

$\psi$ 的单位是 $\mathrm{Hz^{-1}}$。Wiener–Khintchine 关系给出归一化一阶时间相干函数

$$
g^{(1)}(\tau)=
\int_{-\infty}^{\infty}
\psi(\nu)e^{2\pi i\nu\tau}\,d\nu.
$$

在本文使用的傅里叶变换约定下，Parseval 定理给出

$$
\int_{-\infty}^{\infty}|g^{(1)}(\tau)|^2\,d\tau
=
\int_{-\infty}^{\infty}|\psi(\nu)|^2\,d\nu.
$$

左边单位是 s；右边的 $|\psi|^2d\nu$ 单位是 $\mathrm{Hz^{-1}}=\mathrm{s}$，量纲一致。

### A.2 矩形频带的直接计算

对中心频率 $\nu_0$、宽度 $\Delta\nu$ 的理想矩形频带，归一化频谱为

$$
\psi(\nu)=
\begin{cases}
1/\Delta\nu,
& |\nu-\nu_0|\le\Delta\nu/2,\\
0,&\text{其他频率}.
\end{cases}
$$

直接代入 Parseval 关系：

$$
\int |g^{(1)}(\tau)|^2d\tau
=\int |\psi(\nu)|^2d\nu
=\left(\frac{1}{\Delta\nu}\right)^2\Delta\nu
=\frac{1}{\Delta\nu}.
$$

也可以先做傅里叶变换：

$$
g^{(1)}(\tau)=
e^{2\pi i\nu_0\tau}\,
\operatorname{sinc}(\Delta\nu\tau),
$$

其中 $\operatorname{sinc}(x)=\sin(\pi x)/(\pi x)$。相位因子取绝对值后为 1，而 $\operatorname{sinc}^2$ 峰的宽度与 $1/\Delta\nu$ 成正比，所以积分仍为 $1/\Delta\nu$。

### A.3 非矩形频谱和当前 $s_{\rm spec}$

一般频谱不必是矩形。本文定义无量纲谱形因子

$$
s_{\rm spec}\equiv
\Delta\nu\int_{-\infty}^{\infty}|\psi(\nu)|^2d\nu.
$$

于是

$$
\int|g^{(1)}(\tau)|^2d\tau
=\frac{s_{\rm spec}}{\Delta\nu},
$$

再考虑两种未分辨偏振模式的稀释 $p_{\rm pol}$，便得到正文第 5.2 节使用的

$$
\boxed{
\tau_{\rm eff}
=\frac{p_{\rm pol}s_{\rm spec}}{\Delta\nu}
}.
$$

窄带近似 $\Delta\nu\simeq c\Delta\lambda/\lambda_0^2$ 给出

$$
\tau_{\rm eff}\simeq
p_{\rm pol}s_{\rm spec}
\frac{\lambda_0^2}{c\,\Delta\lambda}.
$$

因此人们有时口头说“相干时间与 $1/\Delta\lambda$ 成正比”，但严格的傅里叶共轭变量是频率，基本关系是 $1/\Delta\nu$；只有在 $\Delta\lambda\ll\lambda_0$ 时才可使用上面的波长近似。

---

## 附录 B：长曝光 SNR 公式的自洽检查

### B.1 量纲检查

正文第 8.2 节 SNR 公式的率和相干面积组合为

$$
\frac{r_{\star,i}r_{\star,j}\tau_{\rm eff}}
{\sqrt{r_{{\rm tot},i}r_{{\rm tot},j}}}.
$$

分子的量纲是 $\mathrm{s^{-2}}\times\mathrm{s}=\mathrm{s^{-1}}$，分母也是 $\mathrm{s^{-1}}$，所以该项无量纲；$B_{\rm eff}T_{\rm seg}$ 也是 $\mathrm{s^{-1}}\times\mathrm{s}$，因此最终 SNR 无量纲。

### B.2 理想情况下为什么单通道 SNR 对光学带宽近似不敏感

对两台相同望远镜、忽略 NSB 和暗计数时，$r_{\rm tot}=r_\star$，SNR 中与光学频带有关的部分化为

$$
\frac{r_\star^2\tau_{\rm eff}}{r_\star}
=r_\star\tau_{\rm eff}.
$$

由第 4.2 节和附录 A，平坦连续谱下

$$
r_\star=A_{\rm col}\eta_{\rm det}n_\nu\Delta\nu,
\qquad
\tau_{\rm eff}=\frac{p_{\rm pol}s_{\rm spec}}{\Delta\nu}.
$$

两式相乘后 $\Delta\nu$ 抵消：

$$
r_\star\tau_{\rm eff}
=A_{\rm col}\eta_{\rm det}n_\nu
\,p_{\rm pol}s_{\rm spec}.
$$

这正是“窄带降低光子率，却提高单对光子的相干对比度”的定量表达。若恒星和连续谱 NSB 都经过同一理想滤光片并保持线性，它们随带宽的共同缩放也会产生类似抵消。真正改变最优带宽的是谱线和非平坦光谱、波长依赖的 PDE/透过率、暗计数与附加电子噪声、探测器饱和/恢复、色散、滤光片入射角漂移，以及是否能同时读取多个独立光谱通道。

### B.3 时间展宽怎样传到曝光时间

光学时间核不改变 $\tau_{\rm eff}$ 所表示的 HBT 总面积，而是通过 $H_{\rm opt}(f)$ 改变 $B_{\rm eff}$。由于

$$
\mathrm{SNR}\propto\sqrt{B_{\rm eff}T_{\rm seg}},
$$

在其他条件不变时

$$
\frac{\mathrm{SNR}_2}{\mathrm{SNR}_1}
=\sqrt{\frac{B_{{\rm eff},2}}{B_{{\rm eff},1}}},
\qquad
\frac{T_2}{T_1}
=\frac{B_{{\rm eff},1}}{B_{{\rm eff},2}}
$$

分别给出 SNR 比和达到相同 SNR 所需的曝光时间比。当前 DC 时间核把 125 MHz 降到 110.91 MHz，因此匹配 SNR 降低约 5.8%，等效曝光增加约 12.7%。这就是第 6 节光学时间展宽进入第 8.2 节观测结果的明确路径。

---

## 参考物理关系

- [Pupil plane intensity interferometry with imaging air Cherenkov telescopes](https://academic.oup.com/mnras/article/538/2/867/8015798)：给出时间相干面积、电子带宽和 SNR 的统一表达。
- [Intensity interferometry with more than two detectors?](https://academic.oup.com/mnras/article/437/1/798/1007746)：解释窄带下光子率与相干时间的抵消，以及独立光谱通道的收益。
- [Intensity interferometer results on Sirius with 0.25 m telescopes](https://academic.oup.com/mnras/article/537/3/2527/8003771)：讨论实际系统中探测时间分辨远大于光学相干时间时的相关对比度展宽。
