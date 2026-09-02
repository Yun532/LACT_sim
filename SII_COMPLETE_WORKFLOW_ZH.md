# LACT 单像素恒星强度干涉：完整计算流程

本文说明当前统一版本实际执行的完整流程：

```text
源的空间亮度模型
  → 理论复可见度与 |V|²
  → 32 台望远镜随时间变化的 UVW
  → 单像素光学、SiPM 与短波形标定
  → 两台望远镜的时间相关
  → 全阵列长曝光、带噪声的稀疏 UV 测量
  → 无相位图像重建与可靠性检验
```

对应实现：

- 完整可执行 notebook：[`notebooks/lact_sii_paper_simulation.ipynb`](notebooks/lact_sii_paper_simulation.ipynb)
- 统一模拟模块：[`python/sii_unified.py`](python/sii_unified.py)
- 独立重建模块：[`python/sii_reconstruction.py`](python/sii_reconstruction.py)
- 自动测试：[`tests/test_sii_unified.py`](tests/test_sii_unified.py)

本文和图对应 `intensity-interferometry-unified` 分支提交 `46fe0ed`。图由上述 notebook 使用固定随机种子生成，不是手工绘制的示意图。

## 0. 首先明确三个容易混淆的问题

### 0.1 强度干涉测量的不是普通光学干涉条纹

振幅干涉测量电场的一阶相干，而强度干涉测量两个探测器强度涨落的二阶相关：

$$
g_{ij}^{(2)}(\tau)-1=
\frac{\langle\delta I_i(t)\,\delta I_j(t+\tau)\rangle}
{\langle I_i\rangle\langle I_j\rangle}.
$$

对热光，Siegert 关系将它和空间复相干度联系起来：

$$
g_{ij}^{(2)}(\tau)-1
=\beta\,|V(u,v)|^2\,|g^{(1)}(\tau)|^2.
$$

因此最终测量量是 $|V|^2$，不是复可见度 $V$ 本身。

### 0.2 不能把整段波形积分成一个 PE 总数后再做相关

如果只保留

$$
Q_i=\int x_i(t)\,dt,
$$

纳秒尺度的涨落结构就被丢掉了。正确做法是校正几何时延后，对两条去均值的时间序列逐采样相关，或者使用与之等价的细时间格 PE 序列相关。

### 0.3 不需要为一整夜逐光子保存原始 ADC 波形

250 MS/s 连续记录 6 小时，每台望远镜约有

$$
250\times10^6\times21600\simeq5.4\times10^{12}
$$

个采样。32 台望远镜会超过 $10^{14}$ 个采样。当前程序采用两层模拟：

1. 用微秒级完整波形验证“光子/PE → SiPM → SPE → ADC → 相关峰”；
2. 用经该响应校准的长曝光充分统计量模拟小时级 $|V|^2$ 和误差。

在系统近似平稳、线性且传递函数已校准时，这比保存整夜波形更实用，同时保留强度干涉需要的二阶统计量。

---

## 1. 生成源的空间强度模型

第一步定义天空切平面上的归一化亮度分布

$$
I(l,m),
$$

其中 $l,m$ 是相对于视场中心的小角度方向余弦。小视场下可以直接理解为赤经和赤纬方向的角偏移。

当前默认科学案例是两颗均匀圆盘组成的双星：

$$
I(l,m)=I_1(l-l_1,m-m_1)+rI_2(l-l_2,m-m_2),
$$

其中 $r=F_2/F_1$ 是次星与主星的流量比。输入参数包括：

- 总 AB 星等；
- 主星和次星角直径；
- 双星分离角；
- 位置角；
- 两星流量比。

默认案例为 $m_{AB}=2$、分离角 $0.20\,\mathrm{mas}$ 的不等亮双星。离散天图用于画图和重建验证；理论可见度使用解析圆盘公式计算，不依赖天图像素分辨率。

![理论双星源天图](docs/sii_workflow_figures/01_theoretical_source_sky.png)

**图 1。** 理论源的天空强度分布。横轴是 $\Delta\alpha\cos\delta$，纵轴是 $\Delta\delta$，单位均为 mas；色标为归一化像素亮度。两颗星的直径、分离角和流量比均进入后续理论可见度。

---

## 2. 从源天图计算理论 UV 功率平面

van Cittert–Zernike 定理给出归一化复可见度：

$$
V(u,v)=
\frac{
\iint I(l,m)e^{-2\pi i(ul+vm)}\,dl\,dm
}{
\iint I(l,m)\,dl\,dm
}.
$$

$u,v$ 是投影基线除以波长后的无量纲空间频率。

### 2.1 均匀圆盘

角直径为 $\theta$ 的均匀圆盘具有解析可见度：

$$
V_{\rm disk}(q)=\frac{2J_1(\pi\theta q)}{\pi\theta q},
\qquad q=\sqrt{u^2+v^2}.
$$

### 2.2 双星

两颗均匀圆盘的复可见度为

$$
V(u,v)=
\frac{
V_1(q)e^{-2\pi i(ul_1+vm_1)}
+rV_2(q)e^{-2\pi i(ul_2+vm_2)}
}{1+r}.
$$

强度干涉仪直接测量

$$
P(u,v)=|V(u,v)|^2,
$$

而不是 $V$ 的相位。

对于实值天空亮度，

$$
V(-u,-v)=V^*(u,v),
$$

因此

$$
|V(-u,-v)|^2=|V(u,v)|^2.
$$

这是真正的中心对称。图中若视觉上更像左右对称，是阵列投影、源位置角和绘图范围共同造成的外观，不改变数学上的中心共轭关系。

![理论 UV 功率平面](docs/sii_workflow_figures/02_theoretical_uv_power.png)

**图 2。** 双星的连续理论 $|V(u,v)|^2$。左图使用线性色标，右图使用对数色标显示弱条纹和零点。条纹方向由双星位置角决定，条纹间隔主要由双星分离角决定，圆盘直径控制长基线包络衰减。

---

## 3. 生成 32 台望远镜随时间变化的 UVW 覆盖

对望远镜 $i,j$，由本地东、北、天顶坐标得到物理基线：

$$
\mathbf B_{ij}=\mathbf r_j-\mathbf r_i=(B_E,B_N,B_U).
$$

源时角为

$$
H(t)=\mathrm{LST}(t)-\alpha,
$$

其中 $\alpha$ 是源赤经。地球自转改变 $H$，因此一条固定物理基线会在 UV 平面上画出轨迹。

在站点纬度 $\phi$、源赤纬 $\delta$ 下，源方向为

$$
\hat{\mathbf s}=
\begin{bmatrix}
\cos\delta\sin H\\
\sin\delta\cos\phi-\cos\delta\cos H\sin\phi\\
\sin\delta\sin\phi+\cos\delta\cos H\cos\phi
\end{bmatrix}.
$$

天空切平面的两个基向量为

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

投影坐标为

$$
u_m=\mathbf B\cdot\hat{\mathbf u},\qquad
v_m=\mathbf B\cdot\hat{\mathbf v},\qquad
w_m=\mathbf B\cdot\hat{\mathbf s}.
$$

除以观测波长：

$$
u=\frac{u_m}{\lambda},\qquad
v=\frac{v_m}{\lambda}.
$$

两台望远镜之间的几何时延为

$$
\tau_g=\frac{w_m}{c}.
$$

32 台望远镜共有

$$
N_{\rm baseline}=\frac{32\times31}{2}=496
$$

条独立物理基线。默认 6 小时观测、每 20 分钟一个时间段，共 18 个时间格，因此一晚产生

$$
496\times18=8928
$$

条独立基线—时间记录。

![随观测时间累积的 UV 覆盖](docs/sii_workflow_figures/03_uv_coverage_cumulative_time.png)

**图 3。** 以中天为中心，从一个 20 分钟时间格逐步累积到 6 小时的 UV 覆盖。负基线共轭点以相同显著度显示，但不计作新的独立测量。延长同一晚的观测既增加积分时间，也延长 UV 轨迹；重复相同恒星时角的多晚观测主要提高 SNR，不会持续产生全新的轨迹。

当前案例已用每个 20 分钟时间段内 41 个子时刻检查中点近似。段平均与中点 $|V|^2$ 的最大绝对差约为 0.0072，RMS 约为 0.0011；参考情形单点统计误差约为 0.099，因此当前近似足够。更小的源、更长的基线或更长时间格应改为段内积分。

---

## 4. 单台望远镜的光学、光子、SiPM 和波形

这一部分分成“光学响应标定”和“短波形生成”，不能理解为给整夜每个光子都重复完整光追。

### 4.1 单像素光学标定

当前使用 `main` 的光学过程，对 100 万条轴上、400 nm、纯平行光进行标定，包含：

- LACT2 镜面几何和镜面误差；
- 结构遮挡；
- 单像素集光器；
- 滤光片；
- SiPM PDE；
- 光子到达时间分布。

标定输出被压缩为：

1. 有效面积和总通光效率；
2. 单像素光学到达时间核 $K_{\rm opt}(t)$。

该时间核是完整光学模拟的输出，不是实测到达时间分布，也只代表当前轴上、400 nm 和当前姿态。

### 4.2 恒星探测率

AB 星等转换为光谱流量密度：

$$
F_\nu=3631\,\mathrm{Jy}\times10^{-0.4m_{AB}}.
$$

单光子能量为

$$
E_\gamma=\frac{hc}{\lambda},
$$

光子谱密度为

$$
n_\nu=\frac{F_\nu}{E_\gamma}.
$$

窄带宽下

$$
\Delta\nu\simeq\frac{c\Delta\lambda}{\lambda^2},
$$

因此探测恒星光子率为

$$
r_\star=A_{\rm eff}\eta n_\nu\Delta\nu.
$$

### 4.3 p.e. 点过程

在短时间窗 $T$ 内，恒星、NSB 和暗计数分别按 Poisson 过程生成：

$$
N_\star\sim\mathrm{Poisson}(r_\star T),
$$

$$
N_{\rm NSB}\sim\mathrm{Poisson}(r_{\rm NSB}T),
$$

$$
N_{\rm dark}\sim\mathrm{Poisson}(r_{\rm dark}T).
$$

`main` 兼容的逐 p.e. 表包含到达时间、望远镜、像素、SiPM 位置、波长和来源标签。当前是轴上单像素、400 nm 案例，因此探测后的 p.e. 不再保存天空方向；方向信息已经进入有效面积、光学时间核和几何时延。

### 4.4 p.e. 到 ADC 波形

每个 p.e. 随机分配到一个 SiPM 微单元。同一微单元两次击中的间隔为 $\Delta t$ 时，恢复电荷比例为

$$
f_{\rm rec}=1-\exp\left(-\frac{\Delta t}{\tau_{\rm rec}}\right),
$$

当前按要求采用

$$
\tau_{\rm rec}=10\,\mathrm{ns}.
$$

随后从实测单 PE 电荷样本抽取涨落，并和实测约 170 ns 的 SPE 模板卷积：

$$
x_i(t)=\sum_k q_k f_{{\rm rec},k}
h_{\rm SPE}(t-t_k)+n_i(t).
$$

最后加入电子学噪声并进行 ADC 采样和量化。当前 notebook 展示 2 µs 记录，并在分析窗两侧生成足够长的额外 p.e.，避免 170 ns SPE 尾部在边缘被人为截断。

![实测 SPE、两台望远镜短波形和单记录相关](docs/sii_workflow_figures/04_short_waveform_and_recovery.png)

**图 4。** 左：仓库中的实测公共 SPE 模板。中：两台望远镜的一段 2 µs ADC 等效波形，包含恒星、NSB、SiPM 恢复和电荷涨落。右：单条短记录的互相关。单条 2 µs 记录明显由随机噪声主导，因此它用于验证电子学流程，不能用来声称已经探测到恒星 HBT 信号。

---

## 5. 两台望远镜如何得到一个 UV 点的 $|V|^2$

### 5.1 先校正几何时延

对基线 $ij$，先由 $w$ 计算

$$
\tau_g=\frac{w}{c},
$$

然后把两条波形对齐：

$$
x_j'(t)=x_j(t+\tau_g).
$$

### 5.2 对去均值时间序列做相关

离散估计量为

$$
\widehat C_{ij}(\tau_k)=
\frac{1}{N_k}
\sum_n
[x_i(n)-\bar x_i]
[x_j(n+k)-\bar x_j].
$$

归一化后得到 $g^{(2)}(\tau)-1$。对校正后零延迟附近的相关响应做拟合或匹配积分，即可估计 $|V(u,v)|^2$。

### 5.3 波形相关与细时间格 PE 相关

若电子学近似线性时不变，

$$
x_i(t)=h_i(t)*n_i(t)+e_i(t),
$$

波形相关就是 PE 时间序列相关经过电子学传递函数后的结果。只要 SPE、光学时间核、ADC 采样和电子噪声被正确标定，两种做法在二阶统计量上等价。把整段波形积分为一个总 PE 数则不等价。

### 5.4 模拟中的相关光子对

程序用稀疏相关 p.e. 对或共享 Poisson 分量保存热光正确的二阶矩，超额相关对率为

$$
R_{\rm pair}=r_{\star,i}r_{\star,j}\tau_c|V(u,v)|^2.
$$

这里的“相关对”是 Monte Carlo 表示，不表示一颗光子同时被两台望远镜探测。HBT 相关来自热光强度涨落和光子聚束。

![几何时延校正前后的理论电子学相关峰](docs/sii_workflow_figures/05_two_telescope_g2_peak.png)

**图 5。** 左：几何时延校正前，相关峰位于非零 lag；校正后被移到零延迟。右：零延迟附近的电子学展宽相关峰。峰高约为 $10^{-6}$ 量级，因此真实短波形通常看不到稳定峰，需要大量时间样本累积。

### 5.5 小时级测量的充分统计量

整夜不生成原始 ADC 数组，而使用短波形和仪器响应校准后的强度干涉 SNR：

$$
\mathrm{SNR}_{ij}=
\frac{r_{\star,i}r_{\star,j}}
{\sqrt{r_{{\rm tot},i}r_{{\rm tot},j}}\,\Delta\nu_{\rm opt}}
|V_{ij}|^2
\sqrt{\frac{B_{\rm eff}T}{2}}
\frac{\eta_{\rm elec}\eta_{\rm timing}}{F_{\rm EN}},
$$

其中

$$
r_{\rm tot}=r_\star+r_{\rm NSB}+r_{\rm dark}.
$$

所以单位可见度的统计误差为

$$
\sigma(|V|^2)=
\left[
\frac{r_{\star,i}r_{\star,j}}
{\sqrt{r_{{\rm tot},i}r_{{\rm tot},j}}\,\Delta\nu_{\rm opt}}
\sqrt{\frac{B_{\rm eff}T}{2}}
\frac{\eta_{\rm elec}\eta_{\rm timing}}{F_{\rm EN}}
\right]^{-1}.
$$

NSB、暗计数和独立电子噪声一般不会产生恒星的跨镜相关峰，但会增加方差。光学时间展宽和残余时钟误差会降低相关峰幅度或有效带宽。

当前长曝光模拟为每个“时间段 × 望远镜”生成一次恒星、NSB 和暗计数 Poisson 计数；同一台望远镜的计数、增益和时钟状态由它参与的全部基线共享。这样 TEL.1 的状态会同时影响 1–2、1–3 等基线，而不是把每条基线错误地当成完全独立仪器。

---

## 6. 组合所有望远镜和时间段，得到真实的稀疏 UV 测量

对每个“望远镜对 × 时间段”，程序保存

$$
(u_k,v_k,w_k,\tau_{g,k},
\widehat{|V_k|^2},\sigma_k).
$$

模拟测量可写成

$$
\widehat{|V_k|^2}
=|V_k|^2G_iG_jA_{{\rm timing},ij}
+Z_{ij}+\epsilon_k,
$$

其中：

- $G_i,G_j$ 是望远镜增益和透明度误差；
- $A_{{\rm timing},ij}$ 是残余时延引起的相关衰减；
- $Z_{ij}$ 是每条基线的零点误差；
- $\epsilon_k$ 是统计噪声。

真实观测不是连续填满的 UV 图片，而是一组稀疏、带误差的采样点。低 SNR 时可能得到

$$
\widehat{|V|^2}<0
\quad\text{或}\quad
\widehat{|V|^2}>1.
$$

程序不会把它们裁剪到 $[0,1]$，因为裁剪会制造低 SNR 正偏差。

![模拟的稀疏 UV 测量和测量闭合](docs/sii_workflow_figures/06_simulated_uv_measurements.png)

**图 6。** 左：默认 $m_{AB}=2$ 双星在实际 32 镜、6 小时覆盖上的带噪声 $|V|^2$；颜色值未裁剪。右：模拟测量与注入真值的闭合关系，参考情形每 20 分钟段的单位可见度 SNR 约 10.1，对应 $\sigma(|V|^2)\simeq0.099$。

重建前，对落入同一 UV 单元的点进行逆方差加权：

$$
\bar y=
\frac{\sum_k y_k/\sigma_k^2}
{\sum_k1/\sigma_k^2},
$$

$$
\sigma_{\bar y}=
\left(\sum_k\frac{1}{\sigma_k^2}\right)^{-1/2}.
$$

当前默认合并单元宽度为 $120\,\mathrm{M}\lambda$。添加的 $(-u,-v)$ 共轭点只用于对称显示，不作为第二次独立观测。

---

## 7. 从无相位 $|V|^2$ 重建天空图像

直接对 $|V|^2$ 做逆傅里叶变换得到的是源亮度的自相关，而不是源天图本身。因此主方法直接在图像空间拟合稀疏、带符号、带不确定度的功率测量。

### 7.1 图像约束

设支撑区域内像素亮度为 $x_p$。用 softmax 参数化保证

$$
x_p\ge0,\qquad\sum_px_p=1.
$$

同时限制图像只存在于有限圆形支撑区域。

### 7.2 正向模型

第 $k$ 个 UV 点的复场预测为

$$
F_k=\sum_px_p
e^{-2\pi i(u_kl_p+v_km_p)},
$$

功率预测为

$$
P_k=|F_k|^2.
$$

### 7.3 优化目标

$$
\mathcal L=
\frac{\sum_kw_k\rho_{\rm Huber}(P_k-y_k)}
{\sum_kw_k}
+\lambda\sum_{\langle p,q\rangle}(x_p-x_q)^2,
$$

其中 $w_k\propto1/\sigma_k^2$。Huber 损失降低异常 UV 点的影响，相邻像素正则抑制噪声尖峰。

优化采用多起点 L-BFGS-B，选择目标函数最小的结果。重建程序只读取 $(u,v,|V|^2,\sigma)$，不读取真实双星参数或真实傅里叶相位。

### 7.4 固有歧义

只测量傅里叶模平方时存在：

1. 绝对平移歧义；
2. 180° 镜像歧义。

程序将重建图质心移动到中心以固定平移规范。180° 镜像只能在注入真值验证时对齐，不能由强度干涉数据自身消除。

![无相位图像重建结果](docs/sii_workflow_figures/07_phaseless_reconstruction.png)

**图 7。** 上排依次为合并后的 UV 功率、自相关脏图和正向闭合；下排为统一角分辨率下的真实双星、歧义对齐后的无相位重建和残差。本例正向加权 RMSE 为 0.042，重建与真值的相关系数为 0.972。真值只用于最后验证和歧义对齐，不进入优化器。

Fienup HIO+ER 在 notebook 中作为独立基准，但主结果采用上面的正值、有限支撑、误差加权直接拟合。

---

## 8. 星等极限和角分辨率怎样得到

### 8.1 角分辨率

最简单的尺度估计为

$$
\theta_{\rm nominal}\sim\frac{\lambda}{B_{\rm max,proj}}.
$$

但真正能否分开双星还取决于：

- UV 方向覆盖；
- $|V|^2$ 条纹是否被实际基线采到；
- 星等、NSB 和有效带宽；
- 有效积分时间；
- 标定误差；
- 重建先验和检测判据。

所以不能仅凭最长基线给出可靠成像结论。

### 8.2 星等极限

星等极限不是根据某一张“看起来像双星”的图决定。notebook 对不同星等和积分时间运行多次独立噪声实现，并同时模拟：

- 双星注入样本；
- 单星零假设样本；
- 双星与单星参数模型的 BIC；
- 不读取真值的非参数重建；
- 真阳性率和假阳性率及其置信区间。

只有当双星恢复率足够高、单星假阳性率足够低，并且两者置信区间能够区分时，才把该条件标为可靠。单张漂亮图不能作为星等极限证据。

---

## 9. 当前包含和仍缺失的实测化因素

| 项目 | 当前状态 | 如何进入模拟 |
|---|---|---|
| 32 镜坐标 | 已有 | 生成 496 条物理基线和随时间 UVW |
| 400 nm 单像素完整光学 | 已有模拟标定 | 有效面积、通光效率、到达时间核 |
| 实测 SPE 模板 | 已有 | 逐 p.e. 波形卷积和有效带宽 |
| 实测单 PE 电荷涨落 | 已有 | 逐 p.e. 电荷抽样和 excess-noise factor |
| 4 ns ADC 采样 | 已有 | 波形采样和相关器带宽限制 |
| NSB 光谱/探测率 | 已有 | 独立 Poisson 背景和统计误差 |
| 微单元恢复 | 接口已有 | 当前固定为 10 ns |
| 光子数 Poisson 涨落 | 已有 | 逐时间段、逐望远镜抽样 |
| 望远镜共享增益/时钟状态 | 已有接口 | 同一望远镜参与的基线共享 |
| 串扰、后脉冲 | 暂无可信输入 | 接口保留，当前为 0 |
| 实测暗计数率 | 暂无可信输入 | 接口保留，当前为 0；非零时按段/望远镜 Poisson 抽样 |
| 实测电子学加性噪声 | 暂无可信输入 | 接口保留，当前为 0 |
| 实测 ADC 非线性/量化参数 | 暂无可信输入 | 接口保留，当前为 0 或关闭 |
| 实测时钟稳定性 | 暂无可信输入 | 接口保留，当前残差为 0 |
| 滤光片角响应 | 暂无可信输入 | 接口保留，当前为 0 |
| 不同仰角/离轴角/波长时间核 | 尚未标定 | 当前结果只代表轴上 400 nm 当前姿态 |

因此，当前结果应表述为“使用现有 LACT 光学、SPE、NSB 和阵列输入得到的、带明确缺失项的研究预测”，不能表述为 LACT 已实测最终性能。

---

## 10. 一键运行入口与复现关系

Python 的最小闭环入口为：

```python
result = run_sii_pipeline(
    layout=layout,
    source=source,
    observation=observation,
    instrument=instrument,
    do_reconstruction=True,
)
```

内部顺序固定为：

```python
uvw = generate_uvw(layout, observation, instrument)
measurements, metadata = simulate_uv_observation(
    uvw, source, observation, instrument
)
reconstruction = reconstruct_uv(measurements)
```

其中：

- `source` 控制星等、直径、分离角、位置角和流量比；
- `observation` 控制赤纬、站点纬度、每晚时间、时间格和夜数；
- `instrument` 控制波长、带宽、面积、效率、NSB、SiPM、SPE、ADC 和标定误差；
- `layout` 可以替换为未来更新后的望远镜坐标；
- 重建器与模拟器独立，只依赖带误差的 UV 测量表。

论文级 notebook 还额外生成不同时间 UV 覆盖、不同 NSB、不同观测时间、不同星等、多次噪声实现、单星零假设以及参数/非参数重建对照。

## 11. 本文七张图的复现来源

| 本文图 | notebook 输出文件 | 对应步骤 |
|---|---|---|
| 图 1 | `theoretical_source_sky.png` | 源模型 |
| 图 2 | `theoretical_uv_power.png` | 理论 $|V|^2$ |
| 图 3 | `uv_coverage_cumulative_time.png` | 32 镜 UVW 覆盖 |
| 图 4 | `short_waveform_and_recovery.png` | 单像素波形与 SiPM 响应 |
| 图 5 | `two_telescope_g2_peak.png` | 两镜几何时延与相关峰 |
| 图 6 | `simulated_uv_measurements.png` | 带噪声稀疏 UV 测量 |
| 图 7 | `phaseless_reconstruction.png` | 无相位重建与闭合 |

运行 [`notebooks/lact_sii_paper_simulation.ipynb`](notebooks/lact_sii_paper_simulation.ipynb) 后，原始输出写入 `run_logs/sii_paper_notebook/`。本文 `docs/sii_workflow_figures/` 中的七张图是这些输出的同版本副本，目的是保证 Markdown 在 GitHub 或本地预览时可以直接显示。
