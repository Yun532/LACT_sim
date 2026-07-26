# Photon CSV 从零运行指南

这份指南完成一条最短但完整的用户流程：安装环境、拉取 LACT_sim 与 pyLAST、编译、读取 Photon CSV、生成 LACT ROOT，并画出相机 p.e. 图。以下命令面向 Ubuntu x86_64。

## 1. 安装环境

已有 `conda` 时跳过 Miniforge 安装，只执行创建环境的命令。

```bash
sudo apt-get update
sudo apt-get install -y git curl

curl -L -o Miniforge3.sh \
  https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3.sh -b -p "$HOME/miniforge3"
source "$HOME/miniforge3/etc/profile.d/conda.sh"

conda create -y -n lact-photoncsv -c conda-forge \
  python=3.11 cmake make cxx-compiler "root>=6.24" hdf5 zlib pip
conda activate lact-photoncsv
```

以后重新登录服务器，只需执行：

```bash
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate lact-photoncsv
```

## 2. 下载并编译

```bash
mkdir -p lact-photoncsv-work
cd lact-photoncsv-work

git clone --branch user_v2.0 --single-branch \
  https://github.com/Yun532/LACT_sim.git
git clone --branch lact_sim --single-branch \
  https://github.com/Yun532/pylast.git

python -m pip install ./pylast
python -m pip install -r ./LACT_sim/requirements.txt

cd LACT_sim
make
```

确认 LACT ROOT 与 pyLAST 都可用：

```bash
test -x build/run_corsika_trace
python -c "from pylast.io import LactEventSource; print('pyLAST OK')"
```

编译输出中应出现 `LACT ROOT output enabled`。如果显示 `disabled`，说明当前环境没有被 CMake 找到的 ROOT 6.24+，此时不能生成本流程需要的 `.root` 文件。

## 3. 直接运行仓库自带示例

仓库已经提供输入 CSV 和完整配置，因此第一次运行不需要改任何文件：

```bash
build/run_corsika_trace configs/examples/photon_csv_source.cfg

python python/plot_photon_csv_root_pylast.py \
  run_logs/photon_csv_source/lact_events.root \
  --event-id 1909 \
  --telescope-id 19 \
  --output run_logs/photon_csv_source/pylast_camera.png
```

成功后得到：

- `run_logs/photon_csv_source/lact_events.root`：LACT ROOT 事件文件；
- `run_logs/photon_csv_source/pylast_camera.png`：pyLAST 绘制的相机 p.e. 图。

## 4. 换成自己的 CSV

复制示例配置，并只修改输入路径；所有命令仍需从 `LACT_sim` 仓库根目录运行。

```bash
cp configs/examples/photon_csv_source.cfg my_photon_csv.cfg
```

编辑 `my_photon_csv.cfg` 中这一行：

```ini
source.csv_path=/path/to/my_photons.csv
```

然后运行：

```bash
build/run_corsika_trace my_photon_csv.cfg
```

输出路径仍由配置中的 `output.lact_root_path` 决定。若修改了 `event_id` 或 `telescope_id`，画图命令中的 `--event-id` 和 `--telescope-id` 也要保持一致。

## 5. CSV 输入格式

CSV 第一行必须是表头。最小格式只有六列，每一行表示一个 photon bunch：

```csv
x_m,y_m,z_m,dir_x,dir_y,dir_z
3.9014,0.1662,0,-0.33790,-0.01472,-0.94107
0.0847,1.3690,0,-0.33653,-0.01435,-0.94156
```

必需列：

| 列名 | 含义 |
| --- | --- |
| `x_m,y_m,z_m` | photon bunch 位置，单位 m |
| `dir_x,dir_y,dir_z` | 光子传播方向；程序会自动归一化，但不能为零向量 |

当前示例配置使用 `source.coordinate_frame=corsika_nwu_relative`：全局方向为 North-West-Up，位置相对 telescope；CORSIKA 二维 bunch 使用 `eventio_2d=true` 和 `source.eventio_reference_z_m=-16`。来自其他坐标系的数据不能只改列名，必须同时修改配置中的 `source.coordinate_frame`。

常用的完整表头可以写成：

```csv
x_m,y_m,z_m,dir_x,dir_y,dir_z,time_ns,wavelength_nm,multiplicity,event_id,telescope_id,eventio_2d
3.9014,0.1662,0,-0.33790,-0.01472,-0.94107,12.5,400,1,1909,19,true
```

所有可选列及当前示例配置下的缺省行为：

| 可选列 | 含义 | 未提供时 |
| --- | --- | --- |
| `time_ns` | 到达时间，ns | `source.time_ns`，默认 0 |
| `wavelength_nm` | 波长，nm，必须大于 0 | `source.wavelength_nm=400` |
| `raw_wavelength_nm` | 原始输入中的波长值 | 等于 `wavelength_nm` |
| `weight` | photon bunch 权重，必须不小于 0 | `source.photon_weight=1` |
| `optical_efficiency_preapplied` | 输入是否已应用前级光学效率 | `false` |
| `multiplicity` | 该 bunch 代表的光子数，必须不小于 0 | `source.multiplicity=1` |
| `event_id` | 输出事件编号 | `source.event_id=1909` |
| `telescope_id` | 望远镜编号 | `source.telescope_id=19` |
| `shower_event_id` | shower 编号 | 等于 `event_id` |
| `array_id` | array 编号 | 0 |
| `source_bunch_index` | 原始 bunch 序号 | 从 0 开始的行序号 |
| `eventio_2d` | 是否来自 EventIO/CORSIKA 二维 bunch | `source.eventio_2d=true` |
| `emission_altitude_km` | 发射高度，km | 留空；仅启用相应大气模型时需要 |

布尔列接受 `1/0`、`true/false`、`yes/no` 或 `on/off`。最小六列格式适合所有行共享相同波长、事件号和望远镜号的情况；这些共享值直接在 cfg 中设置即可。

## 6. 使用自己数据时必须核对

1. `telescope.pointing_el_deg = 90 - CORSIKA zenith_deg`。
2. CORSIKA 二维 bunch 的位置和方向应使用与配置一致的 NWU 约定，下降光子的 `dir_z` 通常为负。
3. CSV 中若包含多个事件或 telescope，建议显式提供 `event_id` 和 `telescope_id` 两列。
4. 第一次先关闭 NSB 和 trigger；示例配置已经设置 `nsb.enabled=false`、`trigger.enabled=false`，便于直接检查光学响应。
5. 示例使用 `response.mode=stochastic_pe`，所以相同输入的 p.e. 接受结果由 `response.seed` 保证可复现。

示例输入文件是 [`configs/sources/event1909_tel19_minimal_photons.csv`](configs/sources/event1909_tel19_minimal_photons.csv)，对应配置是 [`configs/examples/photon_csv_source.cfg`](configs/examples/photon_csv_source.cfg)。
