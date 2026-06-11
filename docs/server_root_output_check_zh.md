# 服务器 ROOT 输出与画图检查流程

本文档用于在服务器上安装/加载 ROOT，编译 LACT_sim 的 `lact_event_root_v1`
输出，并用 quick-look 图检查输出内容。

## 1. 推荐环境方式

优先使用服务器已有 module：

```bash
module avail root
module load root
root-config --version
root-config --prefix
```

如果 `module load root` 报错，但 `root-config --version` 能输出版本，说明 ROOT 已经在当前
`PATH` 中可用，只是服务器没有提供 modulefile。可以继续使用当前 ROOT。

如果服务器没有 module，建议用独立 conda/mamba 环境，避免污染系统环境：

```bash
mamba create -n lact-root -c conda-forge root cmake compilers hdf5 zlib uproot matplotlib numpy
mamba activate lact-root
root-config --version
```

如果服务器只有 conda：

```bash
conda create -n lact-root -c conda-forge root cmake compilers hdf5 zlib uproot matplotlib numpy
conda activate lact-root
```

ROOT 建议版本：`>= 6.24`。

## 2. 编译 LACT_sim

在 LACT_sim 仓库根目录：

```bash
make hessio

cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DLACT_ENABLE_ROOT=ON \
  -DCMAKE_PREFIX_PATH="$(root-config --prefix)"

cmake --build build --target run_corsika_trace -j"${SLURM_CPUS_PER_TASK:-8}"
```

确认 CMake 输出里有类似：

```text
Found ROOT 6.xx; lact_event ROOT output enabled
```

如果看到：

```text
ROOT >= 6.24 not found; lact_event ROOT output disabled
```

说明当前 shell 里 ROOT 没有被正确加载。检查：

```bash
which root
which root-config
root-config --prefix
```

部分 ROOT 版本的 `root-config` 没有 `--cmakedir` 参数，即使版本很新也可能如此。
这时不要用 `ROOT_DIR="$(root-config --cmakedir)"`，改用：

```bash
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DLACT_ENABLE_ROOT=ON \
  -DCMAKE_PREFIX_PATH="$(root-config --prefix)"
```

如果仍找不到 ROOT，再定位 CMake 配置文件：

```bash
find "$(root-config --prefix)" -name ROOTConfig.cmake -o -name root-config.cmake
```

假设输出类似 `/path/to/root/lib/cmake/ROOT/ROOTConfig.cmake`，则：

```bash
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DLACT_ENABLE_ROOT=ON \
  -DROOT_DIR=/path/to/root/lib/cmake/ROOT
```

## 3. 配置 ROOT 输出

在运行 `run_corsika_trace` 的配置中加入：

```ini
output.lact_root_enabled=true
output.lact_root_path=run_logs/my_corsika_run/lact_events.root
output.lact_profile=timeseries_pe
```

如果只需要积分图像：

```ini
output.lact_profile=image_pe
```

如果要输出 p.e. 时间序列，需要同时打开 p.e. waveform：

```ini
waveform.enabled=true
waveform.source=pe
waveform.time_bin_width_ns=1.0
waveform.time_window_start_ns=0.0
waveform.time_window_end_ns=100.0
```

运行示例：

```bash
./build/run_corsika_trace configs/your_run.cfg
```

结束输出里应出现：

```text
lact_root_path=run_logs/my_corsika_run/lact_events.root
```

## 4. 画图检查

仓库已提供可直接运行的 ROOT quick-look 配置。只需要把输入 EventIO/simtel 文件作为第二个参数：

```bash
./build/run_corsika_trace \
  configs/examples/corsika_lact_root_quicklook.cfg \
  /path/to/input.simtel.zst
```

默认只跑 `source.max_shower_events=1`，输出在：

```text
run_logs/lact_root_quicklook/lact_events.root
run_logs/lact_root_quicklook/corsika_trace.h5
run_logs/lact_root_quicklook/corsika_trace_summary.csv
```

如果要用最新 full-response 链路、但不引入 NSB，推荐改用：

```bash
./build/run_corsika_trace \
  configs/examples/corsika_lact_root_full_response.cfg \
  /path/to/input.simtel.zst
```

这个配置基于 `configs/official_tests/corsika_full_response_camera.cfg`，保留结构变形、
光学误差、反射率/滤光片/SiPM PDE、new_camera 和 multiplicity trigger；同时用
`configs/nsb/ideal.cfg` 关闭 NSB，并输出：

```text
run_logs/lact_root_full_response/lact_events.root
run_logs/lact_root_full_response/corsika_trace.h5
run_logs/lact_root_full_response/corsika_trace_summary.csv
```

脚本：

```bash
python/plot_lact_root_output.py run_logs/my_corsika_run/lact_events.root \
  --outdir run_logs/my_corsika_run/root_quicklook
```

指定事件和望远镜：

```bash
python/plot_lact_root_output.py run_logs/my_corsika_run/lact_events.root \
  --event-id 1 \
  --tel-id 1 \
  --outdir run_logs/my_corsika_run/root_quicklook_event1_tel1
```

如果缺 Python 画图库：

```bash
mamba install -c conda-forge uproot matplotlib numpy
# 或
python -m pip install uproot awkward matplotlib numpy
```

输出内容：

```text
summary.txt
event_<event>_tel_<tel>_image_pe.png
event_<event>_tel_<tel>_peak_time.png
event_<event>_tel_<tel>_waveform_sum.png      # 有 waveforms 时
total_pe_hist.png
corsika_truth_hist.png
```

重点检查：

```text
image_pe 图：相机上是否有合理 shower image，不应全黑或全在一个异常 pixel。
peak_time 图：有 waveform 时应有有限值，时间梯度不应明显乱跳。
waveform_sum 图：timeseries_pe 时应有 p.e. 随时间的脉冲结构。
total_pe_hist：总 p.e. 分布不应全 0。
corsika_truth_hist：energy/altitude/Xmax 应与输入样本范围一致。
summary.txt：selected_image_sum_pe 和 waveforms sum_pe 应同量级。
```

## 5. pylast 读取检查

安装/编译 pylast 时同样需要 ROOT 环境。Python 侧预期用法：

```python
from pylast.io import LactEventSource

source = LactEventSource("run_logs/my_corsika_run/lact_events.root", max_events=10)
event = source[0]

print(event.event_id)
print(event.simulation.shower.energy)
print(event.r1.tels.keys())
print(event.dl0.tels.keys())
```

`timeseries_pe` 文件中：

```python
tel_id = sorted(event.r1.tels.keys())[0]
r1 = event.r1.tels[tel_id]
print(r1.waveform.shape)  # n_pixels x n_time_bins
```

完整 notebook 流程见：

```text
notebooks/lact_root_to_pylast_visualize.ipynb
```

它从一个 CORSIKA/EventIO/simtel 输入文件开始，依次完成 LACT_sim ROOT 输出、
`pylast.io.LactEventSource` 读取，以及 `pylast.visualize` 画图。

`image_pe` 文件中：

```python
print(r1.waveform.shape)  # n_pixels x 1
```

## 6. 常见问题

### ROOT 找不到

```bash
source /path/to/root/bin/thisroot.sh
cmake -S . -B build -DROOT_DIR="$(root-config --cmakedir)"
```

### 服务器不能联网

先在能联网的机器创建 conda pack 或请管理员安装 module。不要把 ROOT 放进
`external/`，项目只通过 `root-config`/`ROOT_DIR` 查找系统 ROOT。

### 文件有 ROOT 但没有 waveforms

检查：

```ini
output.lact_profile=timeseries_pe
waveform.enabled=true
waveform.source=pe
```

如果是 `image_pe` profile，只有 `observations.image_pe`，R1 会作为
`n_pixels x 1` 的 single-sample p.e. image 读取。
