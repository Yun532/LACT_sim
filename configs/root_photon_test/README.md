# ROOT 光子文件独立相机测试

这个测试读取一个只有 `photons` tree 的 ROOT 文件，把光子转换成 LACT_sim
`PhotonCsv`，然后复用现有 `run_optical_sim` 得到完美情况下的真实相机像素图像。

## 输入 ROOT 格式

默认读取 tree:

```text
photons
```

默认 branch:

```text
photon_x, photon_y, photon_z
photon_u, photon_v
photon_weight
photon_lambda
```

当前脚本假设：

- `photon_x/y/z` 单位是 cm，位于 CORSIKA-like NWU 全局坐标。
- `photon_u/v` 是水平方向余弦。
- `photon_weight` 是该 photon bunch 的 multiplicity。
- `photon_lambda` 单位是 nm。
- `dir_z = -sqrt(1 - u^2 - v^2)`，即光子向下传播。

## 望远镜配置

望远镜信息来自你给的截图，写在：

```text
configs/root_photon_test/telescope_test_root.cfg
```

截图中原始信息是：

```text
telX=-14464.5 cm, telY=-2873.45 cm, telZ=119.7 cm
telzen=25.802 deg, telazi=100.156 deg
```

```ini
position_m=-144.645,-28.7345,1.197
pointing_az_deg=100.156
pointing_el_deg=64.198
```

这里 `pointing_el_deg = 90 - telzen = 64.198 deg`。

## 坐标处理

`test.root` 本身没有完整的阵列/core 元数据，所以这个独立测试采用以下方式：

1. 把 ROOT 里的全局光子位置和方向投影到该望远镜的本地入瞳平面 `local z=0`。
2. 默认用投影后光子束的中位数作为中心，把光子束平移到望远镜入瞳中心。
3. 只保留入瞳平面半径 `8 m` 内的 photon bunch。
4. 之后交给正常 `PhotonCsv` 光追流程。

因此这个测试适合检查“这个 ROOT 光子角分布/波长/权重进入理想望远镜后，相机图像长什么样”，
不是阵列绝对 core/telescope 几何验证。后者需要 ROOT 文件里提供完整阵列和 core 信息。

## 一键运行

```bash
mkdir -p run_logs/root_photon_test
tools/run_root_photon_test.sh /Users/yun/Downloads/test.root \
  2>&1 | tee run_logs/root_photon_test/full_run.log
```

输出：

```text
run_logs/root_photon_test/photons_local.csv
run_logs/root_photon_test/root_conversion.log
run_logs/root_photon_test/run.log
run_logs/root_photon_test/camera_pixels.csv
run_logs/root_photon_test/camera_dense.h5
run_logs/root_photon_test/camera_pe.png
```

## 手动分步运行

转换 ROOT：

```bash
python3 python/root_photons_to_photon_csv.py \
  /Users/yun/Downloads/test.root \
  run_logs/root_photon_test/photons_local.csv \
  --metadata-output run_logs/root_photon_test/root_conversion.log
```

光追：

```bash
build/run_optical_sim configs/root_photon_test/root_photon_camera.cfg \
  2>&1 | tee run_logs/root_photon_test/run.log
```

转成 dense HDF5：

```bash
python3 python/export_trace_hdf5.py \
  --pixel-csv run_logs/root_photon_test/camera_pixels.csv \
  --config configs/root_photon_test/root_photon_camera.cfg \
  --camera-csv configs/cameras/new_camera_pixels.csv \
  --storage dense \
  --output run_logs/root_photon_test/camera_dense.h5
```

画相机：

```bash
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_hdf5_camera.py \
  run_logs/root_photon_test/camera_dense.h5 \
  --event-id 0 \
  --telescope-id 0 \
  --quantity pe \
  --output run_logs/root_photon_test/camera_pe.png
```
