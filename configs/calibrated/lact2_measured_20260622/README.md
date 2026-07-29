# LACT2 实测约束光学配置（20260622）

这组文件是当前推荐的 LACT2 生产光学参数，不覆盖仓库中保留的理想镜面、仅支架
形变以及人工随机误差基准。

运行时的镜面真值表是 `mirror_elevation_series_20260622.csv`。它包含 11 个仰角
锚点、每个仰角 54 片镜子，共 594 行。生产表只保留 9 个运行字段：

- `elevation_deg`、`id`；
- `center_x/y/z`：该仰角下经过支架形变后的镜片中心；
- `normal_x/y/z`：支架形变法向与固定逐镜指向偏转的合成结果；
- `radius_of_curvature`：20260607 汇总表中的逐镜实测曲率半径。

镜片实验编号、现场位置、单镜最小 D80、波段平均反射率和拟合偏转分量只用于
实验追溯，保留在原始表格与拟合结果中，不放入生产 CSV。球面、六边形口径和
`0.8 m` 尺寸等全镜面常量统一写在 `mirror_20260622.cfg` 中；逐镜反射率缩放
缺省为 `1`，额外随机法向误差缺省为 `0`。全镜面统一的等效粗糙度
`0.00725416907642 deg` 通过 `errors.cfg` 施加。

固定法向偏转不是直接测量值，而是使用 20260622 白板 FITS 数据对 11 个仰角
联合拟合得到的 w30 平衡解。它在各仰角保持不变，再与各仰角的支架形变叠加。
本表于 20260729 使用统一后的相机坐标系（`u = -mirror-local x`，
`v = +mirror-local y`）和实验白板原始 `u/v` 坐标重新联合拟合、重新光线追迹
后覆盖。每个仰角使用 300000 个入射光子时，11 个仰角的平均图像 NRMSE 为
`0.15045`、平均相关系数为 `0.98442`、D80 RMSE 为 `0.24267 mm`；逐镜固定
法向偏转的径向 RMS 为 `0.01631 deg`，最大值为 `0.03526 deg`。

## 反射率

绝对反射率不写成逐镜常数，而由
`configs/efficiency/mirror_reflectivity_dm0113_13point_mean.csv` 提供随波长变化的
全局实测曲线。该文件来自 DM0113 镜面 13 个测点的逐波长平均，已经在此前提交
中上线；本次没有修改其数值或文件名。

## 运行入口

- `configs/examples/lactroot_only.cfg`：CORSIKA/EventIO 生产默认入口，只生成
  pyLAST 所需的 ROOT；
- `configs/examples/lact2_measured_parallel.cfg`：同一套镜面参数的平行光白板测试。

`mirror_20260622.cfg` 直接读取上述单张仰角序列表。`errors.cfg` 只施加全镜面统一
的等效粗糙度，额外随机曲率、随机法向和随机反射率均设为零。这种拆分方式
不需要向仰角序列 CSV 增加全零占位列。逐镜字段和统一参数的回退规则见
`docs/facet_parameter_precedence_zh.md`。

## 与原有配置的关系

- `configs/official_tests/perfect_*`：继续用于理想光学基准；
- `configs/official_tests/deformation_*`：继续用于只含支架形变的基准；
- `configs/errors/full_response_1229.cfg`：继续用于覆盖全部人工误差代码路径的
  回归测试，不代表当前实测参数；
- `configs/examples/lactroot_only.cfg`：改为当前实测约束参数，供正式批量模拟使用。
