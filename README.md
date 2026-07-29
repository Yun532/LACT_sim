# LACT 坐标系诊断网站

这是一个与 `main` 源码隔离的成果分支，只保存坐标诊断网站、绘图数据、说明文档和可复现脚本，不包含 LACT 主程序源码。

## 直接使用

下载或克隆本分支后，在浏览器中打开：

`docs/assets/lact-coordinate-system-3d.html`

页面包含四部分：

1. 全局、望远镜本地、镜片/遮挡及 LACT/pyLAST 相机坐标定义；
2. 上/下/西/东各偏 3° 的真实平行光光路、完整反射/遮挡诊断点，以及 LACT/pyLAST 两种相机重画；
3. 不同天顶角的支架形变、同角度光路及完整输出光斑；
4. CORSIKA event 1909 的阵列、芯位、光子发射点反演，以及由同一 ROOT 输出重画的 LACT/pyLAST 相机图。

## 文档

- [中文坐标系总说明](docs/coordinate_systems_zh.md)
- [可执行代码坐标审计（推荐）](docs/coordinate_code_audit_zh.md)
- [早期坐标转换审计（历史记录）](docs/coordinate_transform_audit_zh.md)
- [English coordinate notes](docs/coordinate_systems.md)

## 数据边界

`docs/assets/data/` 保存网页和校验使用的精简真实输出；`configs/` 只保存镜片、形变、遮挡、相机几何数据和坐标事例配置。网页不保存或嵌入 pyLAST 原始 PNG，而是按最新版 reader 与 EventVisualizer 的代码规则从真实输出重画。原始 `run_logs`、`results` 和 4 GB CORSIKA 文件不在本分支中。

## 重新生成与校验

```bash
python python/build_coordinate_diagnostics_html.py --output docs/assets/lact-coordinate-system-3d.html
node python/test_coordinate_diagnostics_runtime.js
```

生成器会核对完整输出点数、相机命中数、四组 3° 输入、坐标投影和 CORSIKA 发射点反演；运行时测试还核对 pyLAST reader/renderer 映射、四页默认南/北屏幕方向和页面交互。
