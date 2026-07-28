# LACT 坐标系诊断网站

这是一个与 `main` 源码隔离的成果分支，只保存坐标诊断网站、绘图数据、说明文档和可复现脚本，不包含 LACT 主程序源码。

## 直接使用

下载或克隆本分支后，在浏览器中打开：

`docs/assets/lact-coordinate-system-3d.html`

页面包含四部分：

1. 全局、望远镜本地、镜片/遮挡及 LACT/pyLAST 相机坐标定义；
2. 四方向真实平行光光路、完整反射/遮挡诊断点和白板光斑；
3. 不同天顶角的支架形变、同角度光路及完整输出光斑；
4. CORSIKA event 1909 的阵列、芯位、光子发射点反演、LACT ROOT 原值相机和最新 pyLAST 原生相机图。

## 文档

- [中文坐标系总说明](docs/coordinate_systems_zh.md)
- [可执行代码坐标审计（推荐）](docs/coordinate_code_audit_zh.md)
- [早期坐标转换审计（历史记录）](docs/coordinate_transform_audit_zh.md)
- [English coordinate notes](docs/coordinate_systems.md)

## 数据边界

`docs/assets/data/` 保存网页和校验使用的精简真实输出、pyLAST 原生 PNG 及来源审计；`configs/` 只保存镜片、形变、遮挡、相机几何数据和坐标事例配置。原始 `run_logs`、`results` 和 4 GB CORSIKA 文件不在本分支中。

## 重新生成与校验

```bash
python python/build_coordinate_diagnostics_html.py --output docs/assets/lact-coordinate-system-3d.html
node python/test_coordinate_diagnostics_runtime.js
```

生成器会核对完整输出点数、相机命中数、坐标投影、CORSIKA 发射点反演、pyLAST/ROOT 坐标与图像 SHA-256，以及页面运行时交互。
