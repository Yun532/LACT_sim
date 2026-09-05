# SII说明与计算证据入口

当前说明拆分为两份，分别对应物理因果链和代码实现：

- [物理原理：从天体亮度到测量与重建](docs/SII_PHYSICS_ZH.md)
- [代码实现、数据接口与复现](docs/SII_IMPLEMENTATION_ZH.md)
- [完整科学验证Notebook](notebooks/sii_complete_waveform_report.ipynb)
- [当前计算结果汇总](validation/sii_science/summary.json)

当前采用main的仪器输入，缺失电子学项保持0且关闭。短波形是实际随机生成的光电子/SPE记录；长曝光是由波形标定误差生成的统计模拟，不是实测恒星数据。

原理、实现和当前数值以以上入口为准。旧的解析扫描及早期GLS Notebook仅保留为历史参考，不能混用其结果。
