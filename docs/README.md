# LACT_sim 文档目录

这里收录项目的长期使用说明、格式定义和设计文档。建议从用户手册开始；带日期的实施和验收记录统一放在 `archive/`，不再堆放在本目录顶层。

## 快速开始

- [中文用户手册](user_guide_zh.md)
- [English user guide](user_guide_en.md)
- [官方测试说明](official_tests.md)
- [最小 Photon CSV 示例](minimal_photon_csv.md)

## 输入与坐标

- [真实 LACT 三维坐标模型：镜面、相机、天空、输入轴与 u/v](assets/lact-coordinate-system-3d.html)
- [交互式坐标系总图：输入 x/y、天空指向、相机 u/v 与 pyLAST](assets/coordinate-system-explorer.html)
- [坐标系中文说明](coordinate_systems_zh.md)
- [Coordinate systems reference](coordinate_systems.md)
- [Photon CSV 格式](photon_csv_format.md)
- [CORSIKA/EventIO 适配器](corsika_eventio_adapter.md)
- [CORSIKA ROOT 适配器](corsika_root_adapter.md)
- [遮挡模型 CSV](obstruction_primitives_csv_zh.md)
- [逐镜参数优先级](facet_parameter_precedence_zh.md)

## 光学、相机与触发

- [Photon response modes](PHOTON_RESPONSE.md)
- [Trigger and waveform output](TRIGGER_AND_WAVEFORM_OUTPUT.md)
- [相机时间和波形](camera_timing_waveform_zh.md)
- [NSB 光谱模型](nsb_spectral_model_zh.md)
- [效率曲线验证](efficiency_validation.md)
- [LACT ROOT 中的 NSB 与触发流程](lact_root_nsb_trigger_flow_zh.md)

## 输出与 pylast

- [HDF5 输出格式](hdf5_output_format.md)
- [pylast event 数据层级](pylast_event_data_levels_zh.md)
- [LACT 到 pylast 的事件输出设计](pylast_lact_event_output_design_zh.md)
- [服务器 ROOT 输出与画图检查](server_root_output_check_zh.md)

## 架构与接口

- [程序中文总览](program_overview_zh.md)
- [模块接口](module_interfaces.md)

## 历史记录

- [2026-07 归档](archive/2026-07/)

归档用于保存已经落地的实施记录和验收结论。运行生成的 HDF5、ROOT、JSON、图片和日志不属于文档，不应提交到仓库；项目根目录的 `results/` 已由 `.gitignore` 忽略。
