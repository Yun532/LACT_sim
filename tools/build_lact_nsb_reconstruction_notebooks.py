#!/usr/bin/env python3
"""Build LACT event reconstruction notebooks from the standard LACT template."""

from __future__ import annotations

import argparse
import copy
import shutil
from pathlib import Path

import nbformat as nbf


DEFAULT_ROOT_FILE = "/tmp/lact_notebook_no_nsb_all.root"
DEFAULT_EVENT_ID = 1417


def code(source: str):
    return nbf.v4.new_code_cell(source.strip() + "\n")


def markdown(source: str):
    return nbf.v4.new_markdown_cell(source.strip() + "\n")


def replace_cell(notebook, index: int, cell):
    cell["id"] = notebook.cells[index].get("id", cell["id"])
    notebook.cells[index] = cell


def build(template, with_nsb: bool):
    notebook = copy.deepcopy(template)
    notebook["metadata"]["kernelspec"] = {
        "display_name": "Python 3",
        "language": "python",
        "name": "python3",
    }
    notebook["metadata"]["language_info"] = {
        "name": "python",
        "version": "3.11",
    }
    mode_name = "有 NSB" if with_nsb else "无 NSB"
    raw_image_level = "simulation_fake" if with_nsb else "simulation"
    poisson_parameters = (
        """
NSB_RATE_PE_PER_NS_PER_PIXEL = 0.070527457
INTEGRATION_WINDOW_NS = 32.0
POISSON_MEAN_PE = NSB_RATE_PE_PER_NS_PER_PIXEL * INTEGRATION_WINDOW_NS
"""
        if with_nsb
        else "POISSON_MEAN_PE = 0.0\n"
    )

    replace_cell(
        notebook,
        0,
        markdown(
            f"""
# LACT 单事件完整可视化与 Hillas 重建（{mode_name}）

本 notebook 沿用 LACT_sim 标准 notebook 的完整顺序：读取事件、阵列图、原始相机图、
图像清理与 Hillas 参数、方向/芯位重建、重建后相机图和 3D SDP。
两种模式只有 Cell 1 的 NSB 参数不同。
"""
        ),
    )
    replace_cell(
        notebook,
        1,
        code(
            f"""
# Cell 1：所有需要修改的输入和分析参数都集中在这里。
from pathlib import Path
import os

os.environ.setdefault("MPLCONFIGDIR", "/tmp/pylast_mplcache")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

INPUT_KIND = "lact_root"
INPUT_FILE = Path({DEFAULT_ROOT_FILE!r})
EVENT_ID = {DEFAULT_EVENT_ID}
EVENT_INDEX = 0
MAX_EVENTS = -1

PICTURE_THRESHOLD_PE = 10.0
BOUNDARY_THRESHOLD_PE = 5.0
MIN_PICTURE_NEIGHBORS = 2
TRIGGER_THRESHOLD_PE = 8.0
TRIGGER_PIXELS = 3
IMAGE_QUERY = "hillas_intensity > 100 && leakage_intensity_width_2 < 0.3"
{poisson_parameters.rstrip()}

RAW_IMAGE_LEVEL = {raw_image_level!r}
CLEAN_IMAGE_LEVEL = "simulation_fake_clean"

print("mode:", {mode_name!r})
print("input kind:", INPUT_KIND)
print("input file:", INPUT_FILE)
print("exists:", INPUT_FILE.exists())
"""
        ),
    )
    replace_cell(
        notebook,
        2,
        code(
            """
import json
import time
import numpy as np

try:
    import plotly.io as pio
    pio.renderers.default = "notebook_connected"
except Exception:
    pio = None

from pylast.io import LactEventSource, SimtelEventSource, RootEventSource
from pylast.calib import Calibrator
from pylast.image import ImageProcessor
from pylast.reco import ShowerProcessor
from pylast.visualize import (
    EventVisualizer,
    hillas_parameter_rows,
    plot_clean_images,
    plot_event_cameras,
    plot_gathered_images,
    plot_event_cores,
    plot_event_trigger_timing,
    plot_event_sdp_planes,
    plot_event_sdp_planes_3d,
    plot_event_sdp_planes_3d_interactive,
    plot_raw_images,
    reconstruction_summary,
)

CALIBRATOR_CONFIG = {
    "Calibrator": {
        "image_extractor_type": "LocalPeakExtractor",
        "LocalPeakExtractor": {
            "window_width": 7,
            "window_shift": 3,
            "apply_correction": False,
        },
    }
}
IMAGE_PROCESSOR_CONFIG = {
    "poisson_noise": POISSON_MEAN_PE,
    "trigger_pe": TRIGGER_THRESHOLD_PE,
    "trigger_pixels": TRIGGER_PIXELS,
    "image_cleaner_type": "Tailcuts_cleaner",
    "Tailcuts_cleaner": {
        "picture_thresh": PICTURE_THRESHOLD_PE,
        "boundary_thresh": BOUNDARY_THRESHOLD_PE,
        "keep_isolated_pixels": False,
        "min_number_picture_neighbors": MIN_PICTURE_NEIGHBORS,
    },
}
SHOWER_PROCESS_CONFIG = {
    "ShowerProcessor": {
        "GeometryReconstructionTypes": ["HillasReconstructor"],
        "HillasReconstructor": {
            "use_fake_hillas": True,
            "ImageQuery": IMAGE_QUERY,
        },
    }
}


def step(name, function):
    print(f"BEGIN {name}", flush=True)
    start = time.perf_counter()
    result = function()
    print(f"END {name}: {time.perf_counter() - start:.3f}s", flush=True)
    return result


def read_one_event(kind, filename, event_id=None, event_index=0, max_events=-1):
    filename = str(filename)
    if kind == "lact_root":
        source_data = LactEventSource(filename, max_events=max_events)
        if event_id is None:
            event = source_data[event_index]
        else:
            event = next(
                source_data[index]
                for index in range(len(source_data))
                if int(source_data[index].event_id) == int(event_id)
            )
    elif kind == "simtel":
        source_data = SimtelEventSource(filename, max_events=max(max_events, event_index + 1))
        event = list(source_data)[event_index]
    elif kind == "root_event":
        source_data = RootEventSource(filename, max_events=max_events)
        event = source_data[event_index]
    else:
        raise ValueError("INPUT_KIND must be lact_root, simtel, or root_event")
    return source_data, event, EventVisualizer(source_data)


def build_processors(source_data):
    calibrator = Calibrator(source_data.subarray, json.dumps(CALIBRATOR_CONFIG))
    image_processor = ImageProcessor(source_data.subarray, json.dumps(IMAGE_PROCESSOR_CONFIG))
    shower_processor = ShowerProcessor(source_data.subarray, json.dumps(SHOWER_PROCESS_CONFIG))
    return calibrator, image_processor, shower_processor


def prepare_raw_image(event, calibrator):
    if getattr(event, "dl0", None) is None:
        step("提取 raw image", lambda: calibrator(event))


def print_event_summary(event):
    shower = event.simulation.shower
    print("event_id:", event.event_id)
    print("run_id:", event.run_id)
    print("energy [TeV]:", float(shower.energy))
    print("true zenith [deg]:", 90.0 - np.degrees(float(shower.alt)))
    print("true azimuth [deg]:", np.degrees(float(shower.az)))
    print("true core [m]:", float(shower.core_x), float(shower.core_y))
"""
        ),
    )
    replace_cell(notebook, 3, markdown("## 1. 读取事件"))
    replace_cell(
        notebook,
        4,
        code(
            """
source_data, event, visualizer = step(
    "读取事件",
    lambda: read_one_event(INPUT_KIND, INPUT_FILE, EVENT_ID, EVENT_INDEX, MAX_EVENTS),
)
calibrator, image_processor, shower_processor = build_processors(source_data)
print("subarray telescopes:", len(source_data.subarray.tels))
"""
        ),
    )
    replace_cell(
        notebook,
        6,
        code(
            """
prepare_raw_image(event, calibrator)
print_event_summary(event)

# ImageProcessor 同时生成 simulation.fake_image、清理掩膜和 Hillas 参数。
step("图像处理 + Hillas 参数", lambda: image_processor(event))
print("Poisson NSB mean [PE/pixel]:", POISSON_MEAN_PE)
"""
        ),
    )
    replace_cell(
        notebook,
        7,
        code(
            """
# 重新选择事件时，只需修改 Cell 1 后从这里以上重新运行。
# source_data, event, visualizer = read_one_event(
#     INPUT_KIND, INPUT_FILE, EVENT_ID, EVENT_INDEX, MAX_EVENTS
# )
# calibrator, image_processor, shower_processor = build_processors(source_data)
# prepare_raw_image(event, calibrator)
# image_processor(event)
"""
        ),
    )
    replace_cell(
        notebook,
        8,
        markdown(
            """
## 2. 望远镜分布、触发时间、芯位与 SDP 投影

保持 LACT_sim 标准 notebook 的阵列级检查。pyLAST 后加的积分 NSB 没有独立时间轴，
所以触发时间图仍读取 LACT_sim ROOT 中原始保存的触发时间。
"""
        ),
    )
    replace_cell(
        notebook,
        9,
        code(
            """
core_result = plot_event_cores(
    event,
    visualizer=visualizer,
    image_level=RAW_IMAGE_LEVEL,
    include_non_triggered=False,
)
"""
        ),
    )
    replace_cell(
        notebook,
        11,
        code(
            """
sdp_result = plot_event_sdp_planes(
    event,
    visualizer=visualizer,
    image_level=RAW_IMAGE_LEVEL,
    include_non_triggered=False,
)
"""
        ),
    )
    replace_cell(
        notebook,
        12,
        code(
            """
# 与 LACT_sim 标准 notebook 一样，直接使用 pyLAST 原生 EventVisualizer。
_ = plot_event_cameras(
    event,
    visualizer=visualizer,
    image_level="simulation",
    include_non_triggered=False,
)
"""
        ),
    )
    replace_cell(
        notebook,
        13,
        markdown(
            f"""
## 3. Raw image（{mode_name}）

这里画清理前的 `{raw_image_level}`。有 NSB 模式已经加入泊松 NSB，并减去平均基线；
无 NSB 模式直接显示纯切伦科夫 `simulation.true_image`。
"""
        ),
    )
    replace_cell(
        notebook,
        14,
        code(
            """
raw_result = plot_raw_images(
    event,
    visualizer=visualizer,
    image_level=RAW_IMAGE_LEVEL,
    include_non_triggered=False,
)
"""
        ),
    )
    replace_cell(
        notebook,
        15,
        code(
            """
_ = plot_gathered_images(
    event,
    visualizer=visualizer,
    image_type=None,
    image_level=RAW_IMAGE_LEVEL,
    show_hillas=False,
    include_non_triggered=False,
)
"""
        ),
    )
    replace_cell(
        notebook,
        16,
        markdown(
            """
## 4. Clean image 与 Hillas 参数

显示同一次 `ImageProcessor` 产生的清理掩膜、清理后图像和 Hillas 椭圆。
"""
        ),
    )
    replace_cell(
        notebook,
        17,
        code(
            """
rows = hillas_parameter_rows(event, image_level=CLEAN_IMAGE_LEVEL)
print("有 Hillas 参数的望远镜:", [row["tel_id"] for row in rows])
for row in rows:
    print(
        f"Tel {row['tel_id']:2d}: intensity={row['intensity']:.2f}, "
        f"length={row['length_rad']:.5g} rad, width={row['width_rad']:.5g} rad, "
        f"psi={np.degrees(row['psi_rad']):.2f} deg"
    )
"""
        ),
    )
    replace_cell(
        notebook,
        18,
        code(
            """
clean_result = plot_clean_images(
    event,
    visualizer=visualizer,
    image_level=CLEAN_IMAGE_LEVEL,
    show_hillas=True,
    include_non_triggered=False,
    ideal=True,
)
"""
        ),
    )
    replace_cell(notebook, 19, markdown("## 5. 方向与芯位重建"))
    replace_cell(
        notebook,
        20,
        code(
            """
step("方向/芯位重建", lambda: shower_processor(event))
summary = reconstruction_summary(event, "HillasReconstructor")
for key, value in summary.items():
    if isinstance(value, float):
        print(f"{key}: {value:.6g}")
    else:
        print(f"{key}: {value}")
assert summary["is_valid"], "该事件没有得到有效 Hillas 立体重建"
"""
        ),
    )
    replace_cell(
        notebook,
        21,
        code(
            """
_ = plot_clean_images(
    event,
    visualizer=visualizer,
    image_level=CLEAN_IMAGE_LEVEL,
    show_hillas=True,
    include_non_triggered=False,
    ideal=True,
    reco=True,
)
"""
        ),
    )
    replace_cell(
        notebook,
        22,
        code(
            """
_ = plot_gathered_images(
    event,
    visualizer=visualizer,
    image_type=None,
    image_level=CLEAN_IMAGE_LEVEL,
    show_hillas=True,
    include_non_triggered=False,
    ideal=True,
    reco=True,
)
"""
        ),
    )
    replace_cell(notebook, 23, markdown("## 6. 3D SDP 平面"))
    replace_cell(
        notebook,
        24,
        code(
            """
sdp_3d_interactive_figure = None
if pio is None:
    print("未安装 Plotly，跳过交互式 3D SDP；下一 Cell 仍会绘制静态版本。")
else:
    sdp_3d_interactive_figure = plot_event_sdp_planes_3d_interactive(
        event,
        visualizer=visualizer,
        image_level=CLEAN_IMAGE_LEVEL,
        include_non_triggered=False,
        z_max=1200.0,
        show_reco=True,
    )
    sdp_3d_interactive_figure
"""
        ),
    )
    replace_cell(
        notebook,
        25,
        code(
            """
sdp_3d_result = plot_event_sdp_planes_3d(
    event,
    visualizer=visualizer,
    image_level=CLEAN_IMAGE_LEVEL,
    include_non_triggered=False,
    z_max=1200.0,
    show_reco=True,
    show=True,
)
"""
        ),
    )
    replace_cell(notebook, 26, markdown("## 7. 可选：保存图片"))
    replace_cell(
        notebook,
        27,
        code(
            """
# OUTPUT_DIR = Path(INPUT_FILE).parent / "pylast_visualize"
# OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
# for name, result in {
#     "core": core_result,
#     "sdp": sdp_result,
#     "raw": raw_result,
#     "clean": clean_result,
#     "sdp_3d": sdp_3d_result,
# }.items():
#     figure = result.get("figure")
#     if figure is not None:
#         figure.savefig(
#             OUTPUT_DIR / f"event_{event.event_id}_{name}.png",
#             dpi=200,
#             bbox_inches="tight",
#         )
# if sdp_3d_interactive_figure is not None:
#     sdp_3d_interactive_figure.write_html(
#         OUTPUT_DIR / f"event_{event.event_id}_sdp_3d_interactive.html",
#         include_plotlyjs="cdn",
#         full_html=True,
#     )
"""
        ),
    )
    replace_cell(
        notebook,
        28,
        markdown(
            """
## 项目位置

这两本 notebook 在 LACT_sim 和 pyLAST 的 `notebooks/` 中保持完全相同。
修改 Cell 1 的 `INPUT_FILE` 和 `EVENT_ID` 后即可用于其他 LACT ROOT 文件。
"""
        ),
    )

    for cell in notebook.cells:
        if cell.cell_type == "code":
            cell.execution_count = None
            cell.outputs = []
    return notebook


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mirror-notebook-dir", type=Path)
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    notebook_dir = repo_root / "notebooks"
    template_path = notebook_dir / "lact_root_to_pylast_visualize.ipynb"
    template = nbf.read(template_path, as_version=4)
    outputs = {
        "lact_event_reconstruction_no_nsb.ipynb": build(template, False),
        "lact_event_reconstruction_with_nsb.ipynb": build(template, True),
    }
    for filename, notebook in outputs.items():
        output_path = notebook_dir / filename
        nbf.write(notebook, output_path)
        print(output_path)
        if args.mirror_notebook_dir is not None:
            args.mirror_notebook_dir.mkdir(parents=True, exist_ok=True)
            mirror_path = args.mirror_notebook_dir / filename
            shutil.copy2(output_path, mirror_path)
            print(mirror_path)


if __name__ == "__main__":
    main()
