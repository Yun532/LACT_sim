#!/usr/bin/env python3
"""Plot LACT_sim proxy waveform camera frames from HDF5.

The waveform datasets are optional debug/diagnostic outputs. They represent
time-binned photon count or p.e. at the camera/collector output, not a real
electronics waveform.
"""

import argparse
from pathlib import Path

import h5py
import matplotlib.pyplot as plt

from plot_hdf5_camera import (
    draw_camera,
    find_images,
    resolve_event_id,
    telescope_pointing_for_image,
)
from plot_orientation import basis_from_corsika_pointing


def output_frame_path(output_dir, image, quantity, frame_index):
    return output_dir / (
        f"event_{int(image['event_id'])}_tel_{int(image['telescope_id']) + 1:02d}_"
        f"{quantity}_bin_{frame_index:04d}.png"
    )


def output_gif_path(gif_arg, output_dir, image, quantity, multiple):
    if not gif_arg:
        return None
    gif_path = Path(gif_arg)
    if not multiple:
        return gif_path
    gif_dir = gif_path if gif_path.suffix == "" else gif_path.with_suffix("")
    gif_dir.mkdir(parents=True, exist_ok=True)
    return gif_dir / (
        f"event_{int(image['event_id'])}_tel_{int(image['telescope_id']) + 1:02d}_{quantity}.gif"
    )


def main():
    parser = argparse.ArgumentParser(
        description="Plot time-binned proxy waveform frames from LACT_sim HDF5."
    )
    parser.add_argument("h5", help="HDF5 file written by run_corsika_trace")
    parser.add_argument("--event-id", type=int, default=None)
    parser.add_argument("--shower-event-number", "--event-number", type=int, default=None)
    parser.add_argument("--shower-event-id", type=int, default=None)
    parser.add_argument("--array-id", type=int, default=0)
    parser.add_argument("--telescope-id", type=int, default=None)
    parser.add_argument("--image-index", type=int, default=None)
    parser.add_argument(
        "--quantity",
        choices=("photon_count", "cherenkov_pe", "nsb_pe", "pe"),
        default="pe",
    )
    parser.add_argument("--output-dir", default="waveform_frames")
    parser.add_argument("--gif", default=None, help="optional animated GIF path")
    parser.add_argument("--stride", type=int, default=1, help="plot every Nth time bin")
    parser.add_argument("--dpi", type=int, default=180)
    parser.add_argument("--raw-camera-xy", action="store_true")
    args = parser.parse_args()

    if args.stride < 1:
        raise SystemExit("--stride must be >= 1")

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    with h5py.File(args.h5, "r") as h5:
        if "waveforms" not in h5:
            raise SystemExit("This HDF5 file has no /waveforms group.")
        dataset_name = f"waveforms/{args.quantity}"
        if dataset_name not in h5:
            available = ", ".join(k for k in h5["waveforms"].keys())
            raise SystemExit(f"{dataset_name} is missing. Available waveform datasets: {available}")

        event_id = resolve_event_id(
            h5,
            event_id=args.event_id,
            shower_event_id=args.shower_event_id,
            shower_event_number=args.shower_event_number,
            array_id=args.array_id,
        )
        index = h5["images/index"][:]
        images = find_images(index, event_id, args.telescope_id, args.image_index)
        camera = h5["camera/pixels"][:]
        pixel_axis = h5["waveforms/pixel_id_axis"][:]
        time_centers = h5["waveforms/time_centers_ns"][:]
        waveform = h5[dataset_name]

        all_saved = []
        gif_outputs = []
        multiple = len(images) > 1
        for image in images:
            image_index = int(image["image_index"])
            if image_index < 0 or image_index >= waveform.shape[0]:
                raise SystemExit(f"image_index={image_index} is outside waveform dataset shape.")

            display_basis = None
            if not args.raw_camera_xy:
                pointing = telescope_pointing_for_image(h5, image)
                if pointing is not None:
                    display_x, display_y, oriented = basis_from_corsika_pointing(*pointing)
                    if oriented:
                        display_basis = (display_x, display_y)

            image_frame_dir = out_dir if not multiple else out_dir / (
                f"event_{int(image['event_id'])}_tel_{int(image['telescope_id']) + 1:02d}_{args.quantity}"
            )
            image_frame_dir.mkdir(parents=True, exist_ok=True)
            saved = []
            for bin_index in range(0, waveform.shape[1], args.stride):
                values = waveform[image_index, bin_index, :]
                values_by_pixel = {int(pid): float(v) for pid, v in zip(pixel_axis, values)}
                fig = draw_camera(
                    camera,
                    image,
                    values_by_pixel,
                    args.quantity,
                    args.dpi,
                    display_basis,
                )
                ax = fig.axes[0]
                ax.set_title(
                    f"event {int(image['event_id'])}, telescope {int(image['telescope_id']) + 1}, "
                    f"t = {float(time_centers[bin_index]):.2f} ns"
                )
                out = output_frame_path(image_frame_dir, image, args.quantity, bin_index)
                fig.savefig(out, bbox_inches="tight")
                plt.close(fig)
                saved.append(out)
            all_saved.extend(saved)

            gif_path = output_gif_path(args.gif, out_dir, image, args.quantity, multiple)
            if gif_path:
                try:
                    from PIL import Image
                except ImportError as exc:
                    raise SystemExit(
                        "Pillow is required for --gif. Frames were written; install pillow to make GIF."
                    ) from exc
                frames = [Image.open(path) for path in saved]
                if frames:
                    if gif_path.parent:
                        gif_path.parent.mkdir(parents=True, exist_ok=True)
                    frames[0].save(
                        gif_path,
                        save_all=True,
                        append_images=frames[1:],
                        duration=120,
                        loop=0,
                    )
                    gif_outputs.append(gif_path)

    if gif_outputs:
        print(f"Saved {len(all_saved)} frames to {out_dir} and {len(gif_outputs)} GIF files")
    else:
        print(f"Saved {len(all_saved)} waveform frames to {out_dir}")


if __name__ == "__main__":
    main()
