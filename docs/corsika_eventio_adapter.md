# CORSIKA/EventIO Adapter

LACT_sim now has a direct binary EventIO-to-PhotonCsv converter:

```text
CORSIKA/EventIO photon bunches -> PhotonCsv -> LACT_sim raytrace
```

This keeps the optical ray tracer independent from external EventIO libraries,
while avoiding fragile parsing of `read_iact` text output.

The converter is:

```text
tools/eventio_to_photon_csv.c
tools/build_eventio_converter.sh
```

It is a thin wrapper around hessioxxx/eventio APIs:

```text
read_tel_photons()
read_tel_photons3d()
```

so hessioxxx still owns compressed-file reading, EventIO block parsing, and
compact/long/3D photon-bunch decoding.

## Vendored hessioxxx

LACT_sim vendors the official hessioxxx source archive under:

```text
external/hessioxxx/
```

See:

```text
external/hessioxxx/README.md
```

for the official download URL, SHA-256, upstream ChangeLog date, and build
steps. The active source tree is:

```text
external/hessioxxx/source
```

The relevant CORSIKA photon bunch structures are:

```text
struct bunch:
  photons
  x,y       arrival position relative to telescope [cm]
  cx,cy     direction cosines in the CORSIKA horizontal detection frame
  ctime     arrival time [ns]
  zem       emission height above sea level [cm]
  lambda    wavelength [nm] or 0

struct bunch3d:
  photons
  x,y,z     arrival position relative to telescope [cm]
  cx,cy,cz  direction cosines in the CORSIKA horizontal detection frame
  ctime     arrival time [ns]
  dist      emission-to-arrival distance [cm]
  lambda    wavelength [nm] or 0
```

## Direct Binary Conversion

Build the converter:

```bash
./tools/build_eventio_converter.sh
```

Run it:

```bash
DYLD_LIBRARY_PATH=external/hessioxxx/source/lib \
  ./build_eventio/eventio_to_photon_csv \
  /path/to/input.zst \
  run_logs/root_bunch_adapter/photon_E500_eventio_direct.csv \
  --event-id-mode event_array100
```

For this validation file, `event_array100` writes:

```text
event_id = shower_event * 100 + array_id
```

which matches the `runid` convention in the provided ROOT file. For other
workflows, use the default `--event-id-mode event` unless an array block needs
to become part of the event key.

Supported filters:

```bash
--filter-event-id N
--filter-array-id N
--filter-telescope-id N
```

The converter writes standard LACT_sim `PhotonCsv`:

```text
x_m,y_m,z_m,dir_x,dir_y,dir_z,time_ns,wavelength_nm,weight,multiplicity,event_id,telescope_id
```

Unit mapping in the low-level reader:

```text
CORSIKA x,y,z [cm] -> PhotonCsv x_m,y_m,z_m [m]
CORSIKA ctime [ns] -> time_ns
CORSIKA lambda [nm] -> wavelength_nm
CORSIKA photons -> multiplicity
2D bunch dir_z = -sqrt(1 - cx^2 - cy^2)
3D bunch dir_z = cz
```

For optical tracing, the recommended `run_corsika_trace` setting is:

```ini
source.eventio_coordinate_frame=corsika_iact
```

In this mode x/y/z are treated as telescope-relative positions in the CORSIKA
horizontal detection frame, and cx/cy/cz are rotated into the LACT telescope
optical frame with `telescope.pointing_az_deg` and
`telescope.pointing_el_deg`. This is different from `telescope_local`, where no
axis rotation is applied.

See [coordinate_systems.md](coordinate_systems.md) for the full CORSIKA/EventIO
axis, azimuth, ARRANG, telescope optical-frame, and camera-plane convention.

## One-Command zst to Images

`run_optical_sim` can now read EventIO directly when LACT_sim is built with the
project-local hessioxxx library:

```ini
source.mode=EventIO
source.eventio_path=
source.event_id_mode=event_array100
source.eventio_coordinate_frame=corsika_iact
source.use_eventio_telescope_position=true
source.filter_event_id=100
source.filter_telescope_id=0
```

`source.use_eventio_telescope_position` defaults to `true` in EventIO mode.
The selected telescope position is taken from the CORSIKA `IO_TYPE_MC_TELPOS`
table using `source.filter_telescope_id`. Set it to `false` only when you want
to force `telescope.position_m` from the LACT_sim config.

Every EventIO run prints a metadata block to the log:

```text
[EventIO metadata]
  input_lines
  tel[...]
  selected_shower_event
  energy_gev / direction / core_position_m
  selected_array_id
  array[...] offsets
```

`run_corsika_trace` uses a compact log format for full-file optical tracing.
It prints the same component-oriented configuration sections as the synthetic
optical tests:

```text
[Telescope]
[Mirror]
[Source]
[Output plane]
[Camera]
[SiPM]
[Electronics]
[Efficiency]
[Errors]
[Model]
```

The CORSIKA input card is written to the log for provenance, but local-path
lines starting with `TELFIL` or `DIRECT` are hidden. The log records their
count as `input_hidden_path_lines`.

## Full EventIO Optical Trace Test

The clean CORSIKA/EventIO workflow has two separated stages:

1. C++ optical trace stage: read one CORSIKA/EventIO file and trace every
   `(event_id, telescope_id)` stream in one run.
2. Plot stage: read only the C++ trace output, then plot a selected event and
   telescope.

Whiteboard trace:

```bash
DYLD_LIBRARY_PATH=$PWD/external/hessioxxx/source/lib \
build/run_corsika_trace configs/official_tests/corsika_whiteboard.cfg \
  /path/to/photon_E500_th0_run000001.zst \
  2>&1 | tee run_logs/official_tests/corsika/whiteboard_run.log
```

Dense HDF5 pixel-camera trace:

```bash
DYLD_LIBRARY_PATH=$PWD/external/hessioxxx/source/lib \
build/run_corsika_trace configs/official_tests/corsika_new_camera.cfg \
  /path/to/photon_E500_th0_run000001.zst \
  2>&1 | tee run_logs/official_tests/corsika/camera_run.log
```

Dense HDF5 pixel-camera trace with constant-rate NSB and simple multiplicity
trigger:

```bash
DYLD_LIBRARY_PATH=$PWD/external/hessioxxx/source/lib \
build/run_corsika_trace configs/official_tests/corsika_nsb_trigger_camera.cfg \
  /path/to/photon_E500_th0_run000001.zst \
  2>&1 | tee run_logs/official_tests/corsika/camera_nsb_trigger_run.log
```

These write:

```text
run_logs/official_tests/corsika/whiteboard_hits.csv
run_logs/official_tests/corsika/whiteboard_summary.csv
run_logs/official_tests/corsika/camera_dense.h5
run_logs/official_tests/corsika/camera_nsb_trigger_dense.h5
```

Plot one dense camera image by CORSIKA shower-event order:

```bash
python3 python/plot_hdf5_camera.py \
  run_logs/official_tests/corsika/camera_dense.h5 \
  --shower-event-number 1 \
  --array-id 0 \
  --telescope-id 7 \
  --quantity pe \
  --output run_logs/official_tests/corsika/camera_shower1_array0_tel7_pe.png
```

Plot one whiteboard event for one telescope:

```bash
python3 python/plot_corsika_trace_output.py \
  run_logs/official_tests/corsika/whiteboard_hits.csv \
  --event-id 100 \
  --telescope-id 0
```

CORSIKA input-card lines and telescope positions are printed in the trace log.
Extra metadata CSV files are intentionally not written by default.

The C++ trace stage reads EventIO metadata first, so the log reports the total
number of shower events and telescopes before photon tracing starts. Photon
bunches are then streamed directly through ray tracing and camera accumulation;
they are not preloaded as one full-file vector. During streaming, the log prints
`event_start`, `event_done`, `stream_progress`, and `stream_done` lines. This is
the recommended path for large server-side CORSIKA files.

Each trace log also contains a `Per-event summary` and a `Per-stream summary`.
A stream is one `(event_id, telescope_id)` pair. In the default
`event_array100` mode, `event_id = shower_event * 100 + array_id`, so the
single `event_id` is enough to select a simulated event stream later. The
summary includes
input bunches, input photon multiplicity, output-plane hits, weighted signal,
and time mean/RMS. Pixel-camera runs additionally include camera hits,
accepted hits, and the number of unique hit pixels. These fields are intended
as the lightweight handoff point for the future trigger stage.

Pixel-camera output is intentionally compact:

```text
event_id,telescope_id,pixel_id,photon_count,pe,signal,time_mean_ns,time_rms_ns
```

`pe` is evaluated after the light collector by multiplying the collected photon
weight by `electronics.pe_conversion`. The conversion can be a constant or a
two-column wavelength table. `signal` is currently kept as a backward-compatible
alias for `pe`.

## Diagnostic Text Bridge

If `read_iact` is built, it can print bunches as text. Then:

```bash
MAX_PRINT_ARRAY=1000000 /path/to/read_iact -n 1000000 input.corsika.gz \
  > run_logs/eventio/read_iact.txt
```

Convert that text to LACT_sim `PhotonCsv`:

```bash
python3 python/read_iact_text_to_photon_csv.py \
  run_logs/eventio/read_iact.txt \
  run_logs/eventio/photons_from_read_iact.csv \
  --default-event-id 1
```

Then run as usual with:

```ini
source.mode=PhotonCsv
source.csv_path=run_logs/eventio/photons_from_read_iact.csv
source.eventio_coordinate_frame=corsika_iact
```

This parser supports printed normal, compact, and 3D bunch formats. It is now
mainly for validation and debugging; the binary converter above is preferred.

## Current Local Probe

The vendored `external/hessioxxx/source` was successfully built with the
traditional Makefile, producing:

```text
external/hessioxxx/source/bin/read_iact
external/hessioxxx/source/lib/libhessio.dylib
```

A Linux `read_iact` binary from another machine cannot run on macOS. Build the
vendored hessioxxx locally instead.

A real CORSIKA IACT probe was run with:

```text
/path/to/run2.corsika.gz
```

and recorded in:

```text
run_logs/eventio_text_adapter_smoke/README.md
run_logs/eventio_text_adapter_smoke/run2_photons_from_read_iact_n5.csv
```

The converted output contains 64274 photon bunch rows from 60 events and 229
telescopes.

## Future In-Process Source

The next target is to make `run_optical_sim` accept:

```ini
source.mode=EventIO
source.eventio_path=...
source.filter_event_id=...
source.filter_telescope_id=...
```

For now, `source.mode=EventIO` intentionally still fails with a clear message.
The validated path is direct binary EventIO conversion to `PhotonCsv`, followed
by the existing `source.mode=PhotonCsv` raytrace.
