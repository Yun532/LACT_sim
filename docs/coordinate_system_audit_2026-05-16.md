# Coordinate System Audit (2026-05-16)

This audit checks LACT_sim coordinate conventions against the current project
implementation, the vendored hessioxxx/sim_telarray headers, and public
CORSIKA IACT/ATMO documentation.

Initial audit notes were written before any fixes. Follow-up changes now make
`run_optical_sim` reject raw EventIO input and clarify the generic
`buildTelescopeFrame()` convention in the documentation.

## External / Upstream Reference Points

- hessio/sim_telarray telescope positions use the array frame:
  `x` toward North, `y` toward West, `z` upward.
  Local source: `external/hessioxxx/source/include/io_hess.h:690-692`.
  Public Doxygen reference:
  <https://www.mpi-hd.mpg.de/hfm/bernlohr/sim_telarray/Documentation/hessio_refman.pdf>
- hessio/CORSIKA photon bunches store `x,y` arrival position relative to the
  telescope in cm and `cx,cy` direction cosines. The 3D bunch stores
  `x,y,z` relative to the telescope and `cx,cy,cz`, with `cz` usually negative
  for downward photons.
  Local source: `external/hessioxxx/source/include/mc_tel.h:73-101`.
- hessio/sim_telarray azimuth is North-to-East in the NWU frame. The code uses
  a negative horizontal `y` component, for example
  `ay = sin(-azimuth) * cos(altitude)`.
  Local source: `external/hessioxxx/source/src/rec_tools.c:574-581`.
- CORSIKA IACT/ATMO public documentation describes photon bunch quantities as
  x/y in the CORSIKA detection plane and `cx/cy` as direction projections onto
  the X/Y axes.
  Public reference:
  <https://www.mpi-hd.mpg.de/hfm/~bernlohr/iact-atmo/iact_refman/group__iact__atmo__interface.html>

## Confirmed Consistent

### CORSIKA/EventIO metadata frame

Project documentation says CORSIKA/EventIO array metadata is NWU:

```text
+x = magnetic North
+y = West
+z = Up
```

This matches hessio/sim_telarray (`io_hess.h:690-692`). The HDF5 writer stores
explicit names `array_x_north_m`, `array_y_west_m`, and `array_z_up_m`, and its
coordinate metadata describes the same NWU frame.

Relevant local code:

- `apps/run_corsika_trace.cpp:695-713`
- `apps/run_corsika_trace.cpp:861-905`
- `docs/coordinate_systems.md:12-29`
- `docs/hdf5_output_format.md:190-220`

### EventIO photon bunch units and relativity

The EventIO reader converts hessio bunch positions from cm to meters and keeps
them as telescope-relative coordinates:

- 2D bunch: `pos = {b.x * 0.01, b.y * 0.01, 0.0}`
- 2D bunch downward direction: `dir_z = -sqrt(1 - cx^2 - cy^2)`
- 3D bunch: `pos = {b.x * 0.01, b.y * 0.01, b.z * 0.01}`, direction from
  `cx,cy,cz`

This matches `mc_tel.h:73-101` and the CORSIKA IACT/ATMO documentation.

Relevant local code:

- `src/io/EventIOPhotonSource.cpp:29-35`
- `src/io/EventIOPhotonSource.cpp:210-238`
- `tools/eventio_to_photon_csv.c:116-119`
- `tools/eventio_to_photon_csv.c:169-176`
- `tools/eventio_to_photon_csv.c:216-223`
- `docs/coordinate_systems.md:66-92`
- `docs/corsika_eventio_adapter.md:35-71`

### Dedicated CORSIKA IACT transform in `run_corsika_trace`

`run_corsika_trace` has a CORSIKA-specific frame:

```cpp
x_axis = (-sin(el) cos(az),  sin(el) sin(az), cos(el))
y_axis = (-sin(az),         -cos(az),         0)
z_axis = ( cos(el) cos(az), -cos(el) sin(az), sin(el))
```

This matches `docs/coordinate_systems.md:145-162` and the hessio/sim_telarray
North-to-East azimuth in NWU. In this mode, EventIO bunch positions are rotated
with dot products and telescope position is not subtracted, which is correct
because the bunch positions are already relative to each telescope.

Relevant local code:

- `apps/run_corsika_trace.cpp:107-147`
- `docs/coordinate_systems.md:139-167`
- `configs/examples/corsika_new_user_full.cfg:57-65`

### MC_TELOFF sign

hessio writes MC_TELOFF as detector-array offsets with respect to the shower
core. LACT_sim stores the event-level core position as the opposite vector:

```text
core_x_north_m = -MC_TELOFF.xoff
core_y_west_m  = -MC_TELOFF.yoff
```

This matches the vendored `io_simtel.c` comments and the current project docs.

Relevant local code:

- `external/hessioxxx/source/src/io_simtel.c:609-646`
- `apps/run_corsika_trace.cpp:216-229`
- `docs/coordinate_systems.md:190-238`
- `docs/hdf5_output_format.md:114-123`

### Array-layout plotting display conversion

The HDF5 array layout plot converts CORSIKA NWU to conventional display axes:

```text
display east  = -array_y_west_m
display north =  array_x_north_m
```

This matches the project coordinate docs.

Relevant local code:

- `python/plot_hdf5_array_layout.py:90-100`
- `python/plot_hdf5_array_layout.py:540-543`
- `docs/coordinate_systems.md:169-187`

## Findings / Risks

### Addressed P1: `run_optical_sim` EventIO path ignored `source.eventio_coordinate_frame`

`SourceRuntimeConfig` defaults EventIO input to:

```ini
source.eventio_coordinate_frame=corsika_iact
```

but `run_optical_sim` constructs `EventIOPhotonSource` directly and later only
applies the generic telescope frame when `source.local_telescope_frame` is true.
It does not call the CORSIKA-aware `transformEventIOBunchToTraceFrame()` used by
`run_corsika_trace`.

Effect:

- A config using `run_optical_sim` with `source.mode=EventIO` and
  `source.eventio_coordinate_frame=corsika_iact` can be misleading.
- The EventIO bunches are treated via `source.local_telescope_frame` behavior,
  not via the documented CORSIKA NWU-to-local basis.
- Because `source.local_telescope_frame` defaults to true, raw EventIO NWU
  bunches can be transformed as if they were already LACT telescope-local
  photons.

Relevant local code:

- `src/app/OpticalSimCommon.cpp:1395-1401`
- `apps/run_optical_sim.cpp:343-355`
- `apps/run_optical_sim.cpp:391-399`
- CORSIKA-aware path exists only in `apps/run_corsika_trace.cpp:130-160`.

Suggested fix direction:

- Implemented direction: `run_optical_sim` now rejects `source.mode=EventIO`.
  CORSIKA/EventIO input is documented as belonging to `run_corsika_trace`.

### Addressed P1: Generic `buildTelescopeFrame()` uses a +Y azimuth basis, while CORSIKA docs describe +Y as West

`buildTelescopeFrame()` currently uses:

```cpp
z_axis = (cos(el) cos(az), cos(el) sin(az), sin(el))
x_axis = (-sin(az), cos(az), 0)
```

This means `az=90 deg` points along global `+y`. In the documented
CORSIKA/sim_telarray NWU frame, `az=90 deg` should point East, i.e. global
`-y`.

This is not necessarily wrong for synthetic/local workflows, but the function
name and docs currently call this an "array" frame without distinguishing ENU-
like generic array coordinates from CORSIKA NWU. That makes downstream use easy
to misinterpret.

Relevant local code:

- `src/app/OpticalSimCommon.cpp:525-552`
- `python/config_io.py:449-462`
- Contrasting CORSIKA transform:
  `apps/run_corsika_trace.cpp:107-127`
- Contrasting docs:
  `docs/coordinate_systems.md:96-116`

Suggested fix direction:

- Implemented direction: keep `buildTelescopeFrame()` as the generic LACT
  frame where azimuth is measured from `+x` toward `+y`, and document that raw
  CORSIKA/EventIO uses the separate `run_corsika_trace` NWU transform.

### P2: Python visualization/helpers contain both generic and CORSIKA bases

`python/config_io.py:telescope_frame_from_config()` mirrors
`buildTelescopeFrame()` and uses `+y` for `az=90 deg`. Other tools use the
CORSIKA-specific sign:

- `python/root_photons_to_photon_csv.py:108-118`
- `python/plot_orientation.py:49-64`
- `python/plot_hdf5_array_layout.py:217-221`

This is consistent only if `config_io.py` is treated as generic-local/ENU-like
and the other helpers as CORSIKA/NWU. The current docs do not make this split
explicit enough.

Relevant local code:

- `python/config_io.py:449-462`
- `python/root_photons_to_photon_csv.py:108-118`
- `python/plot_orientation.py:49-64`

Suggested fix direction:

- Add explicit names such as `telescope_frame_generic_xy()` and
  `telescope_frame_corsika_nwu()`, or add a config-driven frame selector.
- Update plotting/tool docs to say which helper assumes which frame.

### P2: `docs/photon_csv_format.md` does not state the global-frame axis convention

The PhotonCsv doc says:

```ini
local_telescope_frame=false
```

means rows are already in the global array frame, but it does not define whether
that global frame is generic `+y` East, CORSIKA `+y` West, or something else.
Because `run_optical_sim` uses `buildTelescopeFrame()` for global/local
transforms, this ambiguity matters for non-local CSV files.

Relevant docs/code:

- `docs/photon_csv_format.md:36-50`
- `src/app/OpticalSimCommon.cpp:525-552`
- `apps/run_optical_sim.cpp:391-399`

Suggested fix direction:

- Add a short "Global PhotonCsv frame" section that either binds it to the
  generic `buildTelescopeFrame()` convention or introduces a selectable frame.
- Warn that CORSIKA/EventIO photon CSV rows should not be fed as generic global
  CSV unless they have been explicitly converted to the selected LACT frame.

### P3: CORSIKA template comments say `pointing_az_deg = azimuth` without restating NWU sign

The templates point users to `docs/coordinate_systems.md`, but
`configs/templates/corsika_camera.cfg` only says:

```text
pointing_az_deg = azimuth
```

This is true if the azimuth is the hessio/sim_telarray North-to-East value in
the CORSIKA NWU frame, but it can be misread as a generic ENU azimuth if copied
without the docs.

Relevant docs/configs:

- `configs/templates/corsika_camera.cfg:14-17`
- `configs/templates/minimal_corsika_camera.cfg:22-26`
- `README.md:208-216`

Suggested fix direction:

- In CORSIKA templates, say explicitly:
  `azimuth is hessio/sim_telarray N->E; az=90 means East = -CORSIKA y`.

## Recommended Review Order

1. Decide whether `buildTelescopeFrame()` should remain a generic array frame
   or become frame-selectable.
2. Keep raw EventIO/CORSIKA input on `run_corsika_trace`; `run_optical_sim`
   now rejects that mode.
3. After that decision, update docs and helper names so "generic array",
   "CORSIKA NWU", and "telescope-local" are never conflated.
4. Add regression tests for `az=90 deg` in both the generic and CORSIKA paths:
   one should prove the intended `+y` behavior, the other the intended `-y`
   behavior.
