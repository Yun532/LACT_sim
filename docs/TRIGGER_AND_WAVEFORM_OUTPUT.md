# Trigger and waveform output

Camera triggers use a sliding p.e. coincidence window. Configure it with
`trigger.camera_coincidence_window_ns` (20 ns by default). Array triggers
require `trigger.array_multiplicity` camera triggers. Raw inter-telescope
timing is disabled by default with `trigger.array_coincidence_window_ns=0`
for backward compatibility and for inputs without shower geometry. EventIO
runs can enable plane-wave compensation with:

```ini
trigger.array_coincidence_window_ns=50
trigger.array_time_correction=plane_wave
trigger.array_wavefront_speed_m_per_ns=0
```

The corrected time is
`raw_trigger_time + dot(telescope_position_nwu, viewing_direction_nwu) / speed`.
Here array `x` is North, `y` is West, `z` is Up, and shower azimuth increases
North-to-East. Missing shower direction or telescope position is an error;
the program never silently falls back to raw timing when `plane_wave` is
requested.
Zero speed selects the sim_telarray-compatible observation-level air speed;
a positive value explicitly overrides it.

`trigger.coincidence_window_ns` remains a compatibility alias: when either
new setting is omitted, an explicitly configured legacy value is inherited by
that setting. Prefer the two explicit settings in new configurations.

ROOT and HDF5 call the same camera and array trigger functions. When
`waveform.enabled=true` and `waveform.source=pe`, both formats evaluate the
camera trigger from the time-binned p.e. sequence. Without a p.e. time series,
the integrated image is treated as one time bin.

ROOT `output.lact_profile` controls serialization only. `image_pe` evaluates
the same time-binned trigger as `timeseries_pe` and `debug_full`, but does not
write the `waveforms` tree. Enabling deterministic time-binned NSB evaluates
the complete configured time range before `output.save_only_triggered` is
applied, so output profile selection cannot change accidental-trigger or
array-coincidence decisions.

`trigger_time_ns` is the raw local time of the earliest sliding window that
reaches `camera_multiplicity`. It is also repeated explicitly as
`trigger_first_time_ns`. The diagnostic
`trigger_max_multiplicity_time_ns` is the center of the earliest window having
the largest multiplicity; it can be later than the physical trigger. Array
coincidence always uses `trigger_time_ns`. ROOT observations and HDF5
`/trigger/telescope` also store `geometric_delay_ns` and
`coincidence_time_ns`, making the array timing decision directly auditable.

ROOT `image_time_mean_ns` and `image_time_rms_ns` are the signal-weighted
Cherenkov arrival-time moments accumulated before waveform binning. NSB is
included in `image_pe` and the serialized waveform, but does not dilute these
truth-signal timing moments. `image_time_peak_ns` remains the peak of the
time-binned total p.e. waveform; equal-height peaks select the earliest bin.

HDF5 waveform output defaults to sparse COO storage. The `waveforms/samples`
dataset stores `image_index`, `pixel_id`, `time_bin`, `photon_count`, `pe`,
`cherenkov_pe`, and `nsb_pe`; axes and per-image reference times remain in the
`waveforms` group. Set `output.hdf5_waveform_storage=dense` only when an
explicit dense `(image, time, pixel)` cube is required. Dense storage can be
orders of magnitude larger for mostly empty Cherenkov-camera data.
