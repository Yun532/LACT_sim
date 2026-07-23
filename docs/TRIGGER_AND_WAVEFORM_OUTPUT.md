# Trigger and waveform output

Camera triggers use a sliding p.e. coincidence window. Configure it with
`trigger.camera_coincidence_window_ns` (20 ns by default). Array triggers
require `trigger.array_multiplicity` camera triggers. Raw inter-telescope
timing is disabled by default with `trigger.array_coincidence_window_ns=0`
because geometric propagation-delay compensation is not yet implemented.
Positive array windows are available for controlled timing studies, but use
the uncorrected trigger times and should not be enabled for production data.

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

HDF5 waveform output defaults to sparse COO storage. The `waveforms/samples`
dataset stores `image_index`, `pixel_id`, `time_bin`, `photon_count`, `pe`,
`cherenkov_pe`, and `nsb_pe`; axes and per-image reference times remain in the
`waveforms` group. Set `output.hdf5_waveform_storage=dense` only when an
explicit dense `(image, time, pixel)` cube is required. Dense storage can be
orders of magnitude larger for mostly empty Cherenkov-camera data.
