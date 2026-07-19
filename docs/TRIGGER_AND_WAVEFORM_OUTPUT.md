# Trigger and waveform output

Camera triggers use a sliding p.e. coincidence window. Configure it with
`trigger.camera_coincidence_window_ns`. Array triggers then require
`trigger.array_multiplicity` camera triggers inside
`trigger.array_coincidence_window_ns`.

`trigger.coincidence_window_ns` remains a compatibility alias: when either
new setting is omitted, that setting inherits the legacy value.

ROOT and HDF5 call the same camera and array trigger functions. When
`waveform.enabled=true` and `waveform.source=pe`, both formats evaluate the
camera trigger from the time-binned p.e. sequence. Without a p.e. time series,
the integrated image is treated as one time bin.

HDF5 waveform output defaults to sparse COO storage. The `waveforms/samples`
dataset stores `image_index`, `pixel_id`, `time_bin`, `photon_count`, `pe`,
`cherenkov_pe`, and `nsb_pe`; axes and per-image reference times remain in the
`waveforms` group. Set `output.hdf5_waveform_storage=dense` only when an
explicit dense `(image, time, pixel)` cube is required. Dense storage can be
orders of magnitude larger for mostly empty Cherenkov-camera data.
