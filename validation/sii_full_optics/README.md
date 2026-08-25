# SII full-main validation

These saved outputs were produced from Git commit `52245f3` on 2026-08-25.

```bash
build_sii/run_optical_sim \
  configs/optics/lact2_measured_full_response_400nm.cfg

build_sii/run_camera_electronics \
  configs/examples/corsika_lact_pylast_root_only_measured_waveform.cfg \
  validation/sii_full_optics/two_telescope_primary_pe_main_preview.csv \
  validation/sii_full_optics/main_electronics_2us \
  -C electronics.n_pixels=1 \
  -C electronics.sampling.start_ns=0 \
  -C electronics.sampling.end_ns=2000 \
  -C nsb.enabled=false \
  -C electronics.nsb.enabled=false
```

The input CSV already contains independently generated star and NSB primary
p.e. streams, so standalone NSB generation is disabled to avoid double counting.
The current main electronics configuration includes the measured SPE template
and charge distribution but no measured additive-noise PSD or ADC quantization.
