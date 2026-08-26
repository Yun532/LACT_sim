# Crab NSB2 validation input

`crab_rate_cube_nsb2_0p2p0.npz` was produced with the official `nsb2` 0.2.0
full-sky model for LACT (29.3576667 N, 100.1387778 E, 4410 m), pointed at the
Crab at 2026-12-05T19:20:00 UTC.  It contains 12 wavelength bands from 300 to
900 nm and a 41 x 41 angular grid covering +/-4.5 degrees.  The cached cube
lets CI reproduce PhotonCsv and camera plots without downloading the Gaia and
ESO source catalogues again.

