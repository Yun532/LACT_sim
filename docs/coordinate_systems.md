# LACT_sim Coordinate Systems

This document is the reference convention for CORSIKA/EventIO input, telescope
pointing, optical ray tracing, camera images, and array-layout plots.

## Short Answer

For CORSIKA/EventIO runs, LACT_sim keeps the CORSIKA horizontal coordinate
system as the canonical metadata frame:

```text
CORSIKA IACT frame = NWU
+x : magnetic North
+y : West
+z : Up
```

This is not the usual ENU frame. LACT_sim does not currently apply a magnetic
declination correction from magnetic north to geographic true north. Therefore
all CORSIKA/EventIO telescope positions, core positions, and CORSIKA provenance
tables use **magnetic-North-West-Up** coordinates.

For plots that should look like a conventional map, LACT_sim converts only at
the plotting layer:

```text
plot North =  CORSIKA x
plot East  = -CORSIKA y
plot Up    =  CORSIKA z
```

So the program is internally consistent with CORSIKA. If later you need true
geographic NEU/ENU, add one explicit rotation by the site magnetic declination
before comparing to GPS/geographic coordinates.

## Coordinate Frames

### 1. CORSIKA / EventIO Array Frame

CORSIKA IACT photon bunches and EventIO telescope tables are interpreted in the
CORSIKA horizontal frame:

```text
x_N_m : magnetic North positive
y_W_m : West positive
z_U_m : Up positive
```

In HDF5 these are stored with explicit names:

```text
/telescopes/table:
  array_x_north_m
  array_y_west_m
  array_z_up_m

/events/corsika:
  core_x_north_m
  core_y_west_m
```

The legacy aliases `x_m,y_m,z_m` in `/telescopes/table` are kept for old
scripts, but new code should use the explicit `array_x_north_m`,
`array_y_west_m`, and `array_z_up_m` names.

### 2. EventIO Photon Bunch Frame

The hessio/CORSIKA bunch definitions are:

```text
struct bunch:
  x, y       photon arrival position relative to telescope [cm]
  cx, cy     direction cosines in the CORSIKA horizontal frame

struct bunch3d:
  x, y, z    photon arrival position relative to telescope [cm]
  cx, cy, cz direction cosines in the CORSIKA horizontal frame
```

Important point:

```text
Photon bunch x/y/z are already relative to the selected telescope center.
They are not absolute array coordinates.
```

For the usual 2D CORSIKA bunch format, LACT_sim reconstructs:

```text
cz = -sqrt(1 - cx^2 - cy^2)
```

The negative sign means the photon is travelling downward in the CORSIKA
vertical coordinate.

### 3. Telescope Pointing In The CORSIKA Frame

For CORSIKA/EventIO mode, LACT_sim uses an azimuth convention compatible with
hessio/sim_telarray:

```text
pointing_az_deg = 0    points to +x, magnetic North
pointing_az_deg = 90   points to -y, East
pointing_az_deg = 180  points to -x, South
pointing_az_deg = 270  points to +y, West
pointing_el_deg        elevation above horizon
zenith_deg             90 - pointing_el_deg
```

The telescope optical-axis unit vector in CORSIKA NWU coordinates is:

```text
z_axis = ( cos(el) cos(az),
          -cos(el) sin(az),
           sin(el) )
```

For a CORSIKA file named or configured as `zenith=20 deg, azimuth=0 deg`, use:

```ini
telescope.pointing_az_deg=0
telescope.pointing_el_deg=70
```

Again, this azimuth zero is CORSIKA magnetic north unless an external
declination correction is explicitly applied.

### 4. Telescope-Local Optical Frame

The optical ray tracer works in a telescope-local frame:

```text
local +z : telescope optical axis, from mirror toward sky
local -z : on-axis incoming photon direction, from sky toward mirror
local +x : local camera/image horizontal axis chosen from the pointing basis
local +y : completes a right-handed frame
```

For:

```ini
source.eventio_coordinate_frame=corsika_iact
```

LACT_sim constructs this CORSIKA-NWU-to-local basis:

```text
local x_axis = ( -sin(el) cos(az),  sin(el) sin(az), cos(el) )
local y_axis = ( -sin(az),         -cos(az),         0       )
local z_axis = (  cos(el) cos(az), -cos(el) sin(az), sin(el) )
```

Then each EventIO photon bunch is rotated by dot products:

```text
local_pos = ( dot(pos_NWU, x_axis),
              dot(pos_NWU, y_axis),
              dot(pos_NWU, z_axis) )

local_dir = ( dot(dir_NWU, x_axis),
              dot(dir_NWU, y_axis),
              dot(dir_NWU, z_axis) )
```

Because EventIO photon positions are already telescope-relative, the standard
`corsika_iact` mode does not subtract the EventIO telescope position. Telescope
positions are still written to HDF5 for provenance and array-level plotting.

### 5. Array Layout Plots

Array-layout plots use East on the horizontal axis and North on the vertical
axis:

```text
plot_x = East  = -array_y_west_m
plot_y = North =  array_x_north_m
```

Core positions are plotted the same way:

```text
core_East  = -core_y_west_m
core_North =  core_x_north_m
```

The plot label deliberately says `HDF5 positions: x=N, y=W; plotted as East vs
North` to prevent confusing the stored NWU values with the display axes.

## Core Position And MC_TELOFF

CORSIKA files with `CSCAT` can reuse one shower with several randomly shifted
array positions. EventIO stores those shifts in `MC_TELOFF`.

hessio defines `MC_TELOFF` as:

```text
offset of the detector array with respect to the shower core
```

Therefore `MC_TELOFF` is not the core position itself. To get the core position
in the same coordinate system as the input `TELESCOPE` positions:

```text
xcore_input = -xoff
ycore_input = -yoff
```

LACT_sim writes this corrected value to:

```text
/events/corsika/core_x_north_m
/events/corsika/core_y_west_m
```

The raw shower-header core is still preserved once per CORSIKA shower in:

```text
/events/corsika_showers
```

That table is mainly for provenance and debugging. For event-level analysis and
array plots, use `/events/corsika`.

Example for the default test file:

```text
event_id = 46889801
MC_TELOFF array_id=1: xoff = 586.0925 m, yoff = 92.5612 m

stored core:
  core_x_north_m = -586.0925 m
  core_y_west_m  =  -92.5612 m

plotted:
  North = -586.0925 m
  East  =   92.5612 m
```

## CORSIKA Direction Angles

The CORSIKA event header contains shower `THETA`, `PHI`, and `ARRANG`. LACT_sim
stores:

```text
theta_deg
phi_deg
array_rotation_deg
azimuth_north_to_east_deg
```

The hessio-compatible conversion used in the reader is:

```text
azimuth_north_to_east = ARRANG - PHI + 180 deg
```

This is an arrival/pointing convention in the CORSIKA frame. `ARRANG` is the
CORSIKA IACT array arrangement angle. It is not the telescope optical pointing
by itself. The optical trace uses:

```ini
telescope.pointing_az_deg
telescope.pointing_el_deg
```

to build the telescope-local optical frame.

## Camera / Whiteboard Plane

The whiteboard and pixel camera lie on the output plane in telescope-local
coordinates. The current LACT optical convention is:

```text
mirror facets are near z = -16 m
focal plane is near z = -8 m
camera normal points toward the mirror: (0, 0, -1)
```

Recommended configuration:

```ini
output.plane_point=0,0,-8
output.plane_normal=0,0,-1
output.plane_u_axis=1,0,0
output.plane_v_axis=0,1,0
```

Camera configs normally include `configs/outputs/focal_plane_f8.cfg`; whiteboard
tests normally include `configs/outputs/whiteboard_f8.cfg`. They define the
same 8 m plane, but the names distinguish whether the plane is being used as a
real camera face or as a virtual diagnostic whiteboard.

The normal defines the camera front-face direction. Ray-plane intersection is
two-sided, and optical incidence uses `abs(dot(ray, normal))`, so changing only
the normal sign does not by itself change whether photons hit the plane.

The image coordinates are set by `u/v`:

```text
image x = hit.u_m = dot(hit_position - plane_point, plane_u_axis)
image y = hit.v_m = dot(hit_position - plane_point, plane_v_axis)
```

The camera pixel CSV uses the same output-plane coordinates, stored as:

```text
x_m, y_m
```

These are camera/image coordinates, not CORSIKA array coordinates.

## Synthetic Parallel-Beam Angles

Analytic optical tests do not use CORSIKA coordinates. For:

```ini
source.mode=ParallelBeam
source.beam_theta_deg=...
source.beam_phi_deg=...
```

the angles are defined directly in the telescope-local optical frame. They
describe photon propagation direction:

```text
beam_direction = ( sin(theta) cos(phi),
                   sin(theta) sin(phi),
                  -cos(theta) )
```

So:

```text
theta = 0 deg  -> on-axis photons, direction (0, 0, -1)
phi = 0 deg    -> positive local x photon component
phi = 90 deg   -> positive local y photon component
phi = 180 deg  -> negative local x photon component
phi = 270 deg  -> negative local y photon component
```

The apparent sky offset is opposite the transverse photon propagation
component. The focal-plane image of a reflecting telescope is inverted relative
to the sky offset; with the recommended `u=+x, v=+y`, the image centroid has
the same sign as the transverse photon direction.

Example:

```ini
source.beam_theta_deg=3
source.beam_phi_deg=270
```

Photons travel with a negative local-y component, corresponding to a sky source
at positive local-y offset. The focal-plane centroid is expected at negative
camera `v`.

## Practical Rules

1. For CORSIKA/EventIO HDF5 metadata, read coordinates as NWU:
   `x=N`, `y=W`, `z=U`.
2. For array plots, display East-North using `East=-West`, `North=North`.
3. For event-level core positions, use `/events/corsika`, not
   `/events/corsika_showers`.
4. For CORSIKA files with `CSCAT`, core position is `-MC_TELOFF`.
5. For telescope pointing, use CORSIKA-compatible azimuth:
   `0=N`, `90=E=-y`.
6. For camera images, use output-plane `u/v`; these are not array coordinates.
7. Magnetic north is not automatically converted to true north. If a future
   analysis needs true geographic coordinates, apply a site magnetic
   declination rotation explicitly and document its sign and epoch.
