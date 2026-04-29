# LACT_sim Coordinate Systems

This note fixes the coordinate conventions used when CORSIKA/EventIO photon
bunches are traced through the LACT optical model.

## CORSIKA / sim_telarray Horizontal Frame

CORSIKA IACT photon bunches store positions in a horizontal detection frame.
The hessio definition is:

```text
struct bunch:
  x, y       arrival position relative to telescope [cm]
  cx, cy     direction cosines

struct bunch3d:
  x, y, z    arrival position relative to telescope [cm]
  cx, cy, cz direction cosines
```

Important details:

```text
x/y/z are already relative to the selected telescope center.
They are not absolute array coordinates.
The axes are still the CORSIKA horizontal detection axes.
The direction cosines are also in that CORSIKA horizontal frame.
```

For the usual 2D bunch format:

```text
cz = -sqrt(1 - cx^2 - cy^2)
```

The negative sign means downward-going photons in the CORSIKA vertical axis.

## Azimuth Convention

For CORSIKA/sim_telarray-compatible pointing, LACT_sim uses azimuth increasing
from North to East:

```text
az = 0 deg    points to +X / North
az = 90 deg   points to -Y / East
az = 180 deg  points to -X / South
az = 270 deg  points to +Y / West
```

The sky direction vector in the CORSIKA horizontal frame is:

```text
z_axis = ( cos(el) cos(az), -cos(el) sin(az), sin(el) )
```

This follows the hessio/sim_telarray convention:

```text
cx = cos(altitude) * cos(azimuth)
cy = cos(altitude) * sin(-azimuth)
cz = sin(altitude)
```

Therefore `zenith = 20 deg, azimuth = 0 deg` corresponds to:

```text
telescope.pointing_az_deg = 0
telescope.pointing_el_deg = 70
```

## LACT Telescope Optical Frame

The ray tracer itself works in a telescope-local optical frame:

```text
local +z : telescope pointing direction, from mirror toward sky
local -z : on-axis incoming photon direction, from sky toward mirror
local +x : camera horizontal axis chosen from the CORSIKA az/el convention
local +y : completes a right-handed frame
```

For `source.eventio_coordinate_frame=corsika_iact`, the code constructs:

```text
local z_axis = ( cos(el) cos(az), -cos(el) sin(az), sin(el) )
local x_axis = ( -sin(el) cos(az), sin(el) sin(az), cos(el) )
local y_axis = ( -sin(az), -cos(az), 0 )
```

Then each EventIO photon bunch is transformed by dot products:

```text
local_pos = (dot(pos, x_axis), dot(pos, y_axis), dot(pos, z_axis))
local_dir = (dot(dir, x_axis), dot(dir, y_axis), dot(dir, z_axis))
```

Because the CORSIKA photon position is already telescope-relative, LACT_sim does
not subtract the telescope position in `corsika_iact` mode.

## Synthetic Parallel-Beam Angles

For analytic optical tests with:

```ini
source.mode=ParallelBeam
source.beam_theta_deg=...
source.beam_phi_deg=...
```

the angles are defined in the LACT telescope-local optical frame. They describe
the photon propagation direction, not the sky-source position:

```text
beam_direction = ( sin(theta) cos(phi),
                   sin(theta) sin(phi),
                  -cos(theta) )
```

Therefore:

```text
theta = 0 deg  -> on-axis photons, direction (0, 0, -1)
phi = 0 deg    -> positive photon x component
phi = 90 deg   -> positive photon y component
phi = 180 deg  -> negative photon x component
phi = 270 deg  -> negative photon y component
```

Since `beam_direction` is the direction in which photons travel, the apparent
sky-source offset is the opposite transverse direction:

```text
sky offset x  ~= -beam_direction.x
sky offset y  ~= -beam_direction.y
```

For a normal reflecting optical system, the focal-plane image is inverted
relative to the sky offset. In the recommended camera coordinates
`u=+x, v=+y`, this gives the same sign as the transverse photon direction:

```text
photons with +x propagation component -> image at +u
photons with +y propagation component -> image at +v
sky source at +x                    -> image at -u
sky source at +y                    -> image at -v
```

Example: a source above the optical axis in the local `+y` sky direction should
be configured as photons traveling with a `-y` component:

```ini
source.beam_theta_deg=3
source.beam_phi_deg=270
```

The expected image centroid is then at negative camera `v`.

## Camera / Whiteboard Plane

The whiteboard and camera live on the same output plane. In the current LACT
local frame:

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

The normal defines which side is the camera front face. Ray-plane intersection
is two-sided, and optical incidence uses `abs(dot(ray, normal))`, so changing
the sign of the normal does not by itself change the number of hits. The image
coordinate convention is set by `plane_u_axis` and `plane_v_axis`:

```text
plot x = hit.u_m = dot(hit_position - plane_point, plane_u_axis)
plot y = hit.v_m = dot(hit_position - plane_point, plane_v_axis)
```

The camera pixel CSV uses the same `(u, v)` coordinates, stored as pixel
`x_m,y_m`.

## Array Coordinates And Telescope Positions

The EventIO telescope table gives telescope positions in the CORSIKA array
coordinate frame. LACT_sim logs these positions for provenance. They are used
only when:

```ini
source.eventio_coordinate_frame=array_global
```

In the standard mode:

```ini
source.eventio_coordinate_frame=corsika_iact
```

the telescope table is not subtracted from photon positions, because the photon
bunch coordinates are already relative to each telescope.

## ARRANG

`ARRANG` is the CORSIKA IACT array arrangement angle. hessio exposes it as:

```text
array_rotation_deg
```

In LACT_sim logs it appears as:

```text
first_event_arrang_deg
```

It describes the rotation angle of the CORSIKA IACT array arrangement / array
offset system in degrees. It is not the telescope pointing direction. In the
EventIO event header, LACT_sim also logs:

```text
first_event_theta_deg
first_event_phi_deg
first_event_az_N_to_E_deg
```

The CORSIKA-to-North-East azimuth conversion used by hessio-compatible code is:

```text
az_N_to_E = ARRANG - PHIP + 180 deg
```

Use `telescope.pointing_az_deg` and `telescope.pointing_el_deg` for the actual
optical-axis rotation into the ray-tracing frame.
