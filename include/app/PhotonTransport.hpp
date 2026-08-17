#pragma once

// One traced photon, start to finish. This transport layer owns no run-level
// accumulators or output streams.

#include "app/OpticalSimCommon.hpp"
#include "app/PhotonResponseSampler.hpp"
#include "core/PhotonBunch.hpp"
#include "electronics/DetectorPipeline.hpp"
#include "geometry/CameraGeometry.hpp"
#include "geometry/LightCollectorSquareCone.hpp"
#include "geometry/MirrorLayout.hpp"
#include "optics/AtmosphereTransmission.hpp"
#include "optics/OpticalEfficiency.hpp"
#include "optics/OpticalSurfaceHit.hpp"
#include "optics/OpticalTracer.hpp"
#include "optics/OutputPlane.hpp"

namespace lact {

// Where a photon stopped. Each value is one exit point of tracePhoton(),
// in the order the stages are attempted.
enum class PhotonTraceStage {
    // Absorbed by the atmosphere, or lost to the wavelength-dependent
    // detection probability, before any geometry was considered.
    AbsorbedBeforeOptics,
    // The incoming ray met no mirror facet.
    MissedMirror,
    // Structure blocked the upstream ray between sky and mirror.
    BlockedIncoming,
    // Reflected off a facet but never crossed the output plane.
    ReflectedMissedOutput,
    // Structure blocked the reflected ray between mirror and output plane.
    BlockedReflected,
    // Reached the output plane with no camera configured (whiteboard mode).
    ReachedWhiteboard,
    // Reached the output plane but landed outside every pixel/collector.
    LostBetweenPixels,
    // Survived the geometry but lost the post-geometry acceptance draw.
    RejectedPostGeometry,
    // Detected: counted as signal in a pixel.
    Detected,
};

// Read-only context that never changes while photons are traced.
struct PhotonTraceContext {
    const OpticalTracer* tracer = nullptr;
    const OutputPlane* plane = nullptr;
    const OpticalEfficiency* efficiency = nullptr;
    const AtmosphereTransmission* atmosphere = nullptr;
    const ObstructionMask* obstruction = nullptr;
    const CameraGeometry* camera = nullptr;
    const Cone::SquareCone* light_collector = nullptr;
    const SipmConfig* sipm = nullptr;
    const ElectronicsResponse* electronics = nullptr;
    const electronics::DetectorPipelineConfig* detector = nullptr;
    const PhotonResponseSampler* response_sampler = nullptr;
    double speed_of_light_m_per_ns = 0.299792458;
    bool camera_enabled = false;
    bool normalize_direction = true;
    const TelescopeFrame* obstruction_frame = nullptr;
    bool obstruction_mark_only = false;
    bool incoming_obstruction_uses_ray = false;
    // EventIO 2D records are anchors on a line, not creation points, so the
    // mirror may sit at either sign of the ray parameter.
    bool eventio_2d_backproject = false;
};

// The bunch the candidate came from, resolved once per bunch.
struct PhotonTraceBunch {
    const PhotonBunch* bunch = nullptr;
    const MirrorLayout* mirrors = nullptr;
    // Propagation direction in the world frame; only its up-component is used,
    // by the horizontally symmetric atmosphere model.
    Vec3 global_direction;
};

// Optional per-stage timing, accumulated in seconds. Null disables it.
struct PhotonTraceProfile {
    double trace_to_plane_s = 0.0;
    double obstruction_s = 0.0;
    double camera_response_s = 0.0;
};

struct PhotonTraceResult {
    PhotonTraceStage stage = PhotonTraceStage::AbsorbedBeforeOptics;
    OpticalSurfaceHit hit;
    PreGeometryResponse pre_geometry;
    // True once both obstruction checks passed and the photon is on the
    // output plane, whatever happens to it afterwards.
    bool reached_output_plane = false;
};

// Trace one candidate. `candidate` is taken by value because the acceptance
// stages consume its remaining probability as they go.
PhotonTraceResult tracePhoton(const PhotonTraceContext& context,
                              const PhotonTraceBunch& bunch,
                              PhotonCandidate candidate,
                              PhotonTraceProfile* profile);

} // namespace lact
