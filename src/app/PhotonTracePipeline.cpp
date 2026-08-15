#include "app/PhotonTracePipeline.hpp"

#include <algorithm>
#include <chrono>
#include <stdexcept>

namespace lact {
namespace {

class StageTimer {
public:
    explicit StageTimer(PhotonTraceProfile* profile) : profile_(profile)
    {
        if (profile_) {
            start_ = std::chrono::steady_clock::now();
        }
    }

    void chargeTo(double PhotonTraceProfile::*field)
    {
        if (!profile_) return;
        const auto now = std::chrono::steady_clock::now();
        profile_->*field += std::chrono::duration<double>(now - start_).count();
        start_ = now;
    }

    void restart()
    {
        if (profile_) {
            start_ = std::chrono::steady_clock::now();
        }
    }

private:
    PhotonTraceProfile* profile_ = nullptr;
    std::chrono::steady_clock::time_point start_{};
};

// The datasheet PDE for the tiled S17351 is averaged over the whole package,
// gaps included. When the explicit channel-gap geometry is also applied the
// gap loss would be counted twice, so the PDE is conditioned on being inside
// a channel first.
bool conditionsPdeOnChannelArea(const PhotonTraceContext& context,
                                const Photon& photon)
{
    return !photon.optical_efficiency_preapplied &&
           context.detector != nullptr &&
           context.detector->enabled &&
           context.detector->microcell.enabled &&
           context.detector->microcell.layout == "s17351_tiled_2x4" &&
           context.detector->microcell.pde_includes_inter_channel_gaps;
}

} // namespace

PhotonTraceResult runPhotonTrace(const PhotonTraceContext& context,
                                 const PhotonTraceBunch& bunch,
                                 PhotonCandidate candidate,
                                 PhotonTraceProfile* profile)
{
    PhotonTraceResult result;
    Photon& photon = candidate.photon;
    photon.normalizeDirection();

    // ---- 1. Atmospheric transmission from the emission altitude ----------
    double atmosphere_transmission = 1.0;
    if (context.atmosphere->enabled() &&
        !photon.optical_efficiency_preapplied) {
        atmosphere_transmission = context.atmosphere->transmission(
            photon.wavelength_nm,
            bunch.bunch->emission_altitude_km,
            bunch.global_direction);
    }

    // ---- 2. Wavelength-dependent detection probability -------------------
    double detection_probability =
        context.efficiency->preGeometryDetectionProbability(
            photon.wavelength_nm);
    const bool condition_pde = conditionsPdeOnChannelArea(context, photon);
    const double channel_active_fraction =
        condition_pde
            ? electronics::interChannelActiveFraction(
                  context.detector->microcell)
            : 1.0;
    if (condition_pde && candidate.stochastic) {
        detection_probability /= channel_active_fraction;
        if (detection_probability > 1.0 + 1.0e-12) {
            throw std::runtime_error(
                "channel-gap-conditioned pre-geometry detection probability "
                "exceeds one");
        }
        detection_probability = std::min(1.0, detection_probability);
    }

    result.pre_geometry = context.response_sampler->applyPreGeometry(
        candidate, atmosphere_transmission, detection_probability);
    if (!result.pre_geometry.survives) {
        result.stage = PhotonTraceStage::AbsorbedBeforeOptics;
        return result;
    }

    // ---- 3. Mirror intersection and reflection to the output plane -------
    StageTimer timer(profile);
    const bool backproject =
        bunch.bunch->eventio_2d && context.eventio_2d_backproject;
    result.hit = backproject
        ? context.tracer->traceBackprojectedToPlane(
              photon, *bunch.mirrors, *context.plane, *context.efficiency)
        : context.tracer->traceToPlane(
              photon, *bunch.mirrors, *context.plane, *context.efficiency);
    timer.chargeTo(&PhotonTraceProfile::trace_to_plane_s);

    // ---- 4. Structure blocking the upstream ray --------------------------
    if (result.hit.hit_mirror) {
        // A 2D EventIO position is only an anchor on the incoming line, so the
        // upstream leg is checked on the physical ray rather than on whichever
        // side of the mirror holds the record plane.
        const bool blocked = incomingRayBlockedByObstruction(
            result.hit.mirror_point, photon.dir, *context.obstruction, nullptr);
        timer.chargeTo(&PhotonTraceProfile::obstruction_s);
        if (blocked) {
            result.stage = PhotonTraceStage::BlockedIncoming;
            return result;
        }
    }
    if (!result.hit.hit_surface) {
        result.stage = result.hit.hit_mirror
            ? PhotonTraceStage::ReflectedMissedOutput
            : PhotonTraceStage::MissedMirror;
        return result;
    }

    // ---- 5. Structure blocking the reflected ray -------------------------
    timer.restart();
    const bool blocked_reflected = segmentBlockedByObstruction(
        result.hit.mirror_point, result.hit.surface_point,
        *context.obstruction, nullptr);
    timer.chargeTo(&PhotonTraceProfile::obstruction_s);
    if (blocked_reflected) {
        result.stage = PhotonTraceStage::BlockedReflected;
        return result;
    }
    result.reached_output_plane = true;

    // ---- 6a. No camera: the output plane is the whiteboard ---------------
    if (!context.camera_enabled) {
        if (candidate.stochastic &&
            !context.response_sampler->acceptPostGeometry(candidate,
                                                          result.hit)) {
            result.stage = PhotonTraceStage::RejectedPostGeometry;
            return result;
        }
        result.stage = PhotonTraceStage::ReachedWhiteboard;
        return result;
    }

    // ---- 6b. Camera: pixel lookup, light collector, SiPM geometry, PDE ---
    timer.restart();
    applyCameraResponse(*context.camera, context.light_collector,
                        *context.plane, *context.sipm, *context.electronics,
                        result.hit, context.speed_of_light_m_per_ns,
                        context.detector->enabled ? &context.detector->microcell
                                                  : nullptr,
                        condition_pde && !candidate.stochastic
                            ? 1.0 / channel_active_fraction
                            : 1.0);
    timer.chargeTo(&PhotonTraceProfile::camera_response_s);
    if (!result.hit.hit_camera) {
        result.stage = PhotonTraceStage::LostBetweenPixels;
        return result;
    }

    // ---- 7. Post-geometry acceptance draw --------------------------------
    if (candidate.stochastic &&
        !context.response_sampler->acceptPostGeometry(candidate, result.hit)) {
        result.stage = PhotonTraceStage::RejectedPostGeometry;
        return result;
    }

    result.stage = PhotonTraceStage::Detected;
    return result;
}

} // namespace lact
