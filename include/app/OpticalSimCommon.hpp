#pragma once

#include <chrono>
#include <cstdint>
#include <functional>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "core/Photon.hpp"
#include "electronics/DetectorPipeline.hpp"
#include "geometry/CameraGeometry.hpp"
#include "geometry/FacetFactory.hpp"
#include "geometry/FacetLayoutUtils.hpp"
#include "geometry/LightCollectorSquareCone.hpp"
#include "geometry/MirrorFacet.hpp"
#ifdef LACT_HAS_HESSIO
#include "io/EventIOPhotonSource.hpp"
#endif
#include "io/PhotonCsvSource.hpp"
#include "io/SurfaceHitCsvWriter.hpp"
#include "io/SyntheticPhotonSource.hpp"
#include "optics/OpticalEfficiency.hpp"
#include "optics/AtmosphereTransmission.hpp"
#include "optics/OpticalSurfaceHit.hpp"
#include "optics/OpticalTracer.hpp"
#include "optics/OutputPlane.hpp"

namespace lact {

extern const double DEG_TO_RAD;

struct WaveformOutputConfig;

struct ComponentConfigPaths {
    std::string telescope;
    std::string mirror;
    std::string source;
    std::string output;
    std::string camera;
    std::string sipm;
    std::string electronics;
    std::string efficiency;
    std::string atmosphere;
    std::string nsb;
    std::string trigger;
    std::string error;
    std::string obstruction;
};

struct ConfigCommandLine {
    std::vector<std::string> positional;
    std::vector<std::pair<std::string, std::string>> overrides;
    bool help = false;
};

struct PropagationConfig {
    double speed_of_light_m_per_ns = 0.299792458;
};

struct TelescopeConfig {
    int id = 0;
    std::string name = "telescope";
    Vec3 position_m{0.0, 0.0, 0.0};
    double pointing_az_deg = 0.0;
    double pointing_el_deg = 90.0;
    double focal_length_m = 8.0;
    std::string coordinate_system = "array";
};

struct TelescopeFrame {
    Vec3 origin{0.0, 0.0, 0.0};
    Vec3 x_axis{1.0, 0.0, 0.0};
    Vec3 y_axis{0.0, 1.0, 0.0};
    Vec3 z_axis{0.0, 0.0, 1.0};

    Vec3 rotateVector(const Vec3& local) const;
    Vec3 pointToGlobal(const Vec3& local) const;
    Vec3 rotateVectorToLocal(const Vec3& global) const;
    Vec3 pointToLocal(const Vec3& global) const;
};

struct ErrorConfig {
    std::uint64_t random_seed = 987654321ULL;
    double facet_radial_position_sigma_m = 0.0;
    double facet_normal_sigma_deg = 0.0;
    double reflect_direction_sigma_deg = 0.0;
    double radius_of_curvature_sigma_m = 0.0;
    double reflectivity_scale_sigma = 0.0;
    std::string structural_deformation_config;
};

struct ObstructionMask {
    bool enabled = false;
    std::string mode = "mask";
    std::string mask_csv;
    std::string primitives_csv;
    bool check_incoming = true;
    bool check_reflected = false;
    bool mark_only = false;
    double plane_z_m = -16.0;
    double x_min_m = -4.5;
    double y_min_m = -4.5;
    double cell_size_m = 0.02;
    int nx = 0;
    int ny = 0;
    std::vector<std::uint8_t> blocked;

    struct Primitive {
        std::string type;
        std::string name;
        std::string role = "default";
        std::string material_id = "default";
        Vec3 p0;
        Vec3 p1;
        double radius_m = 0.0;
        Vec3 half_size;
        Vec3 bbox_min;
        Vec3 bbox_max;
        Vec3 center;
        double height_m = 0.0;
        double rotation_rad = 0.0;
        int sides = 0;
        double hole_radius_m = 0.0;
        double hole_rotation_rad = 0.0;
        int hole_sides = 0;
        bool has_hole = false;
    };

    std::vector<Primitive> primitives;

    bool contains(double x_m, double y_m) const;
};

struct CameraConfig {
    bool enabled = false;
    std::string mode = "none";
    std::string csv_path;
    std::string pixel_shape = "hexagonal";
    double pixel_size_m = 0.05;
    double pixel_pitch_m = 0.05;
    double radius_m = 0.5;
    std::string collector = "none";
    std::string collector_material = "true_reflect";
    std::string collector_reflectivity_csv;
    double collector_entrance_size_m = 0.0;
    double collector_exit_size_m = 0.0130;
    double collector_height_m = 0.0297;
};

struct SipmConfig {
    double size_m = 0.0130;
};

struct ElectronicsConfig {
    // Integrated SiPM microcell saturation.  This operates directly on the
    // detected p.e. image and therefore does not require waveform output.
    int channels_per_pixel = 8;
    int microcells_per_channel = 33792;
    bool saturation_enabled = true;
    std::string saturation_model = "hard_no_recovery";
};

struct NsbConfig {
    bool enabled = false;
    std::string model = "constant_rate";
    double rate_pe_per_ns_per_pixel = 0.0;
    // NSB integration gate used by the saved camera image. The gate is
    // [0, window_ns) relative to the event/waveform reference time. Waveform
    // generation and triggering may use a wider interval.
    double window_ns = 32.0;
    std::uint64_t seed = 12345ULL;
    std::string spectrum_csv;
    std::string spectrum_unit = "ph_s_nm_sr_m2";
    double effective_area_m2 = 0.0;
    // Angular/spectral average probability for a photon at the collector
    // entrance to reach the SiPM exit. Applied only when spectral_flux is
    // converted into a detected-p.e. rate; constant_rate is already in p.e.
    double collector_mean_transmission = 1.0;
    bool collector_mean_transmission_configured = false;
    double microcell_geometric_acceptance = 1.0;
    std::string pixel_solid_angle = "auto";
    double pixel_solid_angle_sr = 0.0;
    bool computed_from_spectrum = false;
    double spectral_integral_pe_s_sr_m2 = 0.0;
};

struct TriggerConfig {
    bool enabled = false;
    // Camera triggering remains the default retained-event decision.  Array
    // coincidence is an independent optional stage and is disabled by
    // default for current LACT production.
    bool array_enabled = false;
    double pixel_threshold_pe = 5.0;
    int camera_multiplicity = 3;
    int array_multiplicity = 2;
    double coincidence_window_ns = 20.0;
    double camera_coincidence_window_ns = 20.0;
    double array_coincidence_window_ns = 0.0;
    // none: compare raw local camera-trigger times.
    // plane_wave: add the CORSIKA-NWU planar wavefront delay before the
    // array coincidence test.
    std::string array_time_correction = "none";
    // Zero selects the sim_telarray-compatible air speed derived from the
    // EventIO observation altitude. A positive value is an explicit override.
    double array_wavefront_speed_m_per_ns = 0.0;
};

class ElectronicsResponse {
public:
    ElectronicsResponse() = default;
    explicit ElectronicsResponse(const ElectronicsConfig& cfg);
    double peConversion(double wavelength_nm) const;
    double saturatedPe(double primary_pe) const;
    std::vector<double> saturatedWaveform(
        const std::vector<double>& primary_pe_by_time_bin) const;
    double totalMicrocellsPerPixel() const;

private:
    ElectronicsConfig cfg_;
};

struct SourceRuntimeConfig {
    bool use_photon_csv = false;
    bool use_eventio = false;
    // Explicit input-frame names; these do not redefine trace/output axes:
    //   telescope_local      existing local optical input coordinates
    //   corsika_nwu_relative CORSIKA NWU positions relative to the telescope
    //   corsika_nwu_global   absolute CORSIKA NWU array positions
    //   enu_east_relative    ENU positions relative to the telescope; az=0 East
    //   enu_east_global      absolute ENU array positions; az=0 East
    //   lact_generic_global  LACT global NWU; azimuth North-to-East
    std::string coordinate_frame = "telescope_local";
    bool csv_local_telescope_frame = true;
    std::string csv_path;
    std::string eventio_path;
    std::string eventio_coordinate_frame = "corsika_iact";
    // User-facing scalar telescope-local z offset from the EventIO
    // telescope/fiducial reference to the imported LACT optical origin.
    // Internally this is the translation (0,0,z) and applies to both 2D and
    // 3D EventIO positions after their input frame is converted to local.
    double eventio_reference_z_m = -16.0;
    std::string eventio_2d_plane_mode = "auto";
    std::string event_id_mode = "event";
    bool use_eventio_telescope_position = true;
    bool filter_telescope_id = false;
    int selected_telescope_id = 0;
    bool filter_event_id = false;
    int selected_event_id = 0;
    std::vector<int> selected_event_ids;
    bool filter_shower_event_id = false;
    int selected_shower_event_id = 0;
    int max_shower_events = 0;
};

struct CollectorTraceResult {
    bool hit_sipm = false;
    double intensity = 1.0;
    int reflections = 0;
    double path_length_m = 0.0;
    bool reflection_limit_reached = false;
    Cone::Position exit_position{0.0, 0.0, 0.0};
    Cone::DirectionVecter exit_direction{0.0, 0.0, 0.0};
};

struct PixelAccumulator {
    int event_id = 0;
    int telescope_id = 0;
    int pixel_id = -1;
    std::uint64_t photon_count = 0;
    double pe = 0.0;
    double signal = 0.0;
    double time_sum = 0.0;
    double time2_sum = 0.0;
};

using PixelKey = std::tuple<int, int, int>;

std::string trim(const std::string& s);
std::string lowerCopy(std::string s);
bool startsWith(const std::string& text, const std::string& prefix);
std::string resolveRelativePath(const std::string& base_config_path, const std::string& path);
std::map<std::string, std::string> readKeyValueConfig(const std::string& path);
ConfigCommandLine parseConfigCommandLine(int argc, char** argv);
void applyConfigOverrides(
    std::map<std::string, std::string>& cfg,
    const std::vector<std::pair<std::string, std::string>>& overrides);
std::string scopedComponentKey(const std::string& key, const std::string& prefix);
std::map<std::string, std::string> expandConfig(const std::map<std::string, std::string>& main_cfg,
                                                const std::string& main_config_path,
                                                ComponentConfigPaths& paths);
std::string getString(const std::map<std::string, std::string>& cfg,
                      const std::string& key,
                      const std::string& fallback);
double getDouble(const std::map<std::string, std::string>& cfg,
                 const std::string& key,
                 double fallback);
int getInt(const std::map<std::string, std::string>& cfg,
           const std::string& key,
           int fallback);
std::uint64_t getUInt64(const std::map<std::string, std::string>& cfg,
                        const std::string& key,
                        std::uint64_t fallback);
bool getBool(const std::map<std::string, std::string>& cfg,
             const std::string& key,
             bool fallback);
Vec3 parseVec3(const std::string& text, const std::string& key);
Vec3 getVec3(const std::map<std::string, std::string>& cfg,
             const std::string& key,
             const Vec3& fallback);
DishType parseDishType(const std::string& text);
SurfaceType parseSurfaceType(const std::string& text);
ApertureShape parseApertureShape(const std::string& text);
PixelShape parsePixelShape(const std::string& text);
std::string pixelShapeName(PixelShape shape);
SyntheticMode parseSyntheticMode(const std::string& text);
bool isPhotonCsvMode(const std::string& text);
bool isEventIOMode(const std::string& text);
std::string vec3ToString(const Vec3& v);
std::string sourceModeName(SyntheticMode mode);
TelescopeFrame buildTelescopeFrame(const TelescopeConfig& telescope);
TelescopeFrame buildCorsikaNwuTelescopeFrame(const TelescopeConfig& telescope);
TelescopeFrame buildEnuEastTelescopeFrame(const TelescopeConfig& telescope);
std::string normalizeSourceCoordinateFrame(const std::string& frame_name);
std::string sourceCoordinateFrameDescription(const std::string& frame_name);
PhotonBunch transformBunchToTelescopeLocal(const PhotonBunch& input,
                                           const TelescopeConfig& telescope,
                                           const std::string& frame_name);
void applyEventIOReferenceZOffset(PhotonBunch& bunch,
                                  double eventio_reference_z_m);
Vec3 sourceDirectionInWorld(const PhotonBunch& input,
                            const TelescopeConfig& telescope,
                            const std::string& frame_name);
void applyTelescopeFrame(std::vector<MirrorFacet>& facets,
                         OutputPlane& plane,
                         const TelescopeFrame& frame);
void applyTelescopeFrame(Photon& photon, const TelescopeFrame& frame);
double elapsedSeconds(std::chrono::steady_clock::time_point start,
                      std::chrono::steady_clock::time_point stop);
void printSection(const std::string& title);
void printField(const std::string& label, const std::string& value);
std::string doubleToString(double value, int precision = 6);
std::string intToString(std::uint64_t value);
bool parseDoubleStrict(const std::string& text, double& out);
bool isDisabledText(const std::string& text);
std::string factorDescription(const EfficiencyFactorConfig& factor);
std::vector<double> parseDoubleList(const std::string& text, const std::string& key);
std::vector<int> parseIntList(const std::string& text, const std::string& key);
std::vector<std::string> splitCsvCells(const std::string& line);
CameraGeometry readCameraCsv(const std::string& path);
std::vector<std::pair<double, double>> readCollectorReflectivityCsv(const std::string& path);
CameraGeometry makeGeneratedCamera(const CameraConfig& cfg);
std::vector<MirrorFacet> buildElevationSeriesFacets(const std::map<std::string, std::string>& cfg);
std::vector<MirrorFacet> buildFacetsFromConfig(const std::map<std::string, std::string>& cfg);
SyntheticPhotonConfig buildSourceConfig(const std::map<std::string, std::string>& cfg);
SourceRuntimeConfig buildSourceRuntimeConfig(const std::map<std::string, std::string>& cfg);
#ifdef LACT_HAS_HESSIO
EventIOPhotonConfig buildEventIOPhotonConfig(const std::map<std::string, std::string>& cfg,
                                             const SyntheticPhotonConfig& source_cfg,
                                             const SourceRuntimeConfig& runtime_cfg);
void printEventIOMetadata(const EventIOMetadata& metadata,
                          const SourceRuntimeConfig& runtime_cfg);
#endif
PhotonCsvConfig buildPhotonCsvConfig(const std::map<std::string, std::string>& cfg,
                                     const SyntheticPhotonConfig& source_cfg,
                                     const SourceRuntimeConfig& runtime_cfg);
OutputPlane buildOutputPlane(const std::map<std::string, std::string>& cfg);
CameraConfig buildCameraConfig(const std::map<std::string, std::string>& cfg);
SipmConfig buildSipmConfig(const std::map<std::string, std::string>& cfg);
ElectronicsConfig buildElectronicsConfig(const std::map<std::string, std::string>& cfg);
electronics::DetectorPipelineConfig buildDetectorPipelineConfig(
    const std::map<std::string, std::string>& cfg);
void validateCameraDetectorCompatibility(
    const CameraConfig& camera,
    const electronics::DetectorPipelineConfig& detector);
NsbConfig buildNsbConfig(const std::map<std::string, std::string>& cfg);
void resolveNsbSpectralRate(NsbConfig& nsb,
                            const OpticalEfficiencyConfig& efficiency_cfg,
                            const CameraGeometry& camera,
                            const TelescopeConfig& telescope,
                            const electronics::DetectorPipelineConfig* detector = nullptr);
void generateIntegratedNsbPe(const NsbConfig& nsb,
                             int event_id,
                             int telescope_id,
                             std::size_t n_pixels,
                             double window_ns,
                             const std::function<void(std::size_t, float)>& add_sample);
bool nsbTimeInImageWindow(const NsbConfig& nsb, double relative_time_ns);
void generateTimeBinnedNsbPe(const NsbConfig& nsb,
                             const WaveformOutputConfig& waveform_cfg,
                             int event_id,
                             int telescope_id,
                             std::size_t n_pixels,
                             std::size_t n_bins,
                             const std::function<void(std::size_t, std::size_t, float)>& add_sample);
TriggerConfig buildTriggerConfig(const std::map<std::string, std::string>& cfg);
CameraGeometry buildCameraGeometry(const CameraConfig& cfg);
double cameraPixelSizeForCollector(const CameraConfig& cfg, const CameraGeometry& camera);
std::unique_ptr<Cone::SquareCone> buildLightCollector(const CameraConfig& cfg,
                                                      const CameraGeometry& camera);
const CameraPixel* findContainingPixel(const CameraGeometry& camera, double x, double y);
CollectorTraceResult traceLightCollector(const Cone::SquareCone& cone,
                                         const OutputPlane& plane,
                                         const OpticalSurfaceHit& hit,
                                         const CameraPixel& pixel,
                                         const SipmConfig& sipm);
void applyCameraResponse(const CameraGeometry& camera,
                         const Cone::SquareCone* light_collector,
                         const OutputPlane& plane,
                         const SipmConfig& sipm,
                         const ElectronicsResponse& electronics,
                         OpticalSurfaceHit& hit,
                         double speed_of_light_m_per_ns = 0.299792458,
                         const ::lact::electronics::MicrocellConfig* microcell = nullptr,
                         double post_geometry_pde_scale = 1.0);
void accumulatePixelHit(std::map<PixelKey, PixelAccumulator>& pixels,
                        int event_id,
                        int telescope_id,
                        const OpticalSurfaceHit& hit);
void writePixelCsv(const std::string& path,
                   const std::map<PixelKey, PixelAccumulator>& pixels);
TelescopeConfig buildTelescopeConfig(const std::map<std::string, std::string>& cfg);
ErrorConfig buildErrorConfig(const std::map<std::string, std::string>& cfg);
ObstructionMask buildObstructionMask(const std::map<std::string, std::string>& cfg);
bool photonBlockedByObstruction(const Photon& photon,
                                const ObstructionMask& obstruction,
                                const TelescopeFrame* trace_to_local_frame = nullptr);
bool incomingSegmentBlockedByObstruction(const Vec3& a,
                                         const Vec3& b,
                                         const ObstructionMask& obstruction,
                                         const TelescopeFrame* trace_to_local_frame = nullptr);
bool incomingRayBlockedByObstruction(const Vec3& mirror_point,
                                     const Vec3& incoming_direction,
                                     const ObstructionMask& obstruction,
                                     const TelescopeFrame* trace_to_local_frame = nullptr);
bool segmentBlockedByObstruction(const Vec3& a,
                                 const Vec3& b,
                                 const ObstructionMask& obstruction,
                                 const TelescopeFrame* trace_to_local_frame = nullptr);
void applyStructuralDeformation(std::vector<MirrorFacet>& facets,
                                const ErrorConfig& error,
                                const TelescopeConfig& telescope);
void applyFacetErrors(std::vector<MirrorFacet>& facets, const ErrorConfig& error);
void applyFacetEfficiencyScales(std::vector<MirrorFacet>& facets,
                                const std::map<std::string, std::string>& cfg);
double effectiveReflectDirectionSigmaRad(const std::vector<MirrorFacet>& facets,
                                         const ErrorConfig& error);
OpticalEfficiencyConfig buildEfficiencyConfig(const std::map<std::string, std::string>& cfg);
PropagationConfig buildPropagationConfig(const std::map<std::string, std::string>& cfg);

} // namespace lact
