#include "app/OpticalSimCommon.hpp"

#ifdef LACT_HAS_HDF5
#include <hdf5.h>
#endif

#include <algorithm>
#include <cstring>
#include <filesystem>
#include <random>
#include <tuple>

using namespace lact;

namespace {

struct CorsikaTraceOutputConfig {
    std::string hits_csv = "corsika_whiteboard_hits.csv";
    std::string pixel_csv = "corsika_pixel_image.csv";
    std::string summary_csv = "corsika_trace_summary.csv";
    std::string hdf5_path = "corsika_trace.h5";
    std::string format = "hdf5";
    std::string hdf5_storage = "dense";
    bool hdf5_write_components = false;
    bool save_only_triggered = false;
};

struct TraceSummary {
    int event_id = 0;
    int telescope_id = 0;
    std::uint64_t input_bunches = 0;
    double input_photons = 0.0;
    std::uint64_t hit_mirror = 0;
    std::uint64_t hit_output_plane = 0;
    std::uint64_t hit_camera = 0;
    std::uint64_t accepted_camera = 0;
    std::uint64_t lost_between_pixels = 0;
    std::set<int> unique_pixels;
    double weighted_signal = 0.0;
    double weighted_time_sum = 0.0;
    double weighted_time2_sum = 0.0;
};

struct TelescopeEventAccumulator {
    int telescope_id = 0;
    std::set<int> output_events;
    std::uint64_t input_bunches = 0;
    double input_photons = 0.0;
    std::uint64_t hit_mirror = 0;
    std::uint64_t hit_output_plane = 0;
    std::uint64_t hit_camera = 0;
    std::uint64_t accepted_camera = 0;
    std::uint64_t lost_between_pixels = 0;
    std::set<int> unique_pixels;
    double weighted_signal = 0.0;
    double weighted_time_sum = 0.0;
    double weighted_time2_sum = 0.0;
};

using SummaryKey = std::pair<int, int>;

struct WhiteboardHdf5Row {
    std::int64_t event_id;
    std::int32_t telescope_id;
    std::int64_t photon_index;
    std::int32_t mirror_id;
    float surface_x_m;
    float surface_y_m;
    float surface_z_m;
    float u_m;
    float v_m;
    float dir_x;
    float dir_y;
    float dir_z;
    float time_ns;
    float wavelength_nm;
    float weight;
    float relative_efficiency;
    float signal_weight;
};

TelescopeFrame frameForEventIOTelescope(const TelescopeConfig& base_telescope,
                                        const EventIOMetadata& metadata,
                                        int telescope_id,
                                        bool use_eventio_position)
{
    TelescopeConfig telescope = base_telescope;
    if (use_eventio_position) {
        if (auto tel = metadata.telescopeById(telescope_id)) {
            telescope.position_m = {tel->x_m, tel->y_m, tel->z_m};
        }
    }
    return buildTelescopeFrame(telescope);
}

TelescopeFrame corsikaIactFrame(const TelescopeConfig& telescope)
{
    // CORSIKA/sim_telarray convention used by hessio:
    //   global z is upward, azimuth is N -> E, and the horizontal y component
    //   enters sky direction vectors with a minus sign. Photon bunch x/y are
    //   already relative to a telescope position, but still expressed in this
    //   horizontal detection frame. We rotate them into the LACT optical frame
    //   where local +z points to the sky and photons arrive roughly along -z.
    const double az = telescope.pointing_az_deg * DEG_TO_RAD;
    const double el = telescope.pointing_el_deg * DEG_TO_RAD;
    const double sin_el = std::sin(el);
    const double cos_el = std::cos(el);
    const double sin_az = std::sin(az);
    const double cos_az = std::cos(az);

    TelescopeFrame frame;
    frame.origin = {0.0, 0.0, 0.0};
    frame.x_axis = Vec3{-sin_el * cos_az, sin_el * sin_az, cos_el}.normalized();
    frame.y_axis = Vec3{-sin_az, -cos_az, 0.0}.normalized();
    frame.z_axis = Vec3{cos_el * cos_az, -cos_el * sin_az, sin_el}.normalized();
    return frame;
}

PhotonBunch transformEventIOBunchToTraceFrame(
    const PhotonBunch& input,
    const TelescopeConfig& telescope_cfg,
    const EventIOMetadata& metadata,
    const SourceRuntimeConfig& source_runtime_cfg)
{
    const std::string frame_name = lowerCopy(trim(source_runtime_cfg.eventio_coordinate_frame));
    if (frame_name == "telescope_local" || frame_name == "local") {
        return input;
    }

    PhotonBunch out = input;
    if (frame_name == "corsika_iact" || frame_name == "corsika" ||
        frame_name == "simtelarray") {
        const TelescopeFrame frame = corsikaIactFrame(telescope_cfg);
        out.photon.pos = frame.rotateVectorToLocal(input.photon.pos);
        out.photon.dir = frame.rotateVectorToLocal(input.photon.dir).normalized();
        return out;
    }

    if (frame_name == "array_global" || frame_name == "global") {
        const TelescopeFrame frame = frameForEventIOTelescope(
            telescope_cfg,
            metadata,
            input.telescope_id,
            source_runtime_cfg.use_eventio_telescope_position);
        out.photon.pos = frame.pointToLocal(input.photon.pos);
        out.photon.dir = frame.rotateVectorToLocal(input.photon.dir).normalized();
        return out;
    }

    throw std::runtime_error("unsupported source.eventio_coordinate_frame: " +
                             source_runtime_cfg.eventio_coordinate_frame);
    return out;
}

int showerEventFromOutputEvent(int event_id, const std::string& event_id_mode)
{
    const std::string mode = lowerCopy(trim(event_id_mode));
    if (mode == "event_array100" || mode == "runid") {
        return event_id / 100;
    }
    return event_id;
}

int arrayIdFromOutputEvent(int event_id, const std::string& event_id_mode)
{
    const std::string mode = lowerCopy(trim(event_id_mode));
    if (mode == "event_array100" || mode == "runid") {
        return event_id % 100;
    }
    return 0;
}

struct OutputEventMetadata {
    int event_id = 0;
    int shower_event = 0;
    int array_id = 0;
    double energy_gev = 0.0;
    double core_x_north_m = 0.0;
    double core_y_west_m = 0.0;
    double azimuth_north_to_east_deg = 0.0;
    bool found = false;
    bool used_array_offset = false;
};

OutputEventMetadata outputEventMetadata(int event_id,
                                        const std::string& event_id_mode,
                                        const EventIOMetadata& metadata)
{
    OutputEventMetadata out;
    out.event_id = event_id;
    out.shower_event = showerEventFromOutputEvent(event_id, event_id_mode);
    out.array_id = arrayIdFromOutputEvent(event_id, event_id_mode);

    auto event_it = std::find_if(
        metadata.events.begin(), metadata.events.end(),
        [&out](const EventIOEventHeader& event) {
            return event.shower_event_id == out.shower_event;
        });
    if (event_it == metadata.events.end()) {
        return out;
    }

    out.found = true;
    out.energy_gev = event_it->energy_gev;
    out.core_x_north_m = event_it->core_x_m;
    out.core_y_west_m = event_it->core_y_m;
    out.azimuth_north_to_east_deg = event_it->azimuth_north_to_east_deg;

    if (auto offsets = metadata.arrayOffsetsForShower(out.shower_event)) {
        const std::size_t offset_index = static_cast<std::size_t>(out.array_id);
        if (out.array_id >= 0 && offset_index < offsets->x_m.size() &&
            offset_index < offsets->y_m.size()) {
            // MC_TELOFF stores the offset of the detector array with respect
            // to the shower core.  The core position in the telescope/input
            // array frame is therefore the opposite vector.
            out.core_x_north_m = -offsets->x_m[offset_index];
            out.core_y_west_m = -offsets->y_m[offset_index];
            out.used_array_offset = true;
        }
    }
    return out;
}

bool shouldHideInputCardLine(const std::string& line)
{
    const std::string text = lowerCopy(trim(line));
    return startsWith(text, "telfil") || startsWith(text, "direct");
}

void printEventIOMetadataSummary(const EventIOMetadata& metadata)
{
    printSection("EventIO metadata");
    printField("shower_events", intToString(metadata.events.size()));
    printField("telescopes", intToString(metadata.telescopes.size()));
    if (!metadata.events.empty()) {
        int min_event = metadata.events.front().shower_event_id;
        int max_event = metadata.events.front().shower_event_id;
        for (const auto& event : metadata.events) {
            min_event = std::min(min_event, event.shower_event_id);
            max_event = std::max(max_event, event.shower_event_id);
        }
        printField("shower_event_range",
                   intToString(min_event) + ".." + intToString(max_event));

        const auto& first = metadata.events.front();
        printField("first_event_theta_deg", doubleToString(first.theta_deg));
        printField("first_event_phi_deg", doubleToString(first.phi_deg));
        printField("first_event_arrang_deg", doubleToString(first.array_rotation_deg));
        printField("first_event_az_N_to_E_deg",
                   doubleToString(first.azimuth_north_to_east_deg));
    }

    std::size_t hidden_input_lines = 0;
    printField("input_card_lines", intToString(metadata.input_lines.size()));
    for (std::size_t i = 0; i < metadata.input_lines.size(); ++i) {
        if (shouldHideInputCardLine(metadata.input_lines[i])) {
            ++hidden_input_lines;
            continue;
        }
        printField("input[" + intToString(i) + "]", metadata.input_lines[i]);
    }
    if (hidden_input_lines > 0) {
        printField("input_hidden_path_lines", intToString(hidden_input_lines));
    }
}

CorsikaTraceOutputConfig buildCorsikaTraceOutputConfig(
    const std::map<std::string, std::string>& cfg)
{
    CorsikaTraceOutputConfig out;
    out.hits_csv = getString(cfg, "output.hits_csv",
                             getString(cfg, "output.whiteboard_csv", out.hits_csv));
    out.pixel_csv = getString(cfg, "output.pixel_csv", out.pixel_csv);
    out.summary_csv = getString(cfg, "output.summary_csv", out.summary_csv);
    out.hdf5_path = getString(cfg, "output.hdf5_path",
                              getString(cfg, "output.h5_path", out.hdf5_path));
    out.format = lowerCopy(trim(getString(cfg, "output.format", out.format)));
    out.hdf5_storage =
        lowerCopy(trim(getString(cfg, "output.hdf5_storage", out.hdf5_storage)));
    out.hdf5_write_components =
        getBool(cfg, "output.hdf5_write_components", out.hdf5_write_components);
    out.save_only_triggered =
        getBool(cfg, "output.save_only_triggered", out.save_only_triggered);
    if (out.format.empty()) {
        out.format = "hdf5";
    }
    if (!(out.format == "hdf5" || out.format == "h5" ||
          out.format == "csv" || out.format == "both")) {
        throw std::runtime_error("output.format must be hdf5, csv, or both");
    }
    if (out.hdf5_storage.empty()) {
        out.hdf5_storage = "dense";
    }
    if (!(out.hdf5_storage == "sparse" ||
          out.hdf5_storage == "dense" ||
          out.hdf5_storage == "both")) {
        throw std::runtime_error("output.hdf5_storage must be sparse, dense, or both");
    }
    return out;
}

bool outputWantsCsv(const CorsikaTraceOutputConfig& cfg)
{
    return cfg.format == "csv" || cfg.format == "both";
}

bool outputWantsHdf5(const CorsikaTraceOutputConfig& cfg)
{
    return cfg.format == "hdf5" || cfg.format == "h5" || cfg.format == "both";
}

void writeCorsikaWhiteboardHeader(std::ofstream& ofs)
{
    ofs << std::setprecision(10);
    ofs << "event_id,telescope_id,photon_index,mirror_id,"
        << "surface_x_m,surface_y_m,surface_z_m,"
        << "u_m,v_m,dir_x,dir_y,dir_z,"
        << "time_ns,wavelength_nm,weight,relative_efficiency,"
        << "signal_weight\n";
}

void writeCorsikaWhiteboardHit(std::ofstream& ofs,
                               const PhotonBunch& bunch,
                               std::uint64_t photon_index,
                               const OpticalSurfaceHit& hit)
{
    const double signal = hit.weight * hit.relative_efficiency;
    ofs << bunch.event_id << ","
        << bunch.telescope_id << ","
        << photon_index << ","
        << hit.mirror_id << ","
        << hit.surface_point.x << ","
        << hit.surface_point.y << ","
        << hit.surface_point.z << ","
        << hit.u_m << ","
        << hit.v_m << ","
        << hit.out_dir.x << ","
        << hit.out_dir.y << ","
        << hit.out_dir.z << ","
        << hit.time_ns << ","
        << hit.wavelength_nm << ","
        << hit.weight << ","
        << hit.relative_efficiency << ","
        << signal << "\n";
}

WhiteboardHdf5Row makeWhiteboardHdf5Row(const PhotonBunch& bunch,
                                        std::uint64_t photon_index,
                                        const OpticalSurfaceHit& hit)
{
    const double signal = hit.weight * hit.relative_efficiency;
    return WhiteboardHdf5Row{
        static_cast<std::int64_t>(bunch.event_id),
        static_cast<std::int32_t>(bunch.telescope_id),
        static_cast<std::int64_t>(photon_index),
        static_cast<std::int32_t>(hit.mirror_id),
        static_cast<float>(hit.surface_point.x),
        static_cast<float>(hit.surface_point.y),
        static_cast<float>(hit.surface_point.z),
        static_cast<float>(hit.u_m),
        static_cast<float>(hit.v_m),
        static_cast<float>(hit.out_dir.x),
        static_cast<float>(hit.out_dir.y),
        static_cast<float>(hit.out_dir.z),
        static_cast<float>(hit.time_ns),
        static_cast<float>(hit.wavelength_nm),
        static_cast<float>(hit.weight),
        static_cast<float>(hit.relative_efficiency),
        static_cast<float>(signal),
    };
}

void writeSummaryCsv(const std::string& path,
                     const std::map<SummaryKey, TraceSummary>& summaries)
{
    const std::filesystem::path out_path(path);
    if (out_path.has_parent_path()) {
        std::filesystem::create_directories(out_path.parent_path());
    }
    std::ofstream ofs(path);
    if (!ofs) {
        throw std::runtime_error("failed to write summary CSV: " + path);
    }
    ofs << std::setprecision(10);
    ofs << "event_id,telescope_id,input_bunches,input_photons,"
        << "hit_mirror,hit_output_plane,hit_camera,accepted_camera,"
        << "lost_between_pixels,unique_hit_pixels,pe,signal,"
        << "time_mean_ns,time_rms_ns\n";
    for (const auto& kv : summaries) {
        const auto& s = kv.second;
        const double mean = s.weighted_signal > 0.0
            ? s.weighted_time_sum / s.weighted_signal
            : 0.0;
        const double var = s.weighted_signal > 0.0
            ? std::max(0.0, s.weighted_time2_sum / s.weighted_signal - mean * mean)
            : 0.0;
        ofs << s.event_id << ","
            << s.telescope_id << ","
            << s.input_bunches << ","
            << s.input_photons << ","
            << s.hit_mirror << ","
            << s.hit_output_plane << ","
            << s.hit_camera << ","
            << s.accepted_camera << ","
            << s.lost_between_pixels << ","
            << s.unique_pixels.size() << ","
            << s.weighted_signal << ","
            << s.weighted_signal << ","
            << mean << ","
            << std::sqrt(var) << "\n";
    }
}

#ifdef LACT_HAS_HDF5
int hdf5PixelShapeCode(PixelShape shape)
{
    if (shape == PixelShape::Square) return 1;
    if (shape == PixelShape::Hexagonal) return 2;
    if (shape == PixelShape::Circular) return 3;
    return 0;
}

int hdf5FacetShapeCode(ApertureShape shape)
{
    if (shape == ApertureShape::Square) return 1;
    if (shape == ApertureShape::Hexagon) return 2;
    if (shape == ApertureShape::Circular) return 3;
    return 0;
}

void h5Check(herr_t status, const std::string& message)
{
    if (status < 0) {
        throw std::runtime_error("HDF5 write failed: " + message);
    }
}

void writeStringAttribute(hid_t object, const std::string& name, const std::string& value)
{
    hid_t space = H5Screate(H5S_SCALAR);
    if (space < 0) {
        throw std::runtime_error("HDF5 write failed: create scalar dataspace");
    }
    hid_t type = H5Tcopy(H5T_C_S1);
    if (type < 0) {
        H5Sclose(space);
        throw std::runtime_error("HDF5 write failed: create string type");
    }
    H5Tset_size(type, std::max<std::size_t>(1, value.size() + 1));
    H5Tset_strpad(type, H5T_STR_NULLTERM);
    hid_t attr = H5Acreate2(object, name.c_str(), type, space, H5P_DEFAULT, H5P_DEFAULT);
    if (attr < 0) {
        H5Tclose(type);
        H5Sclose(space);
        throw std::runtime_error("HDF5 write failed: create attribute " + name);
    }
    h5Check(H5Awrite(attr, type, value.c_str()), "write attribute " + name);
    H5Aclose(attr);
    H5Tclose(type);
    H5Sclose(space);
}

void writeStringDataset(hid_t group, const std::string& name, const std::string& value)
{
    hid_t space = H5Screate(H5S_SCALAR);
    if (space < 0) {
        throw std::runtime_error("HDF5 write failed: create scalar dataspace");
    }
    hid_t type = H5Tcopy(H5T_C_S1);
    if (type < 0) {
        H5Sclose(space);
        throw std::runtime_error("HDF5 write failed: create string type");
    }
    H5Tset_size(type, std::max<std::size_t>(1, value.size() + 1));
    H5Tset_strpad(type, H5T_STR_NULLTERM);
    hid_t ds = H5Dcreate2(group, name.c_str(), type, space,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (ds < 0) {
        H5Tclose(type);
        H5Sclose(space);
        throw std::runtime_error("HDF5 write failed: create dataset " + name);
    }
    h5Check(H5Dwrite(ds, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, value.c_str()),
            "write dataset " + name);
    H5Dclose(ds);
    H5Tclose(type);
    H5Sclose(space);
}

std::string readTextIfExists(const std::string& path)
{
    if (path.empty() || !std::filesystem::exists(path)) {
        return "";
    }
    std::ifstream ifs(path);
    std::ostringstream ss;
    ss << ifs.rdbuf();
    return ss.str();
}

template <typename Row>
void writeCompound1D(hid_t group,
                     const std::string& name,
                     hid_t type,
                     const std::vector<Row>& rows)
{
    hsize_t dims[1] = {rows.size()};
    hid_t space = H5Screate_simple(1, dims, nullptr);
    if (space < 0) {
        throw std::runtime_error("HDF5 write failed: create dataspace " + name);
    }
    hid_t ds = H5Dcreate2(group, name.c_str(), type, space,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (ds < 0) {
        H5Sclose(space);
        throw std::runtime_error("HDF5 write failed: create dataset " + name);
    }
    if (!rows.empty()) {
        h5Check(H5Dwrite(ds, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, rows.data()),
                "write dataset " + name);
    }
    H5Dclose(ds);
    H5Sclose(space);
}

template <typename T>
void writePlain1D(hid_t group,
                  const std::string& name,
                  hid_t type,
                  const std::vector<T>& values)
{
    hsize_t dims[1] = {values.size()};
    hid_t space = H5Screate_simple(1, dims, nullptr);
    if (space < 0) {
        throw std::runtime_error("HDF5 write failed: create dataspace " + name);
    }
    hid_t ds = H5Dcreate2(group, name.c_str(), type, space,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (ds < 0) {
        H5Sclose(space);
        throw std::runtime_error("HDF5 write failed: create dataset " + name);
    }
    if (!values.empty()) {
        h5Check(H5Dwrite(ds, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data()),
                "write dataset " + name);
    }
    H5Dclose(ds);
    H5Sclose(space);
}

template <typename T>
void writePlain2D(hid_t group,
                  const std::string& name,
                  hid_t type,
                  const std::vector<T>& values,
                  hsize_t n_rows,
                  hsize_t n_cols)
{
    hsize_t dims[2] = {n_rows, n_cols};
    hid_t space = H5Screate_simple(2, dims, nullptr);
    if (space < 0) {
        throw std::runtime_error("HDF5 write failed: create dataspace " + name);
    }
    hid_t ds = H5Dcreate2(group, name.c_str(), type, space,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (ds < 0) {
        H5Sclose(space);
        throw std::runtime_error("HDF5 write failed: create dataset " + name);
    }
    if (!values.empty()) {
        h5Check(H5Dwrite(ds, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data()),
                "write dataset " + name);
    }
    H5Dclose(ds);
    H5Sclose(space);
}

void writeNativeTraceHdf5(const CorsikaTraceOutputConfig& output_cfg,
                          const std::string& main_config_path,
                          const std::map<std::string, std::string>& cfg,
                          const ComponentConfigPaths& component_paths,
                          const SourceRuntimeConfig& source_runtime_cfg,
                          const TelescopeConfig& telescope_cfg,
                          const EventIOMetadata& metadata,
                          const CameraGeometry& camera,
                          const std::vector<MirrorFacet>& facets,
                          const ElectronicsConfig& electronics_cfg,
                          const NsbConfig& nsb_cfg,
                          const TriggerConfig& trigger_cfg,
                          const std::map<SummaryKey, TraceSummary>& summaries,
                          const std::map<PixelKey, PixelAccumulator>& pixels,
                          const std::vector<WhiteboardHdf5Row>& whiteboard_hits)
{
    const bool write_sparse =
        output_cfg.hdf5_storage == "sparse" || output_cfg.hdf5_storage == "both";
    const bool write_dense =
        output_cfg.hdf5_storage == "dense" || output_cfg.hdf5_storage == "both";

    const std::filesystem::path out_path(output_cfg.hdf5_path);
    if (out_path.has_parent_path()) {
        std::filesystem::create_directories(out_path.parent_path());
    }

    hid_t file = H5Fcreate(output_cfg.hdf5_path.c_str(), H5F_ACC_TRUNC,
                           H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0) {
        throw std::runtime_error("failed to create HDF5 file: " + output_cfg.hdf5_path);
    }

    try {
        writeStringAttribute(file, "format", "LACT_sim trace HDF5");
        writeStringAttribute(file, "format_version", "0.1-cpp");
        writeStringAttribute(file, "image_storage", output_cfg.hdf5_storage);
        writeStringAttribute(file, "hdf5_write_components",
                             output_cfg.hdf5_write_components ? "true" : "false");
        writeStringAttribute(file, "save_only_triggered",
                             output_cfg.save_only_triggered ? "true" : "false");
        writeStringAttribute(file, "event_id_mode", source_runtime_cfg.event_id_mode);
        writeStringAttribute(file, "source_eventio_path", source_runtime_cfg.eventio_path);

        hid_t config_group = H5Gcreate2(file, "config", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        writeStringAttribute(config_group, "main_config_path", main_config_path);
        for (const auto& kv : cfg) {
            writeStringAttribute(config_group, kv.first, kv.second);
        }
        const std::string main_text = readTextIfExists(main_config_path);
        if (!main_text.empty()) {
            writeStringDataset(config_group, "main_config_text", main_text);
        }
        hid_t components = H5Gcreate2(config_group, "components",
                                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        const std::vector<std::pair<std::string, std::string>> component_items = {
            {"telescope", component_paths.telescope},
            {"mirror", component_paths.mirror},
            {"source", component_paths.source},
            {"output", component_paths.output},
            {"camera", component_paths.camera},
            {"sipm", component_paths.sipm},
            {"electronics", component_paths.electronics},
            {"efficiency", component_paths.efficiency},
            {"atmosphere", component_paths.atmosphere},
            {"error", component_paths.error},
        };
        for (const auto& item : component_items) {
            if (item.second.empty()) {
                continue;
            }
            writeStringAttribute(components, item.first, item.second);
            const std::string text = readTextIfExists(item.second);
            if (!text.empty()) {
                writeStringDataset(components, item.first + "_text", text);
            }
        }
        H5Gclose(components);
        H5Gclose(config_group);

        hid_t metadata_group = H5Gcreate2(file, "metadata",
                                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        hid_t coordinates_group = H5Gcreate2(metadata_group, "coordinates",
                                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        writeStringAttribute(coordinates_group, "array_position_frame",
                             "CORSIKA IACT NWU horizontal frame");
        writeStringAttribute(coordinates_group, "array_x_m",
                             "CORSIKA magnetic-North-positive telescope position coordinate");
        writeStringAttribute(coordinates_group, "array_y_m",
                             "West-positive telescope position coordinate");
        writeStringAttribute(coordinates_group, "array_z_m",
                             "Up-positive telescope position coordinate");
        writeStringAttribute(coordinates_group, "pointing_az_deg",
                             "CORSIKA magnetic-North-to-East azimuth; 0=+array_x, 90=East/-array_y");
        writeStringAttribute(coordinates_group, "pointing_el_deg",
                             "Elevation above local horizon; zenith angle = 90 - elevation");
        writeStringAttribute(coordinates_group, "eventio_photon_frame",
                             source_runtime_cfg.eventio_coordinate_frame);
        writeStringAttribute(coordinates_group, "eventio_corsika_iact_positions",
                             "Photon bunch x/y/z are telescope-relative CORSIKA IACT coordinates before rotation to telescope-local optics");
        writeStringAttribute(coordinates_group, "eventio_teloff_core_note",
                             "MC_TELOFF is detector-array offset with respect to shower core; /events/corsika stores core = -MC_TELOFF in NWU coordinates");
        writeStringAttribute(coordinates_group, "eventio_telpos_note",
                             "Telescope positions are hessio MC_TELPOS detector coordinates; array_z_up_m may include the detector sphere/radius convention used by the producer");
        writeStringAttribute(coordinates_group, "camera_plane_coordinates",
                             "camera x/y are output-plane u/v coordinates, normally u=local +x and v=local +y");
        H5Gclose(coordinates_group);

        hid_t electronics_group = H5Gcreate2(metadata_group, "electronics",
                                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        writeStringAttribute(electronics_group, "model", "integrated_pe_placeholder");
        writeStringAttribute(electronics_group, "pe_conversion",
                             factorDescription(electronics_cfg.pe_conversion));
        H5Gclose(electronics_group);

        hid_t nsb_group = H5Gcreate2(metadata_group, "nsb",
                                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        writeStringAttribute(nsb_group, "enabled", nsb_cfg.enabled ? "true" : "false");
        writeStringAttribute(nsb_group, "model", nsb_cfg.model);
        writeStringAttribute(nsb_group, "rate_pe_per_ns_per_pixel",
                             doubleToString(nsb_cfg.rate_pe_per_ns_per_pixel));
        writeStringAttribute(nsb_group, "window_ns", doubleToString(nsb_cfg.window_ns));
        writeStringAttribute(nsb_group, "seed", intToString(nsb_cfg.seed));
        H5Gclose(nsb_group);

        hid_t trigger_group_meta = H5Gcreate2(metadata_group, "trigger",
                                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        writeStringAttribute(trigger_group_meta, "enabled",
                             trigger_cfg.enabled ? "true" : "false");
        writeStringAttribute(trigger_group_meta, "model", "simple_multiplicity");
        writeStringAttribute(trigger_group_meta, "pixel_threshold_pe",
                             doubleToString(trigger_cfg.pixel_threshold_pe));
        writeStringAttribute(trigger_group_meta, "camera_multiplicity",
                             intToString(trigger_cfg.camera_multiplicity));
        writeStringAttribute(trigger_group_meta, "array_multiplicity",
                             intToString(trigger_cfg.array_multiplicity));
        writeStringAttribute(trigger_group_meta, "coincidence_window_ns",
                             doubleToString(trigger_cfg.coincidence_window_ns));
        H5Gclose(trigger_group_meta);
        H5Gclose(metadata_group);

        struct CameraRow {
            std::int32_t pixel_id;
            float x_m;
            float y_m;
            float size_m;
            std::int16_t shape_code;
        };
        std::vector<CameraRow> camera_rows;
        camera_rows.reserve(camera.size());
        for (const auto& pixel : camera.pixels()) {
            camera_rows.push_back(CameraRow{
                static_cast<std::int32_t>(pixel.id),
                static_cast<float>(pixel.center.x),
                static_cast<float>(pixel.center.y),
                static_cast<float>(pixel.size),
                static_cast<std::int16_t>(hdf5PixelShapeCode(pixel.shape)),
            });
        }
        hid_t camera_group = H5Gcreate2(file, "camera", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        writeStringAttribute(camera_group, "shape_code_map",
                             "0=unknown,1=square,2=hexagon,3=circular");
        hid_t camera_type = H5Tcreate(H5T_COMPOUND, sizeof(CameraRow));
        H5Tinsert(camera_type, "pixel_id", HOFFSET(CameraRow, pixel_id), H5T_NATIVE_INT32);
        H5Tinsert(camera_type, "x_m", HOFFSET(CameraRow, x_m), H5T_NATIVE_FLOAT);
        H5Tinsert(camera_type, "y_m", HOFFSET(CameraRow, y_m), H5T_NATIVE_FLOAT);
        H5Tinsert(camera_type, "size_m", HOFFSET(CameraRow, size_m), H5T_NATIVE_FLOAT);
        H5Tinsert(camera_type, "shape_code", HOFFSET(CameraRow, shape_code), H5T_NATIVE_INT16);
        writeCompound1D(camera_group, "pixels", camera_type, camera_rows);
        H5Tclose(camera_type);
        H5Gclose(camera_group);

        struct FacetRow {
            std::int32_t mirror_id;
            float center_x_m;
            float center_y_m;
            float center_z_m;
            float normal_x;
            float normal_y;
            float normal_z;
            float radius_of_curvature_m;
            float size1_m;
            float size2_m;
            float aperture_rotation_rad;
            std::int16_t shape_code;
        };
        std::vector<FacetRow> facet_rows;
        facet_rows.reserve(facets.size());
        for (const auto& facet : facets) {
            facet_rows.push_back(FacetRow{
                static_cast<std::int32_t>(facet.id),
                static_cast<float>(facet.center.x),
                static_cast<float>(facet.center.y),
                static_cast<float>(facet.center.z),
                static_cast<float>(facet.normal.x),
                static_cast<float>(facet.normal.y),
                static_cast<float>(facet.normal.z),
                static_cast<float>(facet.radius_of_curvature),
                static_cast<float>(facet.size1),
                static_cast<float>(facet.size2),
                static_cast<float>(facet.aperture_rotation_rad),
                static_cast<std::int16_t>(hdf5FacetShapeCode(facet.aperture_shape)),
            });
        }
        hid_t mirror_group = H5Gcreate2(file, "mirrors", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        writeStringAttribute(mirror_group, "shape_code_map",
                             "0=unknown,1=square,2=hexagon,3=circular");
        hid_t facet_type = H5Tcreate(H5T_COMPOUND, sizeof(FacetRow));
        H5Tinsert(facet_type, "mirror_id", HOFFSET(FacetRow, mirror_id), H5T_NATIVE_INT32);
        H5Tinsert(facet_type, "center_x_m", HOFFSET(FacetRow, center_x_m), H5T_NATIVE_FLOAT);
        H5Tinsert(facet_type, "center_y_m", HOFFSET(FacetRow, center_y_m), H5T_NATIVE_FLOAT);
        H5Tinsert(facet_type, "center_z_m", HOFFSET(FacetRow, center_z_m), H5T_NATIVE_FLOAT);
        H5Tinsert(facet_type, "normal_x", HOFFSET(FacetRow, normal_x), H5T_NATIVE_FLOAT);
        H5Tinsert(facet_type, "normal_y", HOFFSET(FacetRow, normal_y), H5T_NATIVE_FLOAT);
        H5Tinsert(facet_type, "normal_z", HOFFSET(FacetRow, normal_z), H5T_NATIVE_FLOAT);
        H5Tinsert(facet_type, "radius_of_curvature_m",
                  HOFFSET(FacetRow, radius_of_curvature_m), H5T_NATIVE_FLOAT);
        H5Tinsert(facet_type, "size1_m", HOFFSET(FacetRow, size1_m), H5T_NATIVE_FLOAT);
        H5Tinsert(facet_type, "size2_m", HOFFSET(FacetRow, size2_m), H5T_NATIVE_FLOAT);
        H5Tinsert(facet_type, "aperture_rotation_rad",
                  HOFFSET(FacetRow, aperture_rotation_rad), H5T_NATIVE_FLOAT);
        H5Tinsert(facet_type, "shape_code", HOFFSET(FacetRow, shape_code), H5T_NATIVE_INT16);
        writeCompound1D(mirror_group, "facets", facet_type, facet_rows);
        H5Tclose(facet_type);
        H5Gclose(mirror_group);

        struct TelescopeRow {
            std::int32_t telescope_id;
            double x_m;
            double y_m;
            double z_m;
            double array_x_north_m;
            double array_y_west_m;
            double array_z_up_m;
            double radius_m;
            double pointing_az_deg;
            double pointing_el_deg;
            double focal_length_m;
        };
        std::vector<TelescopeRow> telescope_rows;
        if (!metadata.telescopes.empty()) {
            telescope_rows.reserve(metadata.telescopes.size());
            for (const auto& tel : metadata.telescopes) {
                telescope_rows.push_back(TelescopeRow{
                    static_cast<std::int32_t>(tel.telescope_id),
                    tel.x_m,
                    tel.y_m,
                    tel.z_m,
                    tel.x_m,
                    tel.y_m,
                    tel.z_m,
                    tel.radius_m,
                    telescope_cfg.pointing_az_deg,
                    telescope_cfg.pointing_el_deg,
                    telescope_cfg.focal_length_m,
                });
            }
        } else {
            telescope_rows.push_back(TelescopeRow{
                static_cast<std::int32_t>(telescope_cfg.id),
                telescope_cfg.position_m.x,
                telescope_cfg.position_m.y,
                telescope_cfg.position_m.z,
                telescope_cfg.position_m.x,
                telescope_cfg.position_m.y,
                telescope_cfg.position_m.z,
                0.0,
                telescope_cfg.pointing_az_deg,
                telescope_cfg.pointing_el_deg,
                telescope_cfg.focal_length_m,
            });
        }
        hid_t telescope_group = H5Gcreate2(file, "telescopes",
                                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        writeStringAttribute(telescope_group, "coordinate_frame",
                             "CORSIKA IACT horizontal frame");
        writeStringAttribute(telescope_group, "x_m_compat",
                             "same as array_x_north_m; kept for compatibility");
        writeStringAttribute(telescope_group, "y_m_compat",
                             "same as array_y_west_m; kept for compatibility");
        writeStringAttribute(telescope_group, "z_m_compat",
                             "same as array_z_up_m; kept for compatibility");
        writeStringAttribute(telescope_group, "pointing_convention",
                             "azimuth North-to-East, elevation above horizon");
        hid_t telescope_type = H5Tcreate(H5T_COMPOUND, sizeof(TelescopeRow));
        H5Tinsert(telescope_type, "telescope_id",
                  HOFFSET(TelescopeRow, telescope_id), H5T_NATIVE_INT32);
        H5Tinsert(telescope_type, "x_m", HOFFSET(TelescopeRow, x_m), H5T_NATIVE_DOUBLE);
        H5Tinsert(telescope_type, "y_m", HOFFSET(TelescopeRow, y_m), H5T_NATIVE_DOUBLE);
        H5Tinsert(telescope_type, "z_m", HOFFSET(TelescopeRow, z_m), H5T_NATIVE_DOUBLE);
        H5Tinsert(telescope_type, "array_x_north_m",
                  HOFFSET(TelescopeRow, array_x_north_m), H5T_NATIVE_DOUBLE);
        H5Tinsert(telescope_type, "array_y_west_m",
                  HOFFSET(TelescopeRow, array_y_west_m), H5T_NATIVE_DOUBLE);
        H5Tinsert(telescope_type, "array_z_up_m",
                  HOFFSET(TelescopeRow, array_z_up_m), H5T_NATIVE_DOUBLE);
        H5Tinsert(telescope_type, "radius_m",
                  HOFFSET(TelescopeRow, radius_m), H5T_NATIVE_DOUBLE);
        H5Tinsert(telescope_type, "pointing_az_deg",
                  HOFFSET(TelescopeRow, pointing_az_deg), H5T_NATIVE_DOUBLE);
        H5Tinsert(telescope_type, "pointing_el_deg",
                  HOFFSET(TelescopeRow, pointing_el_deg), H5T_NATIVE_DOUBLE);
        H5Tinsert(telescope_type, "focal_length_m",
                  HOFFSET(TelescopeRow, focal_length_m), H5T_NATIVE_DOUBLE);
        writeCompound1D(telescope_group, "table", telescope_type, telescope_rows);
        H5Tclose(telescope_type);
        H5Gclose(telescope_group);

        std::set<SummaryKey> image_keys;
        for (const auto& kv : summaries) {
            image_keys.insert(kv.first);
        }
        for (const auto& kv : pixels) {
            image_keys.insert({std::get<0>(kv.first), std::get<1>(kv.first)});
        }

        struct SparsePixelRow {
            std::int32_t pixel_id;
            std::int32_t photon_count;
            float pe;
            float signal;
            float time_mean_ns;
            float time_rms_ns;
        };
        struct ImageIndexRow {
            std::int32_t image_index;
            std::int64_t event_id;
            std::int32_t telescope_id;
            std::int64_t start;
            std::int32_t count;
            double total_photons;
            double total_pe;
            double total_signal;
            float time_mean_ns;
            float time_rms_ns;
        };
        std::vector<SparsePixelRow> sparse_rows;
        std::vector<ImageIndexRow> image_rows;
        image_rows.reserve(image_keys.size());

        std::int32_t image_index = 0;
        for (const auto& key : image_keys) {
            const auto start = static_cast<std::int64_t>(sparse_rows.size());
            const int event_id = key.first;
            const int telescope_id = key.second;
            double total_signal = 0.0;
            double total_pe = 0.0;
            double total_photons = 0.0;
            for (const auto& kv : pixels) {
                if (std::get<0>(kv.first) != event_id ||
                    std::get<1>(kv.first) != telescope_id) {
                    continue;
                }
                const auto& p = kv.second;
                const double mean = p.signal > 0.0 ? p.time_sum / p.signal : 0.0;
                const double var = p.signal > 0.0
                    ? std::max(0.0, p.time2_sum / p.signal - mean * mean)
                    : 0.0;
                sparse_rows.push_back(SparsePixelRow{
                    static_cast<std::int32_t>(p.pixel_id),
                    static_cast<std::int32_t>(p.photon_count),
                    static_cast<float>(p.pe),
                    static_cast<float>(p.signal),
                    static_cast<float>(mean),
                    static_cast<float>(std::sqrt(var)),
                });
                total_signal += p.signal;
                total_pe += p.pe;
                total_photons += static_cast<double>(p.photon_count);
            }

            float mean = 0.0f;
            float rms = 0.0f;
            auto summary_it = summaries.find(key);
            if (summary_it != summaries.end()) {
                const auto& s = summary_it->second;
                total_signal = s.weighted_signal;
                total_pe = s.weighted_signal;
                total_photons = static_cast<double>(s.hit_camera);
                if (s.weighted_signal > 0.0) {
                    const double m = s.weighted_time_sum / s.weighted_signal;
                    const double v = std::max(0.0, s.weighted_time2_sum / s.weighted_signal - m * m);
                    mean = static_cast<float>(m);
                    rms = static_cast<float>(std::sqrt(v));
                }
            }
            const auto count = static_cast<std::int32_t>(
                static_cast<std::int64_t>(sparse_rows.size()) - start);
            image_rows.push_back(ImageIndexRow{
                image_index++,
                static_cast<std::int64_t>(event_id),
                static_cast<std::int32_t>(telescope_id),
                start,
                count,
                total_photons,
                total_pe,
                total_signal,
                mean,
                rms,
            });
        }

        std::vector<std::int32_t> pixel_id_axis;
        std::map<std::int32_t, std::size_t> pixel_to_col;
        std::vector<float> dense_signal;
        std::vector<float> dense_pe;
        std::vector<float> dense_cherenkov_pe;
        std::vector<float> dense_nsb_pe;
        std::vector<std::int32_t> dense_photon_count;
        const bool have_dense_images = write_dense && !camera_rows.empty();
        if (have_dense_images) {
            pixel_id_axis.reserve(camera_rows.size());
            for (std::size_t i = 0; i < camera_rows.size(); ++i) {
                pixel_id_axis.push_back(camera_rows[i].pixel_id);
                pixel_to_col[camera_rows[i].pixel_id] = i;
            }

            const std::size_t n_images = image_rows.size();
            const std::size_t n_pixels = camera_rows.size();
            dense_signal.assign(n_images * n_pixels, 0.0f);
            dense_pe.assign(n_images * n_pixels, 0.0f);
            dense_photon_count.assign(n_images * n_pixels, 0);
            for (const auto& image : image_rows) {
                const std::size_t row = static_cast<std::size_t>(image.image_index);
                const std::int64_t begin = image.start;
                const std::int64_t end = image.start + image.count;
                for (std::int64_t i = begin; i < end; ++i) {
                    const auto& pixel = sparse_rows[static_cast<std::size_t>(i)];
                    const auto col_it = pixel_to_col.find(pixel.pixel_id);
                    if (col_it == pixel_to_col.end()) {
                        continue;
                    }
                    const std::size_t index = row * n_pixels + col_it->second;
                    dense_signal[index] = pixel.signal;
                    dense_pe[index] = pixel.pe;
                    dense_photon_count[index] = pixel.photon_count;
                }
            }

            dense_cherenkov_pe = dense_pe;
            dense_nsb_pe.assign(n_images * n_pixels, 0.0f);
            if (nsb_cfg.enabled && nsb_cfg.rate_pe_per_ns_per_pixel > 0.0 &&
                nsb_cfg.window_ns > 0.0) {
                std::mt19937_64 rng(nsb_cfg.seed);
                std::poisson_distribution<int> poisson(
                    nsb_cfg.rate_pe_per_ns_per_pixel * nsb_cfg.window_ns);
                for (std::size_t i = 0; i < dense_nsb_pe.size(); ++i) {
                    const float nsb_pe = static_cast<float>(poisson(rng));
                    dense_nsb_pe[i] = nsb_pe;
                    dense_pe[i] += nsb_pe;
                    dense_signal[i] += nsb_pe;
                }
            }

            for (auto& image : image_rows) {
                const std::size_t row = static_cast<std::size_t>(image.image_index);
                double total_pe = 0.0;
                double total_signal = 0.0;
                for (std::size_t col = 0; col < n_pixels; ++col) {
                    const std::size_t index = row * n_pixels + col;
                    total_pe += dense_pe[index];
                    total_signal += dense_signal[index];
                }
                image.total_pe = total_pe;
                image.total_signal = total_signal;
            }
        }

        struct TelescopeTriggerRow {
            std::int64_t event_id;
            std::int32_t telescope_id;
            std::int8_t triggered;
            std::int32_t n_pixels_above_threshold;
            double total_pe;
            float trigger_time_ns;
        };
        struct ArrayTriggerRow {
            std::int64_t event_id;
            std::int8_t array_triggered;
            std::int32_t n_triggered_telescopes;
        };
        std::vector<TelescopeTriggerRow> telescope_trigger_rows;
        std::map<int, int> triggered_telescopes_by_event;
        telescope_trigger_rows.reserve(image_rows.size());
        for (const auto& image : image_rows) {
            int above_threshold = 0;
            double total_pe = image.total_pe;
            if (have_dense_images) {
                const std::size_t row = static_cast<std::size_t>(image.image_index);
                const std::size_t n_pixels = camera_rows.size();
                total_pe = 0.0;
                for (std::size_t col = 0; col < n_pixels; ++col) {
                    const double pe = dense_pe[row * n_pixels + col];
                    total_pe += pe;
                    if (pe >= trigger_cfg.pixel_threshold_pe) {
                        ++above_threshold;
                    }
                }
            } else {
                const std::int64_t begin = image.start;
                const std::int64_t end = image.start + image.count;
                for (std::int64_t i = begin; i < end; ++i) {
                    if (sparse_rows[static_cast<std::size_t>(i)].pe >=
                        trigger_cfg.pixel_threshold_pe) {
                        ++above_threshold;
                    }
                }
            }
            const bool telescope_triggered =
                trigger_cfg.enabled && above_threshold >= trigger_cfg.camera_multiplicity;
            if (telescope_triggered) {
                triggered_telescopes_by_event[static_cast<int>(image.event_id)] += 1;
            }
            telescope_trigger_rows.push_back(TelescopeTriggerRow{
                image.event_id,
                image.telescope_id,
                static_cast<std::int8_t>(telescope_triggered ? 1 : 0),
                static_cast<std::int32_t>(above_threshold),
                total_pe,
                image.time_mean_ns,
            });
        }
        std::set<int> trigger_event_ids;
        for (const auto& key : image_keys) {
            trigger_event_ids.insert(key.first);
        }
        std::vector<ArrayTriggerRow> array_trigger_rows;
        array_trigger_rows.reserve(trigger_event_ids.size());
        for (const int event_id : trigger_event_ids) {
            const int n_triggered = triggered_telescopes_by_event[event_id];
            array_trigger_rows.push_back(ArrayTriggerRow{
                static_cast<std::int64_t>(event_id),
                static_cast<std::int8_t>(
                    trigger_cfg.enabled && n_triggered >= trigger_cfg.array_multiplicity ? 1 : 0),
                static_cast<std::int32_t>(n_triggered),
            });
        }

        struct EventRow {
            std::int32_t event_index;
            std::int64_t event_id;
        };
        struct CorsikaEventRow {
            std::int64_t event_id;
            std::int32_t shower_event_id;
            std::int32_t array_id;
            std::int32_t primary_type;
            double energy_gev;
            double theta_deg;
            double phi_deg;
            double azimuth_north_to_east_deg;
            double core_x_north_m;
            double core_y_west_m;
            double array_rotation_deg;
        };
        struct CorsikaShowerRow {
            std::int32_t shower_event_id;
            std::int32_t primary_type;
            double energy_gev;
            double theta_deg;
            double phi_deg;
            double azimuth_north_to_east_deg;
            double core_x_north_m;
            double core_y_west_m;
            double array_rotation_deg;
        };
        std::vector<EventRow> event_rows;
        std::vector<CorsikaEventRow> corsika_event_rows;
        std::set<int> event_ids;
        for (const auto& key : image_keys) {
            event_ids.insert(key.first);
        }
        int event_index = 0;
        for (const int event_id : event_ids) {
            event_rows.push_back(EventRow{
                event_index++,
                static_cast<std::int64_t>(event_id),
            });
            const int shower_event =
                showerEventFromOutputEvent(event_id, source_runtime_cfg.event_id_mode);
            const int array_id =
                arrayIdFromOutputEvent(event_id, source_runtime_cfg.event_id_mode);
            auto event_it = std::find_if(
                metadata.events.begin(), metadata.events.end(),
                [shower_event](const EventIOEventHeader& event) {
                    return event.shower_event_id == shower_event;
                });
            if (event_it != metadata.events.end()) {
                const OutputEventMetadata event_meta =
                    outputEventMetadata(event_id,
                                        source_runtime_cfg.event_id_mode,
                                        metadata);
                corsika_event_rows.push_back(CorsikaEventRow{
                    static_cast<std::int64_t>(event_id),
                    static_cast<std::int32_t>(shower_event),
                    static_cast<std::int32_t>(array_id),
                    static_cast<std::int32_t>(event_it->primary_type),
                    event_it->energy_gev,
                    event_it->theta_deg,
                    event_it->phi_deg,
                    event_it->azimuth_north_to_east_deg,
                    event_meta.core_x_north_m,
                    event_meta.core_y_west_m,
                    event_it->array_rotation_deg,
                });
            }
        }
        hid_t events_group = H5Gcreate2(file, "events", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        hid_t event_type = H5Tcreate(H5T_COMPOUND, sizeof(EventRow));
        H5Tinsert(event_type, "event_index", HOFFSET(EventRow, event_index), H5T_NATIVE_INT32);
        H5Tinsert(event_type, "event_id", HOFFSET(EventRow, event_id), H5T_NATIVE_INT64);
        writeCompound1D(events_group, "table", event_type, event_rows);
        H5Tclose(event_type);
        if (!corsika_event_rows.empty()) {
            hid_t corsika_event_type =
                H5Tcreate(H5T_COMPOUND, sizeof(CorsikaEventRow));
            H5Tinsert(corsika_event_type, "event_id",
                      HOFFSET(CorsikaEventRow, event_id), H5T_NATIVE_INT64);
            H5Tinsert(corsika_event_type, "shower_event_id",
                      HOFFSET(CorsikaEventRow, shower_event_id), H5T_NATIVE_INT32);
            H5Tinsert(corsika_event_type, "array_id",
                      HOFFSET(CorsikaEventRow, array_id), H5T_NATIVE_INT32);
            H5Tinsert(corsika_event_type, "primary_type",
                      HOFFSET(CorsikaEventRow, primary_type), H5T_NATIVE_INT32);
            H5Tinsert(corsika_event_type, "energy_gev",
                      HOFFSET(CorsikaEventRow, energy_gev), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_event_type, "theta_deg",
                      HOFFSET(CorsikaEventRow, theta_deg), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_event_type, "phi_deg",
                      HOFFSET(CorsikaEventRow, phi_deg), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_event_type, "azimuth_north_to_east_deg",
                      HOFFSET(CorsikaEventRow, azimuth_north_to_east_deg),
                      H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_event_type, "core_x_north_m",
                      HOFFSET(CorsikaEventRow, core_x_north_m), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_event_type, "core_y_west_m",
                      HOFFSET(CorsikaEventRow, core_y_west_m), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_event_type, "array_rotation_deg",
                      HOFFSET(CorsikaEventRow, array_rotation_deg), H5T_NATIVE_DOUBLE);
            writeCompound1D(events_group, "corsika", corsika_event_type,
                            corsika_event_rows);
            H5Tclose(corsika_event_type);
        }
        if (!metadata.events.empty()) {
            std::vector<CorsikaShowerRow> corsika_shower_rows;
            corsika_shower_rows.reserve(metadata.events.size());
            for (const auto& event : metadata.events) {
                corsika_shower_rows.push_back(CorsikaShowerRow{
                    static_cast<std::int32_t>(event.shower_event_id),
                    static_cast<std::int32_t>(event.primary_type),
                    event.energy_gev,
                    event.theta_deg,
                    event.phi_deg,
                    event.azimuth_north_to_east_deg,
                    event.core_x_m,
                    event.core_y_m,
                    event.array_rotation_deg,
                });
            }
            hid_t corsika_shower_type =
                H5Tcreate(H5T_COMPOUND, sizeof(CorsikaShowerRow));
            H5Tinsert(corsika_shower_type, "shower_event_id",
                      HOFFSET(CorsikaShowerRow, shower_event_id), H5T_NATIVE_INT32);
            H5Tinsert(corsika_shower_type, "primary_type",
                      HOFFSET(CorsikaShowerRow, primary_type), H5T_NATIVE_INT32);
            H5Tinsert(corsika_shower_type, "energy_gev",
                      HOFFSET(CorsikaShowerRow, energy_gev), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_shower_type, "theta_deg",
                      HOFFSET(CorsikaShowerRow, theta_deg), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_shower_type, "phi_deg",
                      HOFFSET(CorsikaShowerRow, phi_deg), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_shower_type, "azimuth_north_to_east_deg",
                      HOFFSET(CorsikaShowerRow, azimuth_north_to_east_deg),
                      H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_shower_type, "core_x_north_m",
                      HOFFSET(CorsikaShowerRow, core_x_north_m), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_shower_type, "core_y_west_m",
                      HOFFSET(CorsikaShowerRow, core_y_west_m), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_shower_type, "array_rotation_deg",
                      HOFFSET(CorsikaShowerRow, array_rotation_deg), H5T_NATIVE_DOUBLE);
            writeCompound1D(events_group, "corsika_showers", corsika_shower_type,
                            corsika_shower_rows);
            H5Tclose(corsika_shower_type);
        }
        H5Gclose(events_group);

        hid_t images_group = H5Gcreate2(file, "images", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        hid_t image_type = H5Tcreate(H5T_COMPOUND, sizeof(ImageIndexRow));
        H5Tinsert(image_type, "image_index",
                  HOFFSET(ImageIndexRow, image_index), H5T_NATIVE_INT32);
        H5Tinsert(image_type, "event_id", HOFFSET(ImageIndexRow, event_id), H5T_NATIVE_INT64);
        H5Tinsert(image_type, "telescope_id",
                  HOFFSET(ImageIndexRow, telescope_id), H5T_NATIVE_INT32);
        H5Tinsert(image_type, "start", HOFFSET(ImageIndexRow, start), H5T_NATIVE_INT64);
        H5Tinsert(image_type, "count", HOFFSET(ImageIndexRow, count), H5T_NATIVE_INT32);
        H5Tinsert(image_type, "total_photons",
                  HOFFSET(ImageIndexRow, total_photons), H5T_NATIVE_DOUBLE);
        H5Tinsert(image_type, "total_pe", HOFFSET(ImageIndexRow, total_pe), H5T_NATIVE_DOUBLE);
        H5Tinsert(image_type, "total_signal",
                  HOFFSET(ImageIndexRow, total_signal), H5T_NATIVE_DOUBLE);
        H5Tinsert(image_type, "time_mean_ns",
                  HOFFSET(ImageIndexRow, time_mean_ns), H5T_NATIVE_FLOAT);
        H5Tinsert(image_type, "time_rms_ns",
                  HOFFSET(ImageIndexRow, time_rms_ns), H5T_NATIVE_FLOAT);
        writeCompound1D(images_group, "index", image_type, image_rows);
        H5Tclose(image_type);

        if (write_sparse) {
            hid_t sparse_group = H5Gcreate2(images_group, "sparse",
                                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            hid_t sparse_type = H5Tcreate(H5T_COMPOUND, sizeof(SparsePixelRow));
            H5Tinsert(sparse_type, "pixel_id",
                      HOFFSET(SparsePixelRow, pixel_id), H5T_NATIVE_INT32);
            H5Tinsert(sparse_type, "photon_count",
                      HOFFSET(SparsePixelRow, photon_count), H5T_NATIVE_INT32);
            H5Tinsert(sparse_type, "pe", HOFFSET(SparsePixelRow, pe), H5T_NATIVE_FLOAT);
            H5Tinsert(sparse_type, "signal", HOFFSET(SparsePixelRow, signal), H5T_NATIVE_FLOAT);
            H5Tinsert(sparse_type, "time_mean_ns",
                      HOFFSET(SparsePixelRow, time_mean_ns), H5T_NATIVE_FLOAT);
            H5Tinsert(sparse_type, "time_rms_ns",
                      HOFFSET(SparsePixelRow, time_rms_ns), H5T_NATIVE_FLOAT);
            writeCompound1D(sparse_group, "pixels", sparse_type, sparse_rows);
            H5Tclose(sparse_type);
            H5Gclose(sparse_group);
        }

        if (have_dense_images) {
            const std::size_t n_images = image_rows.size();
            const std::size_t n_pixels = camera_rows.size();
            hid_t dense_group = H5Gcreate2(images_group, "dense",
                                           H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            writePlain1D(dense_group, "pixel_id_axis", H5T_NATIVE_INT32, pixel_id_axis);
            writePlain2D(dense_group,
                         "signal",
                         H5T_NATIVE_FLOAT,
                         dense_signal,
                         static_cast<hsize_t>(n_images),
                         static_cast<hsize_t>(n_pixels));
            writePlain2D(dense_group,
                         "pe",
                         H5T_NATIVE_FLOAT,
                         dense_pe,
                         static_cast<hsize_t>(n_images),
                         static_cast<hsize_t>(n_pixels));
            writePlain2D(dense_group,
                         "photon_count",
                         H5T_NATIVE_INT32,
                         dense_photon_count,
                         static_cast<hsize_t>(n_images),
                         static_cast<hsize_t>(n_pixels));
            if (output_cfg.hdf5_write_components) {
                writePlain2D(dense_group,
                             "cherenkov_pe",
                             H5T_NATIVE_FLOAT,
                             dense_cherenkov_pe,
                             static_cast<hsize_t>(n_images),
                             static_cast<hsize_t>(n_pixels));
                writePlain2D(dense_group,
                             "nsb_pe",
                             H5T_NATIVE_FLOAT,
                             dense_nsb_pe,
                             static_cast<hsize_t>(n_images),
                             static_cast<hsize_t>(n_pixels));
            }
            H5Gclose(dense_group);
        }
        H5Gclose(images_group);

        hid_t trigger_group = H5Gcreate2(file, "trigger",
                                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        hid_t telescope_trigger_type =
            H5Tcreate(H5T_COMPOUND, sizeof(TelescopeTriggerRow));
        H5Tinsert(telescope_trigger_type, "event_id",
                  HOFFSET(TelescopeTriggerRow, event_id), H5T_NATIVE_INT64);
        H5Tinsert(telescope_trigger_type, "telescope_id",
                  HOFFSET(TelescopeTriggerRow, telescope_id), H5T_NATIVE_INT32);
        H5Tinsert(telescope_trigger_type, "triggered",
                  HOFFSET(TelescopeTriggerRow, triggered), H5T_NATIVE_SCHAR);
        H5Tinsert(telescope_trigger_type, "n_pixels_above_threshold",
                  HOFFSET(TelescopeTriggerRow, n_pixels_above_threshold),
                  H5T_NATIVE_INT32);
        H5Tinsert(telescope_trigger_type, "total_pe",
                  HOFFSET(TelescopeTriggerRow, total_pe), H5T_NATIVE_DOUBLE);
        H5Tinsert(telescope_trigger_type, "trigger_time_ns",
                  HOFFSET(TelescopeTriggerRow, trigger_time_ns), H5T_NATIVE_FLOAT);
        writeCompound1D(trigger_group, "telescope", telescope_trigger_type,
                        telescope_trigger_rows);
        H5Tclose(telescope_trigger_type);

        hid_t array_trigger_type = H5Tcreate(H5T_COMPOUND, sizeof(ArrayTriggerRow));
        H5Tinsert(array_trigger_type, "event_id",
                  HOFFSET(ArrayTriggerRow, event_id), H5T_NATIVE_INT64);
        H5Tinsert(array_trigger_type, "array_triggered",
                  HOFFSET(ArrayTriggerRow, array_triggered), H5T_NATIVE_SCHAR);
        H5Tinsert(array_trigger_type, "n_triggered_telescopes",
                  HOFFSET(ArrayTriggerRow, n_triggered_telescopes), H5T_NATIVE_INT32);
        writeCompound1D(trigger_group, "array", array_trigger_type, array_trigger_rows);
        H5Tclose(array_trigger_type);
        H5Gclose(trigger_group);

        if (!whiteboard_hits.empty()) {
            hid_t whiteboard_group = H5Gcreate2(file, "whiteboard",
                                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            hid_t hit_type = H5Tcreate(H5T_COMPOUND, sizeof(WhiteboardHdf5Row));
            H5Tinsert(hit_type, "event_id",
                      HOFFSET(WhiteboardHdf5Row, event_id), H5T_NATIVE_INT64);
            H5Tinsert(hit_type, "telescope_id",
                      HOFFSET(WhiteboardHdf5Row, telescope_id), H5T_NATIVE_INT32);
            H5Tinsert(hit_type, "photon_index",
                      HOFFSET(WhiteboardHdf5Row, photon_index), H5T_NATIVE_INT64);
            H5Tinsert(hit_type, "mirror_id",
                      HOFFSET(WhiteboardHdf5Row, mirror_id), H5T_NATIVE_INT32);
            H5Tinsert(hit_type, "surface_x_m",
                      HOFFSET(WhiteboardHdf5Row, surface_x_m), H5T_NATIVE_FLOAT);
            H5Tinsert(hit_type, "surface_y_m",
                      HOFFSET(WhiteboardHdf5Row, surface_y_m), H5T_NATIVE_FLOAT);
            H5Tinsert(hit_type, "surface_z_m",
                      HOFFSET(WhiteboardHdf5Row, surface_z_m), H5T_NATIVE_FLOAT);
            H5Tinsert(hit_type, "u_m", HOFFSET(WhiteboardHdf5Row, u_m), H5T_NATIVE_FLOAT);
            H5Tinsert(hit_type, "v_m", HOFFSET(WhiteboardHdf5Row, v_m), H5T_NATIVE_FLOAT);
            H5Tinsert(hit_type, "dir_x",
                      HOFFSET(WhiteboardHdf5Row, dir_x), H5T_NATIVE_FLOAT);
            H5Tinsert(hit_type, "dir_y",
                      HOFFSET(WhiteboardHdf5Row, dir_y), H5T_NATIVE_FLOAT);
            H5Tinsert(hit_type, "dir_z",
                      HOFFSET(WhiteboardHdf5Row, dir_z), H5T_NATIVE_FLOAT);
            H5Tinsert(hit_type, "time_ns",
                      HOFFSET(WhiteboardHdf5Row, time_ns), H5T_NATIVE_FLOAT);
            H5Tinsert(hit_type, "wavelength_nm",
                      HOFFSET(WhiteboardHdf5Row, wavelength_nm), H5T_NATIVE_FLOAT);
            H5Tinsert(hit_type, "weight",
                      HOFFSET(WhiteboardHdf5Row, weight), H5T_NATIVE_FLOAT);
            H5Tinsert(hit_type, "relative_efficiency",
                      HOFFSET(WhiteboardHdf5Row, relative_efficiency), H5T_NATIVE_FLOAT);
            H5Tinsert(hit_type, "signal_weight",
                      HOFFSET(WhiteboardHdf5Row, signal_weight), H5T_NATIVE_FLOAT);
            writeCompound1D(whiteboard_group, "hits", hit_type, whiteboard_hits);
            H5Tclose(hit_type);
            H5Gclose(whiteboard_group);
        }

        H5Fclose(file);
    } catch (...) {
        H5Fclose(file);
        throw;
    }
}
#endif

std::string formatStreamSummaryLine(const TraceSummary& s,
                                    bool camera_enabled,
                                    const std::string& event_id_mode)
{
    const int array_id = arrayIdFromOutputEvent(s.event_id, event_id_mode);
    const double mean = s.weighted_signal > 0.0
        ? s.weighted_time_sum / s.weighted_signal
        : 0.0;
    const double var = s.weighted_signal > 0.0
        ? std::max(0.0, s.weighted_time2_sum / s.weighted_signal - mean * mean)
        : 0.0;
    std::ostringstream value;
    value << "array_id=" << array_id
          << " telescope=" << s.telescope_id
          << " event_id=" << s.event_id
          << " input_bunches=" << s.input_bunches
          << " input_photons=" << doubleToString(s.input_photons, 3)
          << " hit_mirror=" << s.hit_mirror
          << " hit_output=" << s.hit_output_plane;
    if (camera_enabled) {
        value << " hit_camera=" << s.hit_camera
              << " accepted=" << s.accepted_camera
              << " lost=" << s.lost_between_pixels
              << " unique_pixels=" << s.unique_pixels.size();
    }
    value << (camera_enabled ? " pe=" : " signal=")
          << doubleToString(s.weighted_signal, 3)
          << " time_mean_ns=" << doubleToString(mean, 3)
          << " time_rms_ns=" << doubleToString(std::sqrt(var), 3);
    return value.str();
}

std::string formatTelescopeEventLine(const TelescopeEventAccumulator& s,
                                     bool camera_enabled)
{
    const double mean = s.weighted_signal > 0.0
        ? s.weighted_time_sum / s.weighted_signal
        : 0.0;
    const double var = s.weighted_signal > 0.0
        ? std::max(0.0, s.weighted_time2_sum / s.weighted_signal - mean * mean)
        : 0.0;
    std::ostringstream value;
    value << "telescope=" << s.telescope_id
          << " output_events=" << s.output_events.size()
          << " input_bunches=" << s.input_bunches
          << " input_photons=" << doubleToString(s.input_photons, 3)
          << " hit_mirror=" << s.hit_mirror
          << " hit_output=" << s.hit_output_plane;
    if (camera_enabled) {
        value << " hit_camera=" << s.hit_camera
              << " accepted=" << s.accepted_camera
              << " lost=" << s.lost_between_pixels
              << " unique_pixels=" << s.unique_pixels.size();
    }
    value << (camera_enabled ? " pe=" : " signal=")
          << doubleToString(s.weighted_signal, 3)
          << " time_mean_ns=" << doubleToString(mean, 3)
          << " time_rms_ns=" << doubleToString(std::sqrt(var), 3);
    return value.str();
}

void printTraceStreamSummary(const std::map<SummaryKey, TraceSummary>& summaries,
                             bool camera_enabled,
                             const std::string& event_id_mode)
{
    printSection("Per-stream summary");
    printField("columns",
               camera_enabled
                   ? "event_id shower_event array_id telescope input_bunches input_photons hit_output hit_camera accepted unique_pixels pe time_mean_ns time_rms_ns"
                   : "event_id shower_event array_id telescope input_bunches input_photons hit_output signal time_mean_ns time_rms_ns");

    for (const auto& kv : summaries) {
        const auto& s = kv.second;
        const int shower_event = showerEventFromOutputEvent(s.event_id, event_id_mode);
        std::ostringstream value;
        value << "shower_event=" << shower_event << " "
              << formatStreamSummaryLine(s, camera_enabled, event_id_mode);
        printField("stream", value.str());
    }
}

void printEventSummary(const std::map<SummaryKey, TraceSummary>& summaries,
                       bool camera_enabled,
                       const std::string& event_id_mode,
                       const EventIOMetadata& metadata)
{
    struct EventAggregate {
        int event_id = 0;
        int shower_event = 0;
        int array_id = 0;
        std::set<int> output_events;
        std::set<int> telescopes;
        std::uint64_t input_bunches = 0;
        double input_photons = 0.0;
        std::uint64_t hit_output_plane = 0;
        std::uint64_t hit_camera = 0;
        std::uint64_t accepted_camera = 0;
        double weighted_signal = 0.0;
        double weighted_time_sum = 0.0;
        double weighted_time2_sum = 0.0;
    };

    std::map<int, EventAggregate> events;
    for (const auto& kv : summaries) {
        const auto& s = kv.second;
        const int shower_event = showerEventFromOutputEvent(s.event_id, event_id_mode);
        const int array_id = arrayIdFromOutputEvent(s.event_id, event_id_mode);
        auto& e = events[s.event_id];
        e.event_id = s.event_id;
        e.shower_event = shower_event;
        e.array_id = array_id;
        e.output_events.insert(s.event_id);
        e.telescopes.insert(s.telescope_id);
        e.input_bunches += s.input_bunches;
        e.input_photons += s.input_photons;
        e.hit_output_plane += s.hit_output_plane;
        e.hit_camera += s.hit_camera;
        e.accepted_camera += s.accepted_camera;
        e.weighted_signal += s.weighted_signal;
        e.weighted_time_sum += s.weighted_time_sum;
        e.weighted_time2_sum += s.weighted_time2_sum;
    }

    printSection("Per-event summary");
    printField("columns",
               camera_enabled
                   ? "event_id shower_event array_id energy_gev core_N_m core_E_m core_source telescopes input_bunches input_photons hit_output hit_camera accepted pe time_mean_ns time_rms_ns"
                   : "event_id shower_event array_id energy_gev core_N_m core_E_m core_source telescopes input_bunches input_photons hit_output signal time_mean_ns time_rms_ns");
    for (const auto& kv : events) {
        const auto& e = kv.second;
        const OutputEventMetadata meta = outputEventMetadata(e.event_id, event_id_mode, metadata);
        const double mean = e.weighted_signal > 0.0
            ? e.weighted_time_sum / e.weighted_signal
            : 0.0;
        const double var = e.weighted_signal > 0.0
            ? std::max(0.0, e.weighted_time2_sum / e.weighted_signal - mean * mean)
            : 0.0;
        std::ostringstream value;
        value << "event_id=" << e.event_id
              << " shower_event=" << e.shower_event
              << " array_id=" << e.array_id;
        if (meta.found) {
            value << " energy_gev=" << doubleToString(meta.energy_gev, 6)
                  << " core_N_m=" << doubleToString(meta.core_x_north_m, 3)
                  << " core_E_m=" << doubleToString(-meta.core_y_west_m, 3)
                  << " core_source="
                  << (meta.used_array_offset ? "negative_MC_TELOFF_array_offset"
                                              : "shower_header");
        } else {
            value << " energy_gev=unknown core_N_m=unknown core_E_m=unknown"
                  << " core_source=missing_metadata";
        }
        value << " output_events=" << e.output_events.size()
              << " telescopes=" << e.telescopes.size()
              << " input_bunches=" << e.input_bunches
              << " input_photons=" << doubleToString(e.input_photons, 3)
              << " hit_output=" << e.hit_output_plane;
        if (camera_enabled) {
            value << " hit_camera=" << e.hit_camera
                  << " accepted=" << e.accepted_camera;
        }
        value << (camera_enabled ? " pe=" : " signal=")
              << doubleToString(e.weighted_signal, 3)
              << " time_mean_ns=" << doubleToString(mean, 3)
              << " time_rms_ns=" << doubleToString(std::sqrt(var), 3);
        printField("event", value.str());
    }
}

void printCorsikaOpticalConfiguration(
    const std::map<std::string, std::string>& cfg,
    const TelescopeConfig& telescope_cfg,
    const TelescopeFrame& telescope_frame,
    const MirrorLayout& mirrors,
    const SyntheticPhotonConfig& source_cfg,
    const SourceRuntimeConfig& source_runtime_cfg,
    const OutputPlane& plane,
    const CorsikaTraceOutputConfig& output_cfg,
    const CameraConfig& camera_cfg,
    const CameraGeometry& camera,
    const std::unique_ptr<Cone::SquareCone>& light_collector,
    const SipmConfig& sipm_cfg,
    const ElectronicsConfig& electronics_cfg,
    const NsbConfig& nsb_cfg,
    const TriggerConfig& trigger_cfg,
    const OpticalEfficiencyConfig& efficiency_cfg,
    const ErrorConfig& error_cfg,
    const PropagationConfig& propagation_cfg)
{
    const std::string mirror_mode = lowerCopy(getString(cfg, "mirror.mode", "generated"));
    const std::string mirror_csv = getString(cfg, "mirror.csv_path", "");

    printSection("Telescope");
    printField("id", intToString(static_cast<std::uint64_t>(telescope_cfg.id)));
    printField("name", telescope_cfg.name);
    printField("position_m", vec3ToString(telescope_cfg.position_m));
    printField("pointing_az_deg", doubleToString(telescope_cfg.pointing_az_deg));
    printField("pointing_el_deg", doubleToString(telescope_cfg.pointing_el_deg));
    printField("focal_length_m", doubleToString(telescope_cfg.focal_length_m));
    printField("coordinate_system", telescope_cfg.coordinate_system);
    printField("coordinate_transform", "local telescope frame -> global frame");
    printField("frame_x_axis", vec3ToString(telescope_frame.x_axis));
    printField("frame_y_axis", vec3ToString(telescope_frame.y_axis));
    printField("frame_z_axis", vec3ToString(telescope_frame.z_axis));

    printSection("Mirror");
    printField("mode", mirror_mode);
    if (!mirror_csv.empty()) {
        printField("csv", mirror_csv);
    }
    if (mirror_mode == "elevation_series" || mirror_mode == "series") {
        printField("series_elevation_deg",
                   doubleToString(getDouble(cfg, "mirror.series_elevation_deg",
                                            telescope_cfg.pointing_el_deg)));
        printField("series_angles_deg", getString(cfg, "mirror.series_angles_deg", ""));
        printField("series_csv_pattern", getString(cfg, "mirror.series_csv_pattern", ""));
    }
    printField("facets", intToString(mirrors.size()));

    printSection("Source");
    printField("mode", "EventIO");
    printField("eventio_path", source_runtime_cfg.eventio_path);
    printField("event_id_mode", source_runtime_cfg.event_id_mode);
    printField("eventio_coordinate_frame", source_runtime_cfg.eventio_coordinate_frame);
    printField("coordinate_interpretation",
               lowerCopy(trim(source_runtime_cfg.eventio_coordinate_frame)) == "corsika_iact"
                   ? "CORSIKA IACT x/y are telescope-relative horizontal coordinates; cx/cy/cz are rotated to telescope-local optical coordinates"
                   : (lowerCopy(trim(source_runtime_cfg.eventio_coordinate_frame)) == "telescope_local" ||
                              lowerCopy(trim(source_runtime_cfg.eventio_coordinate_frame)) == "local"
                          ? "EventIO photon positions/directions are already in telescope-local optical coordinates"
                          : "EventIO photon positions/directions are global array coordinates and will be transformed to telescope-local coordinates"));
    printField("eventio_telescope_position",
               source_runtime_cfg.use_eventio_telescope_position
                   ? "use EventIO telescope table"
                   : "use telescope.position_m for every telescope");
    printField("filter_telescope_id", "off");
    printField("filter_event_id", "off");
    printField("filter_shower_event_id",
               source_runtime_cfg.filter_shower_event_id
                   ? intToString(source_runtime_cfg.selected_shower_event_id)
                   : "off");
    printField("max_shower_events",
               source_runtime_cfg.max_shower_events > 0
                   ? intToString(source_runtime_cfg.max_shower_events)
                   : "off");
    printField("default_wavelength_nm", doubleToString(source_cfg.wavelength_nm));
    printField("default_weight", doubleToString(source_cfg.photon_weight));
    printField("default_multiplicity", doubleToString(source_cfg.multiplicity));

    printSection("Output plane");
    printField("point", vec3ToString(plane.point));
    printField("normal", vec3ToString(plane.normal));
    printField("format", output_cfg.format);
    if (outputWantsHdf5(output_cfg)) {
        printField("hdf5_path", output_cfg.hdf5_path);
        printField("hdf5_storage", output_cfg.hdf5_storage);
    }
    if (camera_cfg.enabled && outputWantsCsv(output_cfg)) {
        printField("pixel_csv", output_cfg.pixel_csv);
    } else if (!camera_cfg.enabled && outputWantsCsv(output_cfg)) {
        printField("hits_csv", output_cfg.hits_csv);
    }
    if (outputWantsCsv(output_cfg)) {
        printField("summary_csv", output_cfg.summary_csv);
    }

    printSection("Camera");
    printField("enabled", camera_cfg.enabled ? "true" : "false");
    if (camera_cfg.enabled) {
        printField("mode", camera_cfg.mode);
        if (!camera_cfg.csv_path.empty()) {
            printField("csv", camera_cfg.csv_path);
        }
        printField("pixels", intToString(static_cast<std::uint64_t>(camera.size())));
        if (!camera.empty()) {
            double min_size = std::numeric_limits<double>::max();
            double max_size = 0.0;
            for (const auto& pixel : camera.pixels()) {
                min_size = std::min(min_size, pixel.size);
                max_size = std::max(max_size, pixel.size);
            }
            printField("pixel_shape", pixelShapeName(camera.pixels().front().shape));
            printField("pixel_size_range_m",
                       doubleToString(min_size) + " .. " + doubleToString(max_size));
        } else {
            printField("pixel_shape", camera_cfg.pixel_shape);
            printField("pixel_size_m", doubleToString(camera_cfg.pixel_size_m));
        }
        if (lowerCopy(trim(camera_cfg.mode)) != "csv") {
            printField("pixel_pitch_m", doubleToString(camera_cfg.pixel_pitch_m));
            printField("radius_m", doubleToString(camera_cfg.radius_m));
        }
        printField("coordinates", "output-plane local u/v");
        printField("collector", light_collector ? camera_cfg.collector
                                                : "not set -> direct pixel containment");
        if (light_collector) {
            printField("collector_material", camera_cfg.collector_material);
            if (!isDisabledText(camera_cfg.collector_reflectivity_csv)) {
                printField("collector_reflectivity_csv",
                           camera_cfg.collector_reflectivity_csv);
            }
            printField("collector_entrance_m",
                       doubleToString(cameraPixelSizeForCollector(camera_cfg, camera)));
            printField("collector_exit_m", doubleToString(camera_cfg.collector_exit_size_m));
            printField("collector_height_m", doubleToString(camera_cfg.collector_height_m));
        }
    } else {
        printField("mode", "whiteboard only");
    }

    printSection("SiPM");
    printField("size_m", doubleToString(sipm_cfg.size_m));

    printSection("Electronics");
    printField("pe_conversion", factorDescription(electronics_cfg.pe_conversion));
    printField("model", "integrated_pe_placeholder");

    printSection("NSB");
    printField("enabled", nsb_cfg.enabled ? "true" : "false");
    printField("model", nsb_cfg.model);
    printField("rate_pe_per_ns_per_pixel",
               doubleToString(nsb_cfg.rate_pe_per_ns_per_pixel));
    printField("window_ns", doubleToString(nsb_cfg.window_ns));
    printField("seed", intToString(nsb_cfg.seed));

    printSection("Trigger");
    printField("enabled", trigger_cfg.enabled ? "true" : "false");
    printField("model", "simple_multiplicity");
    printField("pixel_threshold_pe", doubleToString(trigger_cfg.pixel_threshold_pe));
    printField("camera_multiplicity", intToString(trigger_cfg.camera_multiplicity));
    printField("array_multiplicity", intToString(trigger_cfg.array_multiplicity));
    printField("coincidence_window_ns", doubleToString(trigger_cfg.coincidence_window_ns));

    printSection("Efficiency");
    printField("constant_scale", doubleToString(efficiency_cfg.constant_scale));
    printField("mirror_reflectivity",
               factorDescription(efficiency_cfg.mirror_reflectivity));
    printField("filter_transmission",
               factorDescription(efficiency_cfg.filter_transmission));
    printField("sipm_pde", factorDescription(efficiency_cfg.sipm_pde));
    printField("atmosphere", factorDescription(efficiency_cfg.atmosphere_transmission));
    printField("funnel_acceptance",
               efficiency_cfg.use_funnel_acceptance ? "cos(theta)" : "not set -> 1");

    printSection("Errors");
    printField("random_seed", intToString(error_cfg.random_seed));
    printField("structural_deformation",
               isDisabledText(error_cfg.structural_deformation_config) ? "off" : "on");
    if (!isDisabledText(error_cfg.structural_deformation_config)) {
        printField("structural_deformation_config",
                   error_cfg.structural_deformation_config);
        printField("structural_deformation_elevation_deg",
                   doubleToString(telescope_cfg.pointing_el_deg));
    }
    printField("facet_radial_pos_sigma_m",
               doubleToString(error_cfg.facet_radial_position_sigma_m));
    printField("facet_normal_sigma_deg", doubleToString(error_cfg.facet_normal_sigma_deg));
    printField("reflect_dir_sigma_deg",
               doubleToString(error_cfg.reflect_direction_sigma_deg));
    printField("radius_curvature_sigma_m",
               doubleToString(error_cfg.radius_of_curvature_sigma_m));
    printField("reflectivity_scale_sigma",
               doubleToString(error_cfg.reflectivity_scale_sigma));

    printSection("Model");
    printField("optics", "ideal reflection only");
    printField("speed_of_light_m/ns",
               doubleToString(propagation_cfg.speed_of_light_m_per_ns, 9));
    printField("not included",
               light_collector
                   ? "waveforms, SiPM saturation, crosstalk, afterpulse, dark count, hardware trigger board"
                   : "collector, waveforms, SiPM saturation, crosstalk, afterpulse, dark count, hardware trigger board");
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 2 || argc > 3) {
        std::cerr << "usage: run_corsika_trace <config.txt> [corsika_eventio_file]\n";
        return 2;
    }

    try {
#ifndef LACT_HAS_HESSIO
        throw std::runtime_error(
            "run_corsika_trace requires libhessio. Build external/hessioxxx/source "
            "and reconfigure LACT_sim.");
#else
        std::cout.setf(std::ios::unitbuf);
        std::cout << std::fixed << std::setprecision(6);
        const auto t_start = std::chrono::steady_clock::now();
        auto main_cfg = readKeyValueConfig(argv[1]);
        ComponentConfigPaths component_paths;
        auto cfg = expandConfig(main_cfg, argv[1], component_paths);

        cfg["source.mode"] = getString(cfg, "source.mode", "EventIO");
        if (!isEventIOMode(getString(cfg, "source.mode", "EventIO"))) {
            throw std::runtime_error("run_corsika_trace requires source.mode=EventIO/CORSIKA");
        }

        SyntheticPhotonConfig source_cfg = buildSourceConfig(cfg);
        SourceRuntimeConfig source_runtime_cfg = buildSourceRuntimeConfig(cfg);
        source_runtime_cfg.use_eventio = true;
        source_runtime_cfg.filter_event_id = false;
        source_runtime_cfg.filter_telescope_id = false;
        if (argc == 3) {
            source_runtime_cfg.eventio_path = argv[2];
            cfg["source.eventio_path"] = argv[2];
        }
        if (source_runtime_cfg.eventio_path.empty()) {
            source_runtime_cfg.eventio_path =
                getString(cfg, "corsika.input", getString(cfg, "eventio.input", ""));
        }
        if (source_runtime_cfg.eventio_path.empty()) {
            throw std::runtime_error("source.eventio_path or corsika.input is required");
        }
        TelescopeConfig telescope_cfg = buildTelescopeConfig(cfg);
        std::vector<MirrorFacet> facets = buildFacetsFromConfig(cfg);
        ErrorConfig error_cfg = buildErrorConfig(cfg);
        applyStructuralDeformation(facets, error_cfg, telescope_cfg);
        applyFacetErrors(facets, error_cfg);
        OutputPlane plane = buildOutputPlane(cfg);
        TelescopeFrame telescope_frame = buildTelescopeFrame(telescope_cfg);
        CameraConfig camera_cfg = buildCameraConfig(cfg);
        SipmConfig sipm_cfg = buildSipmConfig(cfg);
        ElectronicsConfig electronics_cfg = buildElectronicsConfig(cfg);
        ElectronicsResponse electronics(electronics_cfg);
        NsbConfig nsb_cfg = buildNsbConfig(cfg);
        TriggerConfig trigger_cfg = buildTriggerConfig(cfg);
        CameraGeometry camera = buildCameraGeometry(camera_cfg);
        auto light_collector = buildLightCollector(camera_cfg, camera);
        MirrorLayout mirrors = makeMirrorLayoutFromFacets(facets);
        OpticalEfficiencyConfig efficiency_cfg = buildEfficiencyConfig(cfg);
        PropagationConfig propagation_cfg = buildPropagationConfig(cfg);
        OpticalEfficiency eff(efficiency_cfg);
        OpticalTracer tracer(propagation_cfg.speed_of_light_m_per_ns,
                             error_cfg.reflect_direction_sigma_deg * DEG_TO_RAD,
                             error_cfg.random_seed);
        CorsikaTraceOutputConfig output_cfg = buildCorsikaTraceOutputConfig(cfg);
        const bool save_csv = outputWantsCsv(output_cfg);
        const bool save_hdf5 = outputWantsHdf5(output_cfg);
        if (save_hdf5) {
#ifndef LACT_HAS_HDF5
            throw std::runtime_error(
                "output.format requests HDF5, but this build was configured without HDF5. "
                "Install HDF5 and re-run CMake, or set output.format=csv.");
#endif
        }

        std::cout << "========================================\n";
        std::cout << "LACT CORSIKA/EventIO trace\n";
        std::cout << "========================================\n";
        printSection("Configuration files");
        printField("main", argv[1]);
        if (!component_paths.telescope.empty()) printField("telescope", component_paths.telescope);
        if (!component_paths.mirror.empty()) printField("mirror", component_paths.mirror);
        if (!component_paths.source.empty()) printField("source", component_paths.source);
        if (!component_paths.output.empty()) printField("output", component_paths.output);
        if (!component_paths.camera.empty()) printField("camera", component_paths.camera);
        if (!component_paths.sipm.empty()) printField("sipm", component_paths.sipm);
        if (!component_paths.electronics.empty()) printField("electronics", component_paths.electronics);
        if (!component_paths.efficiency.empty()) printField("efficiency", component_paths.efficiency);
        if (!component_paths.atmosphere.empty()) printField("atmosphere", component_paths.atmosphere);
        if (!component_paths.nsb.empty()) printField("nsb", component_paths.nsb);
        if (!component_paths.trigger.empty()) printField("trigger", component_paths.trigger);
        if (!component_paths.error.empty()) printField("error", component_paths.error);
        if (component_paths.source.empty()) printField("source", "inline EventIO settings");

        printCorsikaOpticalConfiguration(cfg,
                                         telescope_cfg,
                                         telescope_frame,
                                         mirrors,
                                         source_cfg,
                                         source_runtime_cfg,
                                         plane,
                                         output_cfg,
                                         camera_cfg,
                                         camera,
                                         light_collector,
                                         sipm_cfg,
                                         electronics_cfg,
                                         nsb_cfg,
                                         trigger_cfg,
                                         efficiency_cfg,
                                         error_cfg,
                                         propagation_cfg);
        auto eventio_cfg = buildEventIOPhotonConfig(cfg, source_cfg, source_runtime_cfg);

        printSection("Run");
        printField("status", "reading EventIO metadata");
        std::cerr << "run_corsika_trace: reading EventIO metadata from "
                  << eventio_cfg.path << "\n";
	    EventIOMetadata metadata = readEventIOMetadata(eventio_cfg);
	    printEventIOMetadataSummary(metadata);
	    printField("status", "streaming EventIO photon bunches and tracing");
        std::cerr << "run_corsika_trace: streaming photon bunches; "
                  << "no full-file photon preload is used\n";

        std::ofstream whiteboard_out;
        if (!camera_cfg.enabled && save_csv) {
            const std::filesystem::path out_path(output_cfg.hits_csv);
            if (out_path.has_parent_path()) {
                std::filesystem::create_directories(out_path.parent_path());
            }
            whiteboard_out.open(output_cfg.hits_csv);
            if (!whiteboard_out) {
                throw std::runtime_error("failed to write whiteboard CSV: " +
                                         output_cfg.hits_csv);
            }
            writeCorsikaWhiteboardHeader(whiteboard_out);
        }

        std::map<SummaryKey, TraceSummary> summaries;
        std::map<PixelKey, PixelAccumulator> pixels;
        std::vector<WhiteboardHdf5Row> whiteboard_hits;

        std::uint64_t photon_index = 0;
        int active_shower_event = -1;
        const auto t_trace_start = std::chrono::steady_clock::now();

        auto printRunningEventSummary = [&](int shower_event) {
            if (shower_event < 0) {
                return;
            }
            printSection("Event shower_event=" + intToString(shower_event));
            printField("event_status", "completed");
            printField("telescope_columns",
                       camera_cfg.enabled
                           ? "telescope output_events input_bunches input_photons hit_mirror hit_output hit_camera accepted lost unique_pixels pe time_mean_ns time_rms_ns"
                           : "telescope output_events input_bunches input_photons hit_mirror hit_output signal time_mean_ns time_rms_ns");
            std::uint64_t input_bunches = 0;
            double input_photons = 0.0;
            std::uint64_t hit_output = 0;
            std::uint64_t hit_camera = 0;
            std::uint64_t accepted = 0;
            double signal = 0.0;
            std::map<int, TelescopeEventAccumulator> telescopes;
            std::set<int> output_events;
            for (const auto& kv : summaries) {
                const auto& s = kv.second;
                if (showerEventFromOutputEvent(s.event_id,
                                               source_runtime_cfg.event_id_mode) != shower_event) {
                    continue;
                }
                input_bunches += s.input_bunches;
                input_photons += s.input_photons;
                hit_output += s.hit_output_plane;
                hit_camera += s.hit_camera;
                accepted += s.accepted_camera;
                signal += s.weighted_signal;
                output_events.insert(s.event_id);
                auto& tel = telescopes[s.telescope_id];
                tel.telescope_id = s.telescope_id;
                tel.output_events.insert(s.event_id);
                tel.input_bunches += s.input_bunches;
                tel.input_photons += s.input_photons;
                tel.hit_mirror += s.hit_mirror;
                tel.hit_output_plane += s.hit_output_plane;
                tel.hit_camera += s.hit_camera;
                tel.accepted_camera += s.accepted_camera;
                tel.lost_between_pixels += s.lost_between_pixels;
                tel.unique_pixels.insert(s.unique_pixels.begin(), s.unique_pixels.end());
                tel.weighted_signal += s.weighted_signal;
                tel.weighted_time_sum += s.weighted_time_sum;
                tel.weighted_time2_sum += s.weighted_time2_sum;
            }
            for (const auto& kv : telescopes) {
                printField("telescope",
                           formatTelescopeEventLine(kv.second, camera_cfg.enabled));
            }
            std::ostringstream value;
            value << "output_events=" << output_events.size()
                  << " telescopes=" << telescopes.size()
                  << " input_bunches=" << input_bunches
                  << " input_photons=" << doubleToString(input_photons, 3)
                  << " hit_output=" << hit_output;
            if (camera_cfg.enabled) {
                value << " hit_camera=" << hit_camera
                      << " accepted=" << accepted
                      << " pe=" << doubleToString(signal, 3);
            } else {
                value << " signal=" << doubleToString(signal, 3);
            }
            printField("event_total", value.str());
        };

        auto processBunch = [&](const PhotonBunch& raw_bunch) {
            const PhotonBunch bunch = transformEventIOBunchToTraceFrame(
                raw_bunch, telescope_cfg, metadata, source_runtime_cfg);
            ++photon_index;
            const int shower_event =
                showerEventFromOutputEvent(bunch.event_id, source_runtime_cfg.event_id_mode);
            if (active_shower_event < 0) {
                active_shower_event = shower_event;
                printSection("Event shower_event=" + intToString(shower_event));
                printField("event_status", "started");
            } else if (shower_event != active_shower_event) {
                printRunningEventSummary(active_shower_event);
                active_shower_event = shower_event;
                printSection("Event shower_event=" + intToString(shower_event));
                printField("event_status", "started");
            }
            auto& summary = summaries[{bunch.event_id, bunch.telescope_id}];
            summary.event_id = bunch.event_id;
            summary.telescope_id = bunch.telescope_id;
            summary.input_bunches += 1;
            summary.input_photons += bunch.multiplicity;

            Photon photon = bunch.photon;
            photon.normalizeDirection();
            photon.weight *= bunch.multiplicity;

            OpticalSurfaceHit hit = tracer.traceToPlane(photon, mirrors, plane, eff);
            if (hit.hit_mirror) {
                summary.hit_mirror += 1;
            }
            if (!hit.hit_surface) {
                return;
            }

            summary.hit_output_plane += 1;
            const double base_signal = hit.weight * hit.relative_efficiency;

            if (!camera_cfg.enabled) {
                if (save_hdf5) {
                    whiteboard_hits.push_back(makeWhiteboardHdf5Row(bunch, photon_index, hit));
                }
                if (save_csv) {
                    writeCorsikaWhiteboardHit(whiteboard_out, bunch, photon_index, hit);
                }
                summary.weighted_signal += base_signal;
                summary.weighted_time_sum += base_signal * hit.time_ns;
                summary.weighted_time2_sum += base_signal * hit.time_ns * hit.time_ns;
                return;
            }

            applyCameraResponse(camera, light_collector.get(), plane, sipm_cfg, electronics, hit);
            if (!hit.hit_camera) {
                summary.lost_between_pixels += 1;
                return;
            }
            summary.hit_camera += 1;
            summary.unique_pixels.insert(hit.pixel_id);
            if (hit.accepted) {
                summary.accepted_camera += 1;
            }

            const double signal = hit.weight * hit.relative_efficiency;
            summary.weighted_signal += signal;
            summary.weighted_time_sum += signal * hit.time_ns;
            summary.weighted_time2_sum += signal * hit.time_ns * hit.time_ns;

            accumulatePixelHit(pixels, bunch.event_id, bunch.telescope_id, hit);
        };

        EventIOStreamStats stream_stats = streamEventIOPhotonBunches(
            eventio_cfg,
            processBunch,
            [](const EventIOStreamProgress& progress) {
                std::ostringstream value;
                value << "photon_bunches=" << progress.photon_bunches
                      << " current_shower_event=" << progress.current_shower_event
                      << " elapsed_s=" << doubleToString(progress.elapsed_s, 3);
                printField(progress.final ? "stream_done" : "stream_progress",
                           value.str());
            });
        printRunningEventSummary(active_shower_event);
        const auto t_trace_done = std::chrono::steady_clock::now();

        if (camera_cfg.enabled && save_csv) {
            writePixelCsv(output_cfg.pixel_csv, pixels);
        }
#ifdef LACT_HAS_HDF5
        if (save_hdf5) {
            writeNativeTraceHdf5(output_cfg,
                                 argv[1],
                                 cfg,
                                 component_paths,
                                 source_runtime_cfg,
                                 telescope_cfg,
                                 metadata,
                                 camera,
                                 facets,
                                 electronics_cfg,
                                 nsb_cfg,
                                 trigger_cfg,
                                 summaries,
                                 pixels,
                                 whiteboard_hits);
        }
#endif
        if (save_csv) {
            writeSummaryCsv(output_cfg.summary_csv, summaries);
        }
        const auto t_done = std::chrono::steady_clock::now();

        printSection("Input");
        printField("config", argv[1]);
        printField("eventio_path", source_runtime_cfg.eventio_path);
        printField("event_id_mode", source_runtime_cfg.event_id_mode);
        printField("output_format", output_cfg.format);
        printField("input_bunches", intToString(stream_stats.photon_bunches));
        printField("telescopes_in_eventio", intToString(metadata.telescopes.size()));
        printField("shower_events", intToString(metadata.events.size()));
        printField("streams", intToString(summaries.size()));
        for (const auto& tel : metadata.telescopes) {
            std::ostringstream label;
            label << "tel[" << tel.telescope_id << "]";
            std::ostringstream value;
            value << "pos_m=(" << doubleToString(tel.x_m) << ", "
                  << doubleToString(tel.y_m) << ", "
                  << doubleToString(tel.z_m) << "), radius_m="
                  << doubleToString(tel.radius_m);
            printField(label.str(), value.str());
        }
        printEventSummary(summaries,
                          camera_cfg.enabled,
                          source_runtime_cfg.event_id_mode,
                          metadata);
        printSection("Detailed stream summary");
        if (save_hdf5) {
            printField("hdf5_path", output_cfg.hdf5_path);
        }
        if (save_csv) {
            printField("summary_csv", output_cfg.summary_csv);
            if (camera_cfg.enabled) {
                printField("pixel_csv", output_cfg.pixel_csv);
            } else {
                printField("hits_csv", output_cfg.hits_csv);
            }
        }
        printField("note", save_csv
                               ? "per-array/per-telescope details are written to CSV and/or HDF5, not expanded in the log"
                               : "per-array/per-telescope details are written to HDF5, not expanded in the log");
        printSection("Timing");
        printField("trace_time_s", doubleToString(elapsedSeconds(t_trace_start, t_trace_done)));
        printField("total_time_s", doubleToString(elapsedSeconds(t_start, t_done)));
        printSection("Machine-readable summary");
        std::cout << "input_bunches=" << stream_stats.photon_bunches << "\n";
        std::cout << "shower_events=" << metadata.events.size() << "\n";
        std::cout << "streams=" << summaries.size() << "\n";
        std::cout << "camera_enabled=" << (camera_cfg.enabled ? 1 : 0) << "\n";
        if (save_hdf5) {
            std::cout << "hdf5_path=" << output_cfg.hdf5_path << "\n";
        }
        if (save_csv) {
            std::cout << "summary_csv=" << output_cfg.summary_csv << "\n";
            if (camera_cfg.enabled) {
                std::cout << "pixel_csv=" << output_cfg.pixel_csv << "\n";
            } else {
                std::cout << "hits_csv=" << output_cfg.hits_csv << "\n";
            }
        }
        return 0;
#endif
    } catch (const std::exception& ex) {
        std::cerr << "run_corsika_trace error: " << ex.what() << "\n";
        return 1;
    }
}
