#include "io/CorsikaTraceHdf5Writer.hpp"

#ifdef LACT_HAS_HDF5

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <sstream>

#include <hdf5.h>

#include "app/ArrayTimingCorrection.hpp"
#include "app/TriggerResponse.hpp"
#include "io/EventIOArrayTiming.hpp"
#include "io/Hdf5WaveformWriter.hpp"

namespace lact {
namespace {

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

template <typename T>
void writePlain3D(hid_t group,
                  const std::string& name,
                  hid_t type,
                  const std::vector<T>& values,
                  hsize_t n0,
                  hsize_t n1,
                  hsize_t n2)
{
    hsize_t dims[3] = {n0, n1, n2};
    hid_t space = H5Screate_simple(3, dims, nullptr);
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

} // namespace

void writeNativeTraceHdf5(const CorsikaTraceOutputConfig& output_cfg,
                          const WaveformOutputConfig& waveform_cfg,
                          const std::string& main_config_path,
                          const std::map<std::string, std::string>& cfg,
                          const ComponentConfigPaths& component_paths,
                          const SourceRuntimeConfig& source_runtime_cfg,
                          const TelescopeConfig& telescope_cfg,
                          const EventIOMetadata& metadata,
                          const CameraGeometry& camera,
                          const std::vector<MirrorFacet>& facets,
                          const SipmConfig& sipm_cfg,
                          const ElectronicsConfig& electronics_cfg,
                          const electronics::DetectorPipelineConfig& detector_cfg,
                          const OpticalEfficiencyConfig& efficiency_cfg,
                          const NsbConfig& nsb_cfg,
                          const TriggerConfig& trigger_cfg,
                          const std::map<SummaryKey, TraceSummary>& summaries,
                          const std::map<PixelKey, PixelAccumulator>& pixels,
                          const std::map<WaveformKey, WaveformPixelAccumulator>& waveforms,
                          const std::vector<RawWaveformHit>& raw_waveform_hits,
                          const CameraElectronicsEventMap& electronics_events,
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
        writeStringAttribute(
            file, "producer_version",
            getString(cfg, "provenance.producer_version", "source-tree"));
        writeStringAttribute(
            file, "source_path",
            getString(cfg, "provenance.source_path", ""));
        writeStringAttribute(
            file, "source_sha256",
            getString(cfg, "provenance.source_sha256", ""));
        writeStringAttribute(file, "image_storage", output_cfg.hdf5_storage);
        writeStringAttribute(file, "hdf5_write_components",
                             output_cfg.hdf5_write_components ? "true" : "false");
        writeStringAttribute(file, "lact_root_write_components",
                             output_cfg.lact_root_write_components ? "true" : "false");
        writeStringAttribute(file, "save_only_triggered",
                             output_cfg.save_only_triggered ? "true" : "false");
        writeStringAttribute(file, "write_pixel_time_stats",
                             output_cfg.write_pixel_time_stats ? "true" : "false");
        writeStringAttribute(file, "waveform_enabled",
                             waveform_cfg.enabled ? "true" : "false");
        writeStringAttribute(file, "hdf5_write_waveforms",
                             output_cfg.hdf5_write_waveforms ? "true" : "false");
        writeStringAttribute(file, "hdf5_waveform_storage",
                             output_cfg.hdf5_waveform_storage);
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
            {"obstruction", component_paths.obstruction},
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
        writeStringAttribute(coordinates_group, "source_coordinate_frame",
                             source_runtime_cfg.coordinate_frame);
        writeStringAttribute(coordinates_group, "eventio_photon_frame",
                             source_runtime_cfg.coordinate_frame);
        writeStringAttribute(coordinates_group, "eventio_corsika_iact_positions",
                             "Photon bunch x/y/z are telescope-relative CORSIKA IACT coordinates before rotation to telescope-local optics");
        writeStringAttribute(coordinates_group, "eventio_teloff_core_note",
                             "MC_TELOFF is detector-array offset with respect to shower core; /events/corsika stores core = -MC_TELOFF in NWU coordinates");
        writeStringAttribute(coordinates_group, "eventio_telpos_note",
                             "Telescope positions are hessio MC_TELPOS detector coordinates; array_z_up_m may include the detector sphere/radius convention used by the producer");
        writeStringAttribute(coordinates_group, "camera_plane_coordinates",
                             "camera x/y are output-plane u/v coordinates; for the LACT focal plane u=-mirror-local x (East at north pointing) and v=+mirror-local y (sky-up)");
        H5Gclose(coordinates_group);

        hid_t sipm_group = H5Gcreate2(metadata_group, "sipm",
                                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        writeStringAttribute(sipm_group, "size_m", doubleToString(sipm_cfg.size_m));
        writeStringAttribute(sipm_group, "pde", factorDescription(efficiency_cfg.sipm_pde));
        H5Gclose(sipm_group);

        hid_t efficiency_group = H5Gcreate2(metadata_group, "efficiency",
                                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        writeStringAttribute(efficiency_group, "constant_scale",
                             doubleToString(efficiency_cfg.constant_scale));
        writeStringAttribute(efficiency_group, "mirror_reflectivity",
                             factorDescription(efficiency_cfg.mirror_reflectivity));
        writeStringAttribute(efficiency_group, "filter_transmission",
                             factorDescription(efficiency_cfg.filter_transmission));
        writeStringAttribute(efficiency_group, "atmosphere",
                             factorDescription(efficiency_cfg.atmosphere_transmission));
        writeStringAttribute(efficiency_group, "funnel_acceptance",
                             efficiency_cfg.use_funnel_acceptance ? "cos(theta)" : "not set -> 1");
        H5Gclose(efficiency_group);

        hid_t electronics_group = H5Gcreate2(metadata_group, "electronics",
                                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        writeStringAttribute(electronics_group, "enabled",
                             detector_cfg.enabled ? "true" : "false");
        writeStringAttribute(
            electronics_group, "model",
            detector_cfg.enabled ? "explicit_microcell_and_single_pe"
                                 : "disabled");
        writeStringAttribute(electronics_group, "microcell_saturation_enabled",
                             detector_cfg.microcell.saturation_enabled
                                 ? "true"
                                 : "false");
        writeStringAttribute(electronics_group, "single_pe_enabled",
                             detector_cfg.single_pe.enabled ? "true" : "false");
        writeStringAttribute(electronics_group, "single_pe_model",
                             detector_cfg.single_pe.model);
        writeStringAttribute(electronics_group, "single_pe_unit",
                             detector_cfg.single_pe.unit);
        writeStringAttribute(electronics_group, "template_time_reference",
                             detector_cfg.single_pe.template_time_reference);
        writeStringAttribute(
            electronics_group, "charge_fluctuation_enabled",
            detector_cfg.single_pe.charge_fluctuation.enabled ? "true" : "false");
        writeStringAttribute(electronics_group, "charge_fluctuation_model",
                             detector_cfg.single_pe.charge_fluctuation.model);
        writeStringAttribute(electronics_group, "time_jitter_enabled",
                             detector_cfg.single_pe.time_jitter.enabled
                                 ? "true"
                                 : "false");
        writeStringAttribute(electronics_group, "response",
                             "SiPM PDE is applied before this detector pipeline");
        H5Gclose(electronics_group);

        hid_t waveform_group_meta = H5Gcreate2(metadata_group, "waveform",
                                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        writeStringAttribute(waveform_group_meta, "enabled",
                             waveform_cfg.enabled ? "true" : "false");
        writeStringAttribute(waveform_group_meta, "source", waveform_cfg.source);
        writeStringAttribute(waveform_group_meta, "time_reference",
                             waveform_cfg.time_reference);
        writeStringAttribute(waveform_group_meta, "time_bin_width_ns",
                             doubleToString(waveform_cfg.time_bin_width_ns));
        writeStringAttribute(waveform_group_meta, "time_window_start_ns",
                             doubleToString(waveform_cfg.time_window_start_ns));
        writeStringAttribute(waveform_group_meta, "time_window_end_ns",
                             doubleToString(waveform_cfg.time_window_end_ns));
        writeStringAttribute(
            waveform_group_meta, "note",
            waveform_cfg.source == "electronics"
                ? "measured single-p.e. pulses are superposed after microcell saturation and sampled at the configured interval"
                : "time-binned camera p.e. output; when NSB is enabled it is sampled per time bin");
        H5Gclose(waveform_group_meta);

        hid_t nsb_group = H5Gcreate2(metadata_group, "nsb",
                                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        writeStringAttribute(nsb_group, "enabled", nsb_cfg.enabled ? "true" : "false");
        writeStringAttribute(nsb_group, "model", nsb_cfg.model);
        writeStringAttribute(nsb_group, "rate_pe_per_ns_per_pixel",
                             doubleToString(nsb_cfg.rate_pe_per_ns_per_pixel));
        writeStringAttribute(nsb_group, "window_ns", doubleToString(nsb_cfg.window_ns));
        writeStringAttribute(nsb_group, "seed", intToString(nsb_cfg.seed));
        writeStringAttribute(nsb_group, "spectrum_csv", nsb_cfg.spectrum_csv);
        writeStringAttribute(nsb_group, "spectrum_unit", nsb_cfg.spectrum_unit);
        writeStringAttribute(nsb_group, "effective_area_m2",
                             doubleToString(nsb_cfg.effective_area_m2));
        writeStringAttribute(nsb_group, "collector_mean_transmission",
                             doubleToString(
                                 nsb_cfg.collector_mean_transmission));
        writeStringAttribute(nsb_group, "microcell_geometric_acceptance",
                             doubleToString(
                                 nsb_cfg.microcell_geometric_acceptance));
        writeStringAttribute(nsb_group, "pixel_solid_angle_sr",
                             doubleToString(nsb_cfg.pixel_solid_angle_sr));
        writeStringAttribute(nsb_group, "computed_from_spectrum",
                             nsb_cfg.computed_from_spectrum ? "true" : "false");
        writeStringAttribute(nsb_group, "spectral_integral_pe_s_sr_m2",
                             doubleToString(nsb_cfg.spectral_integral_pe_s_sr_m2));
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
        writeStringAttribute(trigger_group_meta, "camera_coincidence_window_ns",
                             doubleToString(trigger_cfg.camera_coincidence_window_ns));
        writeStringAttribute(trigger_group_meta, "array_coincidence_window_ns",
                             doubleToString(trigger_cfg.array_coincidence_window_ns));
        writeStringAttribute(trigger_group_meta, "array_time_correction",
                             trigger_cfg.array_time_correction);
        writeStringAttribute(trigger_group_meta, "array_wavefront_speed_m_per_ns",
                             doubleToString(
                                 trigger_cfg.array_time_correction == "plane_wave"
                                     ? resolveEventIOArrayWavefrontSpeedMPerNs(
                                           trigger_cfg, metadata)
                                     : trigger_cfg.array_wavefront_speed_m_per_ns,
                                 15));
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
            float primary_cherenkov_pe;
            float primary_nsb_pe;
            float primary_dark_pe;
            float fired_cherenkov_pe;
            float fired_nsb_pe;
            float fired_dark_pe;
            float gap_lost_pe;
            float saturation_lost_pe;
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
            double time_mean_ns;
            double time_rms_ns;
            double time_first_ns;
            double reference_time_ns;
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
            const auto canonical = electronics_events.find(key);
            if (detector_cfg.enabled && canonical == electronics_events.end()) {
                throw std::runtime_error(
                    "HDF5 writer has no canonical electronics result for event " +
                    std::to_string(event_id) + " telescope " +
                    std::to_string(telescope_id));
            }
            if (detector_cfg.enabled) {
                const auto& result = canonical->second.detector;
                if (result.pixels.size() != camera_rows.size()) {
                    throw std::runtime_error(
                        "canonical electronics pixel axis does not match HDF5 camera");
                }
                for (std::size_t col = 0; col < result.pixels.size(); ++col) {
                    const int pixel_id = camera_rows[col].pixel_id;
                    const auto optical = pixels.find(
                        PixelKey{event_id, telescope_id, pixel_id});
                    const PixelAccumulator* p =
                        optical == pixels.end() ? nullptr : &optical->second;
                    const double mean = p && p->signal > 0.0
                        ? p->time_sum / p->signal : 0.0;
                    const double var = p && p->signal > 0.0
                        ? std::max(0.0, p->time2_sum / p->signal - mean * mean)
                        : 0.0;
                    const auto& detector_pixel = result.pixels[col];
                    const double fired_pe =
                        detector_pixel.fired_cherenkov_pe +
                        detector_pixel.fired_nsb_pe +
                        detector_pixel.fired_dark_pe;
                    const double primary_pe =
                        detector_pixel.primary_cherenkov_pe +
                        detector_pixel.primary_nsb_pe +
                        detector_pixel.primary_dark_pe;
                    if (fired_pe <= 0.0 && primary_pe <= 0.0) {
                        continue;
                    }
                    sparse_rows.push_back(SparsePixelRow{
                        static_cast<std::int32_t>(pixel_id),
                        static_cast<std::int32_t>(p ? p->photon_count : 0),
                        static_cast<float>(fired_pe),
                        static_cast<float>(fired_pe),
                        static_cast<float>(detector_pixel.primary_cherenkov_pe),
                        static_cast<float>(detector_pixel.primary_nsb_pe),
                        static_cast<float>(detector_pixel.primary_dark_pe),
                        static_cast<float>(detector_pixel.fired_cherenkov_pe),
                        static_cast<float>(detector_pixel.fired_nsb_pe),
                        static_cast<float>(detector_pixel.fired_dark_pe),
                        static_cast<float>(detector_pixel.gap_lost_pe),
                        static_cast<float>(detector_pixel.saturation_lost_pe),
                        static_cast<float>(mean),
                        static_cast<float>(std::sqrt(var)),
                    });
                    total_signal += fired_pe;
                    total_pe += fired_pe;
                    total_photons += static_cast<double>(p ? p->photon_count : 0);
                }
            } else {
                const PixelKey pixel_begin{
                    event_id, telescope_id, std::numeric_limits<int>::min()};
                const PixelKey pixel_end{
                    event_id, telescope_id, std::numeric_limits<int>::max()};
                for (auto it = pixels.lower_bound(pixel_begin);
                     it != pixels.end() && it->first <= pixel_end;
                     ++it) {
                    const auto& p = it->second;
                    const double mean = p.signal > 0.0
                        ? p.time_sum / p.signal : 0.0;
                    const double var = p.signal > 0.0
                        ? std::max(0.0, p.time2_sum / p.signal - mean * mean)
                        : 0.0;
                    sparse_rows.push_back(SparsePixelRow{
                        static_cast<std::int32_t>(p.pixel_id),
                        static_cast<std::int32_t>(p.photon_count),
                        static_cast<float>(p.pe),
                        static_cast<float>(p.signal),
                        static_cast<float>(p.cherenkov_pe),
                        static_cast<float>(p.nsb_pe),
                        static_cast<float>(p.dark_pe),
                        static_cast<float>(p.cherenkov_pe),
                        static_cast<float>(p.nsb_pe),
                        static_cast<float>(p.dark_pe), 0.0f, 0.0f,
                        static_cast<float>(mean),
                        static_cast<float>(std::sqrt(var)),
                    });
                    total_signal += p.signal;
                    total_pe += p.pe;
                    total_photons += static_cast<double>(p.photon_count);
                }
            }

            double mean = 0.0;
            double rms = 0.0;
            double first = 0.0;
            auto summary_it = summaries.find(key);
            if (summary_it != summaries.end()) {
                const auto& s = summary_it->second;
                if (!detector_cfg.enabled) {
                    total_signal = s.weighted_signal;
                    total_pe = s.weighted_signal;
                }
                total_photons = static_cast<double>(s.hit_camera);
                if (std::isfinite(s.first_cherenkov_time_ns)) {
                    first = s.first_cherenkov_time_ns;
                }
                if (s.weighted_signal > 0.0) {
                    const double m = s.weighted_time_sum / s.weighted_signal;
                    const double v = std::max(0.0, s.weighted_time2_sum / s.weighted_signal - m * m);
                    mean = m;
                    rms = std::sqrt(v);
                }
            }
            const auto count = static_cast<std::int32_t>(
                static_cast<std::int64_t>(sparse_rows.size()) - start);
            double reference_time_ns = 0.0;
            if (detector_cfg.enabled) {
                reference_time_ns = electronics_events.at(key).reference_time_ns;
            } else if (waveform_cfg.time_reference == "image_first") {
                reference_time_ns = first;
            } else if (waveform_cfg.time_reference == "image_mean") {
                reference_time_ns = mean;
            }
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
                first,
                reference_time_ns,
            });
        }

        std::vector<std::int32_t> pixel_id_axis;
        std::map<std::int32_t, std::size_t> pixel_to_col;
        std::vector<float> dense_signal;
        std::vector<float> dense_pe;
        std::vector<float> dense_cherenkov_pe;
        std::vector<float> dense_nsb_pe;
        std::vector<float> dense_primary_dark_pe;
        std::vector<float> dense_fired_cherenkov_pe;
        std::vector<float> dense_fired_nsb_pe;
        std::vector<float> dense_fired_dark_pe;
        std::vector<float> dense_gap_lost_pe;
        std::vector<float> dense_saturation_lost_pe;
        std::vector<float> dense_time_mean_ns;
        std::vector<float> dense_time_rms_ns;
        std::vector<std::int32_t> dense_photon_count;
        const bool have_dense_images = write_dense && !camera_rows.empty();
        const bool have_camera_axis = !camera_rows.empty();
        if (have_camera_axis) {
            pixel_id_axis.reserve(camera_rows.size());
            for (std::size_t i = 0; i < camera_rows.size(); ++i) {
                pixel_id_axis.push_back(camera_rows[i].pixel_id);
                pixel_to_col[camera_rows[i].pixel_id] = i;
            }
        }

        if (have_dense_images) {
            const std::size_t n_images = image_rows.size();
            const std::size_t n_pixels = camera_rows.size();
            dense_signal.assign(n_images * n_pixels, 0.0f);
            dense_pe.assign(n_images * n_pixels, 0.0f);
            dense_photon_count.assign(n_images * n_pixels, 0);
            if (detector_cfg.enabled || output_cfg.hdf5_write_components) {
                dense_cherenkov_pe.assign(n_images * n_pixels, 0.0f);
                dense_nsb_pe.assign(n_images * n_pixels, 0.0f);
                dense_primary_dark_pe.assign(n_images * n_pixels, 0.0f);
                dense_fired_cherenkov_pe.assign(n_images * n_pixels, 0.0f);
                dense_fired_nsb_pe.assign(n_images * n_pixels, 0.0f);
                dense_fired_dark_pe.assign(n_images * n_pixels, 0.0f);
                dense_gap_lost_pe.assign(n_images * n_pixels, 0.0f);
                dense_saturation_lost_pe.assign(n_images * n_pixels, 0.0f);
            }
            if (output_cfg.write_pixel_time_stats) {
                dense_time_mean_ns.assign(n_images * n_pixels, 0.0f);
                dense_time_rms_ns.assign(n_images * n_pixels, 0.0f);
            }
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
                    if (!dense_cherenkov_pe.empty()) {
                        dense_cherenkov_pe[index] =
                            pixel.primary_cherenkov_pe;
                        dense_nsb_pe[index] = pixel.primary_nsb_pe;
                        dense_primary_dark_pe[index] = pixel.primary_dark_pe;
                        dense_fired_cherenkov_pe[index] =
                            pixel.fired_cherenkov_pe;
                        dense_fired_nsb_pe[index] = pixel.fired_nsb_pe;
                        dense_fired_dark_pe[index] = pixel.fired_dark_pe;
                        dense_gap_lost_pe[index] = pixel.gap_lost_pe;
                        dense_saturation_lost_pe[index] =
                            pixel.saturation_lost_pe;
                    }
                    if (output_cfg.write_pixel_time_stats) {
                        dense_time_mean_ns[index] = pixel.time_mean_ns;
                        dense_time_rms_ns[index] = pixel.time_rms_ns;
                    }
                }
            }

            if (!detector_cfg.enabled && output_cfg.hdf5_write_components) {
                dense_cherenkov_pe = dense_pe;
                dense_nsb_pe.assign(n_images * n_pixels, 0.0f);
                dense_primary_dark_pe.assign(n_images * n_pixels, 0.0f);
                dense_fired_cherenkov_pe = dense_pe;
                dense_fired_nsb_pe.assign(n_images * n_pixels, 0.0f);
                dense_fired_dark_pe.assign(n_images * n_pixels, 0.0f);
                dense_gap_lost_pe.assign(n_images * n_pixels, 0.0f);
                dense_saturation_lost_pe.assign(n_images * n_pixels, 0.0f);
            }
            if (!detector_cfg.enabled && nsb_cfg.enabled &&
                nsb_cfg.rate_pe_per_ns_per_pixel > 0.0 &&
                nsb_cfg.window_ns > 0.0) {
                if (waveform_cfg.enabled && waveform_cfg.source == "pe") {
                    const std::size_t n_bins = waveformBinCount(waveform_cfg);
                    for (const auto& image : image_rows) {
                        const std::size_t row =
                            static_cast<std::size_t>(image.image_index);
                        generateTimeBinnedNsbPe(
                            nsb_cfg,
                            waveform_cfg,
                            static_cast<int>(image.event_id),
                            static_cast<int>(image.telescope_id),
                            n_pixels,
                            n_bins,
                            [&](std::size_t col, std::size_t bin, float nsb_pe) {
                                const std::size_t index = row * n_pixels + col;
                                if (nsbTimeInImageWindow(
                                        nsb_cfg,
                                        waveform_cfg.time_window_start_ns +
                                            (static_cast<double>(bin) + 0.5) *
                                                waveform_cfg.time_bin_width_ns)) {
                                    if (!dense_nsb_pe.empty()) {
                                        dense_nsb_pe[index] += nsb_pe;
                                        dense_fired_nsb_pe[index] += nsb_pe;
                                    }
                                    dense_pe[index] += nsb_pe;
                                    dense_signal[index] += nsb_pe;
                                }
                            });
                    }
                } else {
                    for (const auto& image : image_rows) {
                        const std::size_t row =
                            static_cast<std::size_t>(image.image_index);
                        generateIntegratedNsbPe(
                            nsb_cfg,
                            static_cast<int>(image.event_id),
                            static_cast<int>(image.telescope_id),
                            n_pixels,
                            nsb_cfg.window_ns,
                            [&](std::size_t col, float nsb_pe) {
                                const std::size_t i = row * n_pixels + col;
                                if (!dense_nsb_pe.empty()) {
                                    dense_nsb_pe[i] = nsb_pe;
                                    dense_fired_nsb_pe[i] = nsb_pe;
                                }
                                dense_pe[i] += nsb_pe;
                                dense_signal[i] += nsb_pe;
                            });
                    }
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
            double trigger_time_ns;
            double trigger_first_time_ns;
            double trigger_max_multiplicity_time_ns;
            double geometric_delay_ns = std::numeric_limits<double>::quiet_NaN();
            double coincidence_time_ns = std::numeric_limits<double>::quiet_NaN();
        };
        struct ArrayTriggerRow {
            std::int64_t event_id;
            std::int8_t array_triggered;
            std::int32_t n_triggered_telescopes;
        };
        std::vector<TelescopeTriggerRow> telescope_trigger_rows;
        std::map<SummaryKey, std::size_t> telescope_trigger_row_index;
        std::map<int, std::vector<TelescopeTriggerTime>>
            telescope_trigger_times_by_event;
        telescope_trigger_rows.reserve(image_rows.size());
        for (const auto& image : image_rows) {
            const SummaryKey image_key{
                static_cast<int>(image.event_id),
                static_cast<int>(image.telescope_id)};
            if (detector_cfg.enabled) {
                const auto canonical = electronics_events.find(image_key);
                if (canonical == electronics_events.end()) {
                    throw std::runtime_error(
                        "HDF5 trigger writer has no canonical electronics event");
                }
                const auto& event = canonical->second;
                const auto& camera_trigger = event.detector.camera_trigger;
                if (camera_trigger.triggered) {
                    telescope_trigger_times_by_event[
                        static_cast<int>(image.event_id)].push_back(
                            TelescopeTriggerTime{
                                static_cast<int>(image.telescope_id),
                                event.trigger_time_ns,
                                event.coincidence_time_ns,
                                std::isfinite(event.geometric_delay_ns)
                                    ? event.geometric_delay_ns : 0.0,
                            });
                }
                telescope_trigger_row_index[image_key] =
                    telescope_trigger_rows.size();
                telescope_trigger_rows.push_back(TelescopeTriggerRow{
                    image.event_id,
                    image.telescope_id,
                    static_cast<std::int8_t>(camera_trigger.triggered ? 1 : 0),
                    static_cast<std::int32_t>(
                        camera_trigger.max_pixels_above_threshold),
                    image.total_pe,
                    event.trigger_time_ns,
                    event.trigger_time_ns,
                    event.trigger_time_ns,
                    event.geometric_delay_ns,
                    event.coincidence_time_ns,
                });
                continue;
            }
            double total_pe = image.total_pe;
            std::size_t trigger_pixels = 0;
            std::size_t trigger_bins = 1;
            double trigger_bin_width_ns = 1.0;
            double first_trigger_bin_center_ns =
                std::isfinite(image.time_mean_ns) ? image.time_mean_ns : 0.0;
            std::function<double(std::size_t, std::size_t)> pe_at;
            std::unordered_map<std::size_t, double> trigger_waveform_pe;
            std::vector<double> sparse_trigger_pe;
            const bool use_time_trigger =
                waveform_cfg.enabled && waveform_cfg.source == "pe" &&
                !camera_rows.empty();
            if (use_time_trigger) {
                trigger_pixels = camera_rows.size();
                trigger_bins = waveformBinCount(waveform_cfg);
                trigger_bin_width_ns = waveform_cfg.time_bin_width_ns;
                double reference_time_ns = 0.0;
                if (waveform_cfg.time_reference == "image_first" &&
                    std::isfinite(image.time_first_ns)) {
                    reference_time_ns = image.time_first_ns;
                } else if (waveform_cfg.time_reference == "image_mean" &&
                           std::isfinite(image.time_mean_ns)) {
                    reference_time_ns = image.time_mean_ns;
                }
                first_trigger_bin_center_ns =
                    reference_time_ns + waveform_cfg.time_window_start_ns +
                    0.5 * waveform_cfg.time_bin_width_ns;
                auto add_trigger_pe = [&](int pixel_id, int bin, double pe) {
                    const auto col_it = pixel_to_col.find(pixel_id);
                    if (col_it == pixel_to_col.end() || bin < 0 ||
                        static_cast<std::size_t>(bin) >= trigger_bins || pe == 0.0) {
                        return;
                    }
                    trigger_waveform_pe[
                        col_it->second * trigger_bins + static_cast<std::size_t>(bin)] += pe;
                };
                if (waveformUsesImageReference(waveform_cfg)) {
                    for (const auto& hit : raw_waveform_hits) {
                        if (hit.event_id != image.event_id ||
                            hit.telescope_id != image.telescope_id) {
                            continue;
                        }
                        add_trigger_pe(
                            hit.pixel_id,
                            waveformBinForTime(
                                waveform_cfg, hit.time_ns - reference_time_ns),
                            hit.pe);
                    }
                } else {
                    const WaveformKey begin_key{
                        static_cast<int>(image.event_id),
                        static_cast<int>(image.telescope_id),
                        std::numeric_limits<int>::min(),
                        std::numeric_limits<int>::min()};
                    const WaveformKey end_key{
                        static_cast<int>(image.event_id),
                        static_cast<int>(image.telescope_id),
                        std::numeric_limits<int>::max(),
                        std::numeric_limits<int>::max()};
                    for (auto it = waveforms.lower_bound(begin_key);
                         it != waveforms.end() && it->first <= end_key;
                         ++it) {
                        add_trigger_pe(
                            it->second.pixel_id, it->second.time_bin, it->second.pe);
                    }
                }

                // Use the same deterministic NSB realization as the dense
                // image and serialized waveform.  A separately seeded
                // per-cell sampler would have the same distribution but
                // would make the trigger impossible to reproduce from the
                // saved waveform.
                generateTimeBinnedNsbPe(
                    nsb_cfg,
                    waveform_cfg,
                    static_cast<int>(image.event_id),
                    static_cast<int>(image.telescope_id),
                    trigger_pixels,
                    trigger_bins,
                    [&](std::size_t col, std::size_t bin, float pe) {
                        trigger_waveform_pe[col * trigger_bins + bin] += pe;
                    });
                pe_at = [&](std::size_t col, std::size_t bin) {
                    const auto found = trigger_waveform_pe.find(
                        col * trigger_bins + bin);
                    return found == trigger_waveform_pe.end()
                        ? 0.0
                        : found->second;
                };
            } else if (have_dense_images) {
                const std::size_t row = static_cast<std::size_t>(image.image_index);
                const std::size_t n_pixels = camera_rows.size();
                total_pe = 0.0;
                for (std::size_t col = 0; col < n_pixels; ++col) {
                    total_pe += dense_pe[row * n_pixels + col];
                }
                trigger_pixels = n_pixels;
                pe_at = [&, row, n_pixels](std::size_t col, std::size_t) {
                    return static_cast<double>(dense_pe[row * n_pixels + col]);
                };
            } else {
                const std::int64_t begin = image.start;
                const std::int64_t end = image.start + image.count;
                trigger_pixels = camera_rows.size();
                sparse_trigger_pe.assign(trigger_pixels, 0.0);
                for (std::int64_t i = begin; i < end; ++i) {
                    const auto& pixel = sparse_rows[static_cast<std::size_t>(i)];
                    const auto col_it = pixel_to_col.find(pixel.pixel_id);
                    if (col_it != pixel_to_col.end()) {
                        sparse_trigger_pe[col_it->second] += pixel.pe;
                    }
                }
                pe_at = [&](std::size_t col, std::size_t) {
                    return sparse_trigger_pe[col];
                };
            }
            const auto camera_trigger = evaluateBinnedPeTrigger(
                trigger_pixels,
                trigger_bins,
                trigger_bin_width_ns,
                first_trigger_bin_center_ns,
                trigger_cfg, pe_at);
            if (camera_trigger.triggered) {
                telescope_trigger_times_by_event[static_cast<int>(image.event_id)]
                    .push_back(TelescopeTriggerTime{
                        static_cast<int>(image.telescope_id),
                        camera_trigger.trigger_time_ns});
            }
            telescope_trigger_row_index[{
                static_cast<int>(image.event_id),
                static_cast<int>(image.telescope_id)}] =
                telescope_trigger_rows.size();
            telescope_trigger_rows.push_back(TelescopeTriggerRow{
                image.event_id,
                image.telescope_id,
                static_cast<std::int8_t>(camera_trigger.triggered ? 1 : 0),
                static_cast<std::int32_t>(
                    camera_trigger.n_pixels_above_threshold),
                total_pe,
                camera_trigger.trigger_time_ns,
                camera_trigger.first_trigger_time_ns,
                camera_trigger.max_multiplicity_time_ns,
            });
        }
        std::set<int> trigger_event_ids;
        for (const auto& key : image_keys) {
            trigger_event_ids.insert(key.first);
        }
        std::vector<ArrayTriggerRow> array_trigger_rows;
        std::map<int, ArrayTriggerDecision> array_trigger_decisions;
        array_trigger_rows.reserve(trigger_event_ids.size());
        for (const int event_id : trigger_event_ids) {
            applyEventIOArrayTimingCorrection(
                telescope_trigger_times_by_event[event_id], event_id,
                source_runtime_cfg.event_id_mode, trigger_cfg,
                telescope_cfg, metadata);
            for (const auto& trigger_time :
                 telescope_trigger_times_by_event[event_id]) {
                const auto row_index = telescope_trigger_row_index.find({
                    event_id, trigger_time.telescope_id});
                if (row_index == telescope_trigger_row_index.end()) {
                    throw std::runtime_error(
                        "array timing result has no HDF5 telescope trigger row");
                }
                auto& row = telescope_trigger_rows[row_index->second];
                row.geometric_delay_ns = trigger_time.geometric_delay_ns;
                row.coincidence_time_ns =
                    std::isfinite(trigger_time.coincidence_time_ns)
                        ? trigger_time.coincidence_time_ns
                        : trigger_time.trigger_time_ns;
            }
            const auto decision = evaluateArrayTrigger(
                telescope_trigger_times_by_event[event_id], trigger_cfg);
            array_trigger_decisions[event_id] = decision;
            const int n_triggered = static_cast<int>(
                telescope_trigger_times_by_event[event_id].size());
            array_trigger_rows.push_back(ArrayTriggerRow{
                static_cast<std::int64_t>(event_id),
                static_cast<std::int8_t>(decision.triggered ? 1 : 0),
                static_cast<std::int32_t>(n_triggered),
            });
        }

        std::set<SummaryKey> selected_image_keys;
        bool filter_images = false;
        if (output_cfg.save_only_triggered && trigger_cfg.enabled) {
            filter_images = true;
            for (const auto& row : telescope_trigger_rows) {
                const auto& array_decision =
                    array_trigger_decisions[static_cast<int>(row.event_id)];
                const bool telescope_is_coincident = std::binary_search(
                    array_decision.coincident_telescope_ids.begin(),
                    array_decision.coincident_telescope_ids.end(),
                    static_cast<int>(row.telescope_id));
                if (row.triggered && array_decision.triggered &&
                    telescope_is_coincident) {
                    selected_image_keys.insert({
                        static_cast<int>(row.event_id),
                        static_cast<int>(row.telescope_id),
                    });
                }
            }
        } else if (!output_cfg.save_only_triggered) {
            filter_images = true;
            for (const auto& image : image_rows) {
                if (image.total_pe > 0.0) {
                    selected_image_keys.insert({
                        static_cast<int>(image.event_id),
                        static_cast<int>(image.telescope_id),
                    });
                }
            }
        }

        if (filter_images) {

            std::vector<ImageIndexRow> filtered_image_rows;
            std::vector<SparsePixelRow> filtered_sparse_rows;
            std::vector<float> filtered_dense_signal;
            std::vector<float> filtered_dense_pe;
            std::vector<float> filtered_dense_cherenkov_pe;
            std::vector<float> filtered_dense_nsb_pe;
            std::vector<float> filtered_dense_primary_dark_pe;
            std::vector<float> filtered_dense_fired_cherenkov_pe;
            std::vector<float> filtered_dense_fired_nsb_pe;
            std::vector<float> filtered_dense_fired_dark_pe;
            std::vector<float> filtered_dense_gap_lost_pe;
            std::vector<float> filtered_dense_saturation_lost_pe;
            std::vector<float> filtered_dense_time_mean_ns;
            std::vector<float> filtered_dense_time_rms_ns;
            std::vector<std::int32_t> filtered_dense_photon_count;

            const std::size_t n_pixels = camera_rows.size();
            if (have_dense_images) {
                filtered_dense_signal.reserve(selected_image_keys.size() * n_pixels);
                filtered_dense_pe.reserve(selected_image_keys.size() * n_pixels);
                filtered_dense_photon_count.reserve(selected_image_keys.size() * n_pixels);
                if (!dense_cherenkov_pe.empty()) {
                    filtered_dense_cherenkov_pe.reserve(selected_image_keys.size() * n_pixels);
                    filtered_dense_nsb_pe.reserve(selected_image_keys.size() * n_pixels);
                    filtered_dense_primary_dark_pe.reserve(selected_image_keys.size() * n_pixels);
                    filtered_dense_fired_cherenkov_pe.reserve(selected_image_keys.size() * n_pixels);
                    filtered_dense_fired_nsb_pe.reserve(selected_image_keys.size() * n_pixels);
                    filtered_dense_fired_dark_pe.reserve(selected_image_keys.size() * n_pixels);
                    filtered_dense_gap_lost_pe.reserve(selected_image_keys.size() * n_pixels);
                    filtered_dense_saturation_lost_pe.reserve(selected_image_keys.size() * n_pixels);
                }
                if (output_cfg.write_pixel_time_stats) {
                    filtered_dense_time_mean_ns.reserve(selected_image_keys.size() * n_pixels);
                    filtered_dense_time_rms_ns.reserve(selected_image_keys.size() * n_pixels);
                }
            }

            for (const auto& image : image_rows) {
                const SummaryKey key{
                    static_cast<int>(image.event_id),
                    static_cast<int>(image.telescope_id),
                };
                if (selected_image_keys.find(key) == selected_image_keys.end()) {
                    continue;
                }

                ImageIndexRow filtered = image;
                filtered.image_index = static_cast<std::int32_t>(filtered_image_rows.size());
                filtered.start = static_cast<std::int64_t>(filtered_sparse_rows.size());

                const std::int64_t begin = image.start;
                const std::int64_t end = image.start + image.count;
                for (std::int64_t i = begin; i < end; ++i) {
                    filtered_sparse_rows.push_back(sparse_rows[static_cast<std::size_t>(i)]);
                }
                filtered.count = static_cast<std::int32_t>(
                    static_cast<std::int64_t>(filtered_sparse_rows.size()) - filtered.start);

                if (have_dense_images) {
                    const std::size_t old_row = static_cast<std::size_t>(image.image_index);
                    const std::size_t old_begin = old_row * n_pixels;
                    const std::size_t old_end = old_begin + n_pixels;
                    filtered_dense_signal.insert(filtered_dense_signal.end(),
                                                 dense_signal.begin() + old_begin,
                                                 dense_signal.begin() + old_end);
                    filtered_dense_pe.insert(filtered_dense_pe.end(),
                                             dense_pe.begin() + old_begin,
                                             dense_pe.begin() + old_end);
                    filtered_dense_photon_count.insert(filtered_dense_photon_count.end(),
                                                       dense_photon_count.begin() + old_begin,
                                                       dense_photon_count.begin() + old_end);
                    if (!dense_cherenkov_pe.empty()) {
                        filtered_dense_cherenkov_pe.insert(filtered_dense_cherenkov_pe.end(),
                                                           dense_cherenkov_pe.begin() + old_begin,
                                                           dense_cherenkov_pe.begin() + old_end);
                        filtered_dense_nsb_pe.insert(filtered_dense_nsb_pe.end(),
                                                     dense_nsb_pe.begin() + old_begin,
                                                     dense_nsb_pe.begin() + old_end);
                        filtered_dense_primary_dark_pe.insert(
                            filtered_dense_primary_dark_pe.end(),
                            dense_primary_dark_pe.begin() + old_begin,
                            dense_primary_dark_pe.begin() + old_end);
                        filtered_dense_fired_cherenkov_pe.insert(
                            filtered_dense_fired_cherenkov_pe.end(),
                            dense_fired_cherenkov_pe.begin() + old_begin,
                            dense_fired_cherenkov_pe.begin() + old_end);
                        filtered_dense_fired_nsb_pe.insert(
                            filtered_dense_fired_nsb_pe.end(),
                            dense_fired_nsb_pe.begin() + old_begin,
                            dense_fired_nsb_pe.begin() + old_end);
                        filtered_dense_fired_dark_pe.insert(
                            filtered_dense_fired_dark_pe.end(),
                            dense_fired_dark_pe.begin() + old_begin,
                            dense_fired_dark_pe.begin() + old_end);
                        filtered_dense_gap_lost_pe.insert(
                            filtered_dense_gap_lost_pe.end(),
                            dense_gap_lost_pe.begin() + old_begin,
                            dense_gap_lost_pe.begin() + old_end);
                        filtered_dense_saturation_lost_pe.insert(
                            filtered_dense_saturation_lost_pe.end(),
                            dense_saturation_lost_pe.begin() + old_begin,
                            dense_saturation_lost_pe.begin() + old_end);
                    }
                    if (output_cfg.write_pixel_time_stats) {
                        filtered_dense_time_mean_ns.insert(filtered_dense_time_mean_ns.end(),
                                                           dense_time_mean_ns.begin() + old_begin,
                                                           dense_time_mean_ns.begin() + old_end);
                        filtered_dense_time_rms_ns.insert(filtered_dense_time_rms_ns.end(),
                                                          dense_time_rms_ns.begin() + old_begin,
                                                          dense_time_rms_ns.begin() + old_end);
                    }
                }

                filtered_image_rows.push_back(filtered);
            }

            image_rows.swap(filtered_image_rows);
            sparse_rows.swap(filtered_sparse_rows);
            if (have_dense_images) {
                dense_signal.swap(filtered_dense_signal);
                dense_pe.swap(filtered_dense_pe);
                dense_photon_count.swap(filtered_dense_photon_count);
                if (!dense_cherenkov_pe.empty()) {
                    dense_cherenkov_pe.swap(filtered_dense_cherenkov_pe);
                    dense_nsb_pe.swap(filtered_dense_nsb_pe);
                    dense_primary_dark_pe.swap(filtered_dense_primary_dark_pe);
                    dense_fired_cherenkov_pe.swap(
                        filtered_dense_fired_cherenkov_pe);
                    dense_fired_nsb_pe.swap(filtered_dense_fired_nsb_pe);
                    dense_fired_dark_pe.swap(filtered_dense_fired_dark_pe);
                    dense_gap_lost_pe.swap(filtered_dense_gap_lost_pe);
                    dense_saturation_lost_pe.swap(
                        filtered_dense_saturation_lost_pe);
                }
                if (output_cfg.write_pixel_time_stats) {
                    dense_time_mean_ns.swap(filtered_dense_time_mean_ns);
                    dense_time_rms_ns.swap(filtered_dense_time_rms_ns);
                }
            }
            if (detector_cfg.enabled) {
                telescope_trigger_rows.erase(
                    std::remove_if(
                        telescope_trigger_rows.begin(),
                        telescope_trigger_rows.end(),
                        [&electronics_events](const TelescopeTriggerRow& row) {
                            const auto event = electronics_events.find({
                                static_cast<int>(row.event_id),
                                static_cast<int>(row.telescope_id)});
                            return event == electronics_events.end() ||
                                   !event->second.selected_for_output;
                        }),
                    telescope_trigger_rows.end());
            }
        }

        struct EventRow {
            std::int32_t event_index;
            std::int64_t event_id;
        };
        struct CorsikaEventRow {
            std::int64_t event_id;
            std::int32_t shower_event_id;
            std::int32_t array_id;
            std::int32_t has_explicit_area_weight;
            std::int32_t primary_type;
            double energy_gev;
            double theta_deg;
            double phi_deg;
            double azimuth_north_to_east_deg;
            double core_x_north_m;
            double core_y_west_m;
            double array_rotation_deg;
            double array_time_offset_ns;
            double area_weight_m2;
            double h_first_int_m;
            double x_max_g_cm2;
            double h_max_m;
            double starting_grammage_g_cm2;
            double ground_gammas;
            double ground_electrons;
            double ground_hadrons;
            double ground_muons;
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
            double h_first_int_m;
            double x_max_g_cm2;
            double h_max_m;
            double starting_grammage_g_cm2;
            double ground_gammas;
            double ground_electrons;
            double ground_hadrons;
            double ground_muons;
        };
        std::vector<EventRow> event_rows;
        std::vector<CorsikaEventRow> corsika_event_rows;
        std::set<int> event_ids;
        for (const auto& image : image_rows) {
            event_ids.insert(static_cast<int>(image.event_id));
        }
        int event_index = 0;
        for (const int event_id : event_ids) {
            event_rows.push_back(EventRow{
                event_index++,
                static_cast<std::int64_t>(event_id),
            });
            const OutputEventMetadata event_meta = outputEventMetadata(
                event_id, source_runtime_cfg.event_id_mode, metadata);
            const int shower_event = event_meta.shower_event;
            const int array_id = event_meta.array_id;
            auto event_it = std::find_if(
                metadata.events.begin(), metadata.events.end(),
                [shower_event](const EventIOEventHeader& event) {
                    return event.shower_event_id == shower_event;
                });
            if (event_it != metadata.events.end()) {
                corsika_event_rows.push_back(CorsikaEventRow{
                    static_cast<std::int64_t>(event_id),
                    static_cast<std::int32_t>(shower_event),
                    static_cast<std::int32_t>(array_id),
                    static_cast<std::int32_t>(event_meta.has_explicit_area_weight),
                    static_cast<std::int32_t>(event_it->primary_type),
                    event_it->energy_gev,
                    event_it->theta_deg,
                    event_it->phi_deg,
                    event_it->azimuth_north_to_east_deg,
                    event_meta.core_x_north_m,
                    event_meta.core_y_west_m,
                    event_it->array_rotation_deg,
                    event_meta.array_time_offset_ns,
                    event_meta.area_weight_m2,
                    event_it->h_first_int_m,
                    event_it->x_max_g_cm2,
                    event_it->h_max_m,
                    event_it->starting_grammage_g_cm2,
                    event_it->ground_gammas,
                    event_it->ground_electrons,
                    event_it->ground_hadrons,
                    event_it->ground_muons,
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
            H5Tinsert(corsika_event_type, "has_explicit_area_weight",
                      HOFFSET(CorsikaEventRow, has_explicit_area_weight),
                      H5T_NATIVE_INT32);
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
            H5Tinsert(corsika_event_type, "array_time_offset_ns",
                      HOFFSET(CorsikaEventRow, array_time_offset_ns), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_event_type, "area_weight_m2",
                      HOFFSET(CorsikaEventRow, area_weight_m2), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_event_type, "h_first_int_m",
                      HOFFSET(CorsikaEventRow, h_first_int_m), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_event_type, "x_max_g_cm2",
                      HOFFSET(CorsikaEventRow, x_max_g_cm2), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_event_type, "h_max_m",
                      HOFFSET(CorsikaEventRow, h_max_m), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_event_type, "starting_grammage_g_cm2",
                      HOFFSET(CorsikaEventRow, starting_grammage_g_cm2), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_event_type, "ground_gammas",
                      HOFFSET(CorsikaEventRow, ground_gammas), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_event_type, "ground_electrons",
                      HOFFSET(CorsikaEventRow, ground_electrons), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_event_type, "ground_hadrons",
                      HOFFSET(CorsikaEventRow, ground_hadrons), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_event_type, "ground_muons",
                      HOFFSET(CorsikaEventRow, ground_muons), H5T_NATIVE_DOUBLE);
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
                    event.h_first_int_m,
                    event.x_max_g_cm2,
                    event.h_max_m,
                    event.starting_grammage_g_cm2,
                    event.ground_gammas,
                    event.ground_electrons,
                    event.ground_hadrons,
                    event.ground_muons,
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
            H5Tinsert(corsika_shower_type, "h_first_int_m",
                      HOFFSET(CorsikaShowerRow, h_first_int_m), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_shower_type, "x_max_g_cm2",
                      HOFFSET(CorsikaShowerRow, x_max_g_cm2), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_shower_type, "h_max_m",
                      HOFFSET(CorsikaShowerRow, h_max_m), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_shower_type, "starting_grammage_g_cm2",
                      HOFFSET(CorsikaShowerRow, starting_grammage_g_cm2), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_shower_type, "ground_gammas",
                      HOFFSET(CorsikaShowerRow, ground_gammas), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_shower_type, "ground_electrons",
                      HOFFSET(CorsikaShowerRow, ground_electrons), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_shower_type, "ground_hadrons",
                      HOFFSET(CorsikaShowerRow, ground_hadrons), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_shower_type, "ground_muons",
                      HOFFSET(CorsikaShowerRow, ground_muons), H5T_NATIVE_DOUBLE);
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
                  HOFFSET(ImageIndexRow, time_mean_ns), H5T_NATIVE_DOUBLE);
        H5Tinsert(image_type, "time_rms_ns",
                  HOFFSET(ImageIndexRow, time_rms_ns), H5T_NATIVE_DOUBLE);
        H5Tinsert(image_type, "time_first_ns",
                  HOFFSET(ImageIndexRow, time_first_ns), H5T_NATIVE_DOUBLE);
        H5Tinsert(image_type, "reference_time_ns",
                  HOFFSET(ImageIndexRow, reference_time_ns), H5T_NATIVE_DOUBLE);
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
            H5Tinsert(sparse_type, "primary_cherenkov_pe",
                      HOFFSET(SparsePixelRow, primary_cherenkov_pe), H5T_NATIVE_FLOAT);
            H5Tinsert(sparse_type, "primary_nsb_pe",
                      HOFFSET(SparsePixelRow, primary_nsb_pe), H5T_NATIVE_FLOAT);
            H5Tinsert(sparse_type, "primary_dark_pe",
                      HOFFSET(SparsePixelRow, primary_dark_pe), H5T_NATIVE_FLOAT);
            H5Tinsert(sparse_type, "fired_cherenkov_pe",
                      HOFFSET(SparsePixelRow, fired_cherenkov_pe), H5T_NATIVE_FLOAT);
            H5Tinsert(sparse_type, "fired_nsb_pe",
                      HOFFSET(SparsePixelRow, fired_nsb_pe), H5T_NATIVE_FLOAT);
            H5Tinsert(sparse_type, "fired_dark_pe",
                      HOFFSET(SparsePixelRow, fired_dark_pe), H5T_NATIVE_FLOAT);
            H5Tinsert(sparse_type, "gap_lost_pe",
                      HOFFSET(SparsePixelRow, gap_lost_pe), H5T_NATIVE_FLOAT);
            H5Tinsert(sparse_type, "saturation_lost_pe",
                      HOFFSET(SparsePixelRow, saturation_lost_pe), H5T_NATIVE_FLOAT);
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
            if (!dense_cherenkov_pe.empty()) {
                writePlain2D(dense_group,
                             "primary_cherenkov_pe",
                             H5T_NATIVE_FLOAT,
                             dense_cherenkov_pe,
                             static_cast<hsize_t>(n_images),
                             static_cast<hsize_t>(n_pixels));
                writePlain2D(dense_group,
                             "primary_nsb_pe",
                             H5T_NATIVE_FLOAT,
                             dense_nsb_pe,
                             static_cast<hsize_t>(n_images),
                             static_cast<hsize_t>(n_pixels));
                writePlain2D(dense_group, "primary_dark_pe",
                             H5T_NATIVE_FLOAT, dense_primary_dark_pe,
                             static_cast<hsize_t>(n_images),
                             static_cast<hsize_t>(n_pixels));
                writePlain2D(dense_group, "fired_cherenkov_pe",
                             H5T_NATIVE_FLOAT, dense_fired_cherenkov_pe,
                             static_cast<hsize_t>(n_images),
                             static_cast<hsize_t>(n_pixels));
                writePlain2D(dense_group, "fired_nsb_pe",
                             H5T_NATIVE_FLOAT, dense_fired_nsb_pe,
                             static_cast<hsize_t>(n_images),
                             static_cast<hsize_t>(n_pixels));
                writePlain2D(dense_group, "fired_dark_pe",
                             H5T_NATIVE_FLOAT, dense_fired_dark_pe,
                             static_cast<hsize_t>(n_images),
                             static_cast<hsize_t>(n_pixels));
                writePlain2D(dense_group, "gap_lost_pe",
                             H5T_NATIVE_FLOAT, dense_gap_lost_pe,
                             static_cast<hsize_t>(n_images),
                             static_cast<hsize_t>(n_pixels));
                writePlain2D(dense_group, "saturation_lost_pe",
                             H5T_NATIVE_FLOAT, dense_saturation_lost_pe,
                             static_cast<hsize_t>(n_images),
                             static_cast<hsize_t>(n_pixels));
            }
            if (output_cfg.write_pixel_time_stats) {
                writePlain2D(dense_group,
                             "time_mean_ns",
                             H5T_NATIVE_FLOAT,
                             dense_time_mean_ns,
                             static_cast<hsize_t>(n_images),
                             static_cast<hsize_t>(n_pixels));
                writePlain2D(dense_group,
                             "time_rms_ns",
                             H5T_NATIVE_FLOAT,
                             dense_time_rms_ns,
                             static_cast<hsize_t>(n_images),
                             static_cast<hsize_t>(n_pixels));
            }
            H5Gclose(dense_group);
        }
        H5Gclose(images_group);

        if (waveform_cfg.enabled && output_cfg.hdf5_write_waveforms &&
            have_camera_axis) {
            std::vector<Hdf5WaveformImage> waveform_images;
            waveform_images.reserve(image_rows.size());
            for (const auto& image : image_rows) {
                waveform_images.push_back(Hdf5WaveformImage{
                    image.image_index,
                    static_cast<int>(image.event_id),
                    static_cast<int>(image.telescope_id),
                    static_cast<double>(image.time_first_ns),
                    static_cast<double>(image.time_mean_ns),
                });
            }
            writeHdf5Waveforms(file,
                               output_cfg,
                               waveform_cfg,
                               nsb_cfg,
                               pixel_id_axis,
                               waveform_images,
                               waveforms,
                               raw_waveform_hits,
                               electronics_events);
        }

        if (detector_cfg.enabled) {
            struct PrimaryPeRow {
                std::int64_t event_id;
                std::int32_t telescope_id;
                double reference_time_ns;
                std::int32_t pixel_id;
                double time_ns;
                double global_time_ns;
                double sensor_x_m;
                double sensor_y_m;
                double wavelength_nm;
                double primary_pe;
                std::int32_t origin;
            };
            struct FiredPeRow {
                std::int64_t event_id;
                std::int32_t telescope_id;
                double reference_time_ns;
                std::int32_t pixel_id;
                double time_ns;
                double global_time_ns;
                std::int32_t channel_id;
                std::int32_t microcell_id;
                double fired_pe;
                double charge_factor;
                double time_jitter_ns;
                std::int32_t origin;
            };
            std::vector<PrimaryPeRow> primary_rows;
            std::vector<FiredPeRow> fired_rows;
            for (const auto& item : electronics_events) {
                const auto& event = item.second;
                if (!event.selected_for_output) continue;
                if (detector_cfg.save_primary_sequence) {
                    for (const auto& hit : event.detector.primary_hits) {
                        primary_rows.push_back(PrimaryPeRow{
                            hit.event_id, hit.telescope_id,
                            event.reference_time_ns,
                            pixel_id_axis.at(static_cast<std::size_t>(hit.pixel_id)),
                            hit.time_ns,
                            event.reference_time_ns + hit.time_ns,
                            hit.sensor_x_m, hit.sensor_y_m, hit.wavelength_nm,
                            hit.primary_pe, static_cast<std::int32_t>(hit.origin),
                        });
                    }
                }
                if (detector_cfg.save_fired_sequence) {
                    for (const auto& hit : event.detector.fired_hits) {
                        fired_rows.push_back(FiredPeRow{
                            hit.event_id, hit.telescope_id,
                            event.reference_time_ns,
                            pixel_id_axis.at(static_cast<std::size_t>(hit.pixel_id)),
                            hit.time_ns,
                            event.reference_time_ns + hit.time_ns,
                            hit.channel_id, hit.microcell_id, hit.fired_pe,
                            hit.charge_factor, hit.time_jitter_ns,
                            static_cast<std::int32_t>(hit.origin),
                        });
                    }
                }
            }
            hid_t electronics_group = H5Gcreate2(
                file, "electronics", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            if (detector_cfg.save_primary_sequence) {
                hid_t type = H5Tcreate(H5T_COMPOUND, sizeof(PrimaryPeRow));
                H5Tinsert(type, "event_id", HOFFSET(PrimaryPeRow, event_id), H5T_NATIVE_INT64);
                H5Tinsert(type, "telescope_id", HOFFSET(PrimaryPeRow, telescope_id), H5T_NATIVE_INT32);
                H5Tinsert(type, "reference_time_ns", HOFFSET(PrimaryPeRow, reference_time_ns), H5T_NATIVE_DOUBLE);
                H5Tinsert(type, "pixel_id", HOFFSET(PrimaryPeRow, pixel_id), H5T_NATIVE_INT32);
                H5Tinsert(type, "time_ns", HOFFSET(PrimaryPeRow, time_ns), H5T_NATIVE_DOUBLE);
                H5Tinsert(type, "global_time_ns", HOFFSET(PrimaryPeRow, global_time_ns), H5T_NATIVE_DOUBLE);
                H5Tinsert(type, "sensor_x_m", HOFFSET(PrimaryPeRow, sensor_x_m), H5T_NATIVE_DOUBLE);
                H5Tinsert(type, "sensor_y_m", HOFFSET(PrimaryPeRow, sensor_y_m), H5T_NATIVE_DOUBLE);
                H5Tinsert(type, "wavelength_nm", HOFFSET(PrimaryPeRow, wavelength_nm), H5T_NATIVE_DOUBLE);
                H5Tinsert(type, "primary_pe", HOFFSET(PrimaryPeRow, primary_pe), H5T_NATIVE_DOUBLE);
                H5Tinsert(type, "origin", HOFFSET(PrimaryPeRow, origin), H5T_NATIVE_INT32);
                writeCompound1D(electronics_group, "primary_pe_hits", type, primary_rows);
                H5Tclose(type);
            }
            if (detector_cfg.save_fired_sequence) {
                hid_t type = H5Tcreate(H5T_COMPOUND, sizeof(FiredPeRow));
                H5Tinsert(type, "event_id", HOFFSET(FiredPeRow, event_id), H5T_NATIVE_INT64);
                H5Tinsert(type, "telescope_id", HOFFSET(FiredPeRow, telescope_id), H5T_NATIVE_INT32);
                H5Tinsert(type, "reference_time_ns", HOFFSET(FiredPeRow, reference_time_ns), H5T_NATIVE_DOUBLE);
                H5Tinsert(type, "pixel_id", HOFFSET(FiredPeRow, pixel_id), H5T_NATIVE_INT32);
                H5Tinsert(type, "time_ns", HOFFSET(FiredPeRow, time_ns), H5T_NATIVE_DOUBLE);
                H5Tinsert(type, "global_time_ns", HOFFSET(FiredPeRow, global_time_ns), H5T_NATIVE_DOUBLE);
                H5Tinsert(type, "channel_id", HOFFSET(FiredPeRow, channel_id), H5T_NATIVE_INT32);
                H5Tinsert(type, "microcell_id", HOFFSET(FiredPeRow, microcell_id), H5T_NATIVE_INT32);
                H5Tinsert(type, "fired_pe", HOFFSET(FiredPeRow, fired_pe), H5T_NATIVE_DOUBLE);
                H5Tinsert(type, "charge_factor", HOFFSET(FiredPeRow, charge_factor), H5T_NATIVE_DOUBLE);
                H5Tinsert(type, "time_jitter_ns", HOFFSET(FiredPeRow, time_jitter_ns), H5T_NATIVE_DOUBLE);
                H5Tinsert(type, "origin", HOFFSET(FiredPeRow, origin), H5T_NATIVE_INT32);
                writeCompound1D(electronics_group, "fired_pe_hits", type, fired_rows);
                H5Tclose(type);
            }
            H5Gclose(electronics_group);
        }

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
                  HOFFSET(TelescopeTriggerRow, trigger_time_ns), H5T_NATIVE_DOUBLE);
        H5Tinsert(telescope_trigger_type, "trigger_first_time_ns",
                  HOFFSET(TelescopeTriggerRow, trigger_first_time_ns),
                  H5T_NATIVE_DOUBLE);
        H5Tinsert(telescope_trigger_type, "trigger_max_multiplicity_time_ns",
                  HOFFSET(TelescopeTriggerRow,
                          trigger_max_multiplicity_time_ns),
                  H5T_NATIVE_DOUBLE);
        H5Tinsert(telescope_trigger_type, "geometric_delay_ns",
                  HOFFSET(TelescopeTriggerRow, geometric_delay_ns),
                  H5T_NATIVE_DOUBLE);
        H5Tinsert(telescope_trigger_type, "coincidence_time_ns",
                  HOFFSET(TelescopeTriggerRow, coincidence_time_ns),
                  H5T_NATIVE_DOUBLE);
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
            if (output_cfg.whiteboard_emitter_info) {
                H5Tinsert(hit_type, "has_emitter",
                          HOFFSET(WhiteboardHdf5Row, has_emitter), H5T_NATIVE_UINT8);
                H5Tinsert(hit_type, "emitter_mass_gev",
                          HOFFSET(WhiteboardHdf5Row, emitter_mass_gev), H5T_NATIVE_FLOAT);
                H5Tinsert(hit_type, "emitter_charge",
                          HOFFSET(WhiteboardHdf5Row, emitter_charge), H5T_NATIVE_FLOAT);
                H5Tinsert(hit_type, "emitter_energy_gev",
                          HOFFSET(WhiteboardHdf5Row, emitter_energy_gev), H5T_NATIVE_FLOAT);
                H5Tinsert(hit_type, "emitter_time_ns",
                          HOFFSET(WhiteboardHdf5Row, emitter_time_ns), H5T_NATIVE_FLOAT);
            }
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

} // namespace lact

#endif // LACT_HAS_HDF5
