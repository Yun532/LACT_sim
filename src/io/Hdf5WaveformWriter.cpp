#include "io/Hdf5WaveformWriter.hpp"

#include "app/NsbResponseSampler.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace lact {
namespace {

struct SparseWaveformSample {
    std::int32_t image_index = 0;
    std::int32_t pixel_id = 0;
    std::int32_t time_bin = 0;
    std::int32_t photon_count = 0;
    float pe = 0.0f;
    float cherenkov_pe = 0.0f;
    float nsb_pe = 0.0f;
};

void check(herr_t status, const std::string& message)
{
    if (status < 0) {
        throw std::runtime_error("HDF5 waveform write failed: " + message);
    }
}

void writeStringAttribute(hid_t object,
                          const std::string& name,
                          const std::string& value)
{
    hid_t space = H5Screate(H5S_SCALAR);
    hid_t type = H5Tcopy(H5T_C_S1);
    H5Tset_size(type, std::max<std::size_t>(1, value.size() + 1));
    H5Tset_strpad(type, H5T_STR_NULLTERM);
    hid_t attribute = H5Acreate2(
        object, name.c_str(), type, space, H5P_DEFAULT, H5P_DEFAULT);
    check(H5Awrite(attribute, type, value.c_str()), "attribute " + name);
    H5Aclose(attribute);
    H5Tclose(type);
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
    hid_t dataset = H5Dcreate2(
        group, name.c_str(), type, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (!values.empty()) {
        check(H5Dwrite(dataset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data()),
              "dataset " + name);
    }
    H5Dclose(dataset);
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
    hid_t dataset = H5Dcreate2(
        group, name.c_str(), type, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (!values.empty()) {
        check(H5Dwrite(dataset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data()),
              "dataset " + name);
    }
    H5Dclose(dataset);
    H5Sclose(space);
}

void writeSparseSamples(hid_t group,
                        const std::vector<SparseWaveformSample>& samples)
{
    hid_t type = H5Tcreate(H5T_COMPOUND, sizeof(SparseWaveformSample));
    H5Tinsert(type, "image_index",
              HOFFSET(SparseWaveformSample, image_index), H5T_NATIVE_INT32);
    H5Tinsert(type, "pixel_id",
              HOFFSET(SparseWaveformSample, pixel_id), H5T_NATIVE_INT32);
    H5Tinsert(type, "time_bin",
              HOFFSET(SparseWaveformSample, time_bin), H5T_NATIVE_INT32);
    H5Tinsert(type, "photon_count",
              HOFFSET(SparseWaveformSample, photon_count), H5T_NATIVE_INT32);
    H5Tinsert(type, "pe", HOFFSET(SparseWaveformSample, pe), H5T_NATIVE_FLOAT);
    H5Tinsert(type, "cherenkov_pe",
              HOFFSET(SparseWaveformSample, cherenkov_pe), H5T_NATIVE_FLOAT);
    H5Tinsert(type, "nsb_pe",
              HOFFSET(SparseWaveformSample, nsb_pe), H5T_NATIVE_FLOAT);
    hsize_t dims[1] = {samples.size()};
    hid_t space = H5Screate_simple(1, dims, nullptr);
    hid_t dataset = H5Dcreate2(
        group, "samples", type, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (!samples.empty()) {
        check(H5Dwrite(dataset, type, H5S_ALL, H5S_ALL,
                       H5P_DEFAULT, samples.data()),
              "dataset samples");
    }
    H5Dclose(dataset);
    H5Sclose(space);
    H5Tclose(type);
}

std::size_t binCount(const WaveformOutputConfig& waveform)
{
    return static_cast<std::size_t>(std::ceil(
        (waveform.time_window_end_ns - waveform.time_window_start_ns) /
        waveform.time_bin_width_ns));
}

int binForTime(const WaveformOutputConfig& waveform, double time_ns)
{
    if (time_ns < waveform.time_window_start_ns ||
        time_ns >= waveform.time_window_end_ns) {
        return -1;
    }
    const int bin = static_cast<int>(std::floor(
        (time_ns - waveform.time_window_start_ns) /
        waveform.time_bin_width_ns));
    return bin >= 0 && static_cast<std::size_t>(bin) < binCount(waveform)
               ? bin
               : -1;
}

bool usesImageReference(const WaveformOutputConfig& waveform)
{
    return waveform.time_reference == "image_first" ||
           waveform.time_reference == "image_mean";
}

} // namespace

void writeHdf5Waveforms(
    hid_t file,
    const CorsikaTraceOutputConfig& output,
    const WaveformOutputConfig& waveform,
    const NsbConfig& nsb,
    const std::vector<std::int32_t>& pixel_id_axis,
    const std::vector<Hdf5WaveformImage>& images,
    const std::map<WaveformKey, WaveformPixelAccumulator>& waveforms,
    const std::vector<RawWaveformHit>& raw_hits)
{
    if (!waveform.enabled || !output.hdf5_write_waveforms ||
        pixel_id_axis.empty()) {
        return;
    }
    const std::size_t n_images = images.size();
    const std::size_t n_pixels = pixel_id_axis.size();
    const std::size_t n_bins = binCount(waveform);
    const std::vector<int> nsb_pixel_ids(
        pixel_id_axis.begin(), pixel_id_axis.end());

    std::vector<double> time_edges(n_bins + 1, 0.0);
    std::vector<double> time_centers(n_bins, 0.0);
    for (std::size_t i = 0; i <= n_bins; ++i) {
        time_edges[i] = waveform.time_window_start_ns +
                        static_cast<double>(i) * waveform.time_bin_width_ns;
    }
    for (std::size_t i = 0; i < n_bins; ++i) {
        time_centers[i] = 0.5 * (time_edges[i] + time_edges[i + 1]);
    }

    std::map<SummaryKey, std::size_t> image_row_by_key;
    std::map<std::int32_t, std::size_t> row_by_image_index;
    std::map<int, std::size_t> pixel_to_col;
    std::vector<double> reference_time_ns(n_images, 0.0);
    for (std::size_t col = 0; col < pixel_id_axis.size(); ++col) {
        pixel_to_col[pixel_id_axis[col]] = col;
    }
    for (std::size_t row = 0; row < images.size(); ++row) {
        image_row_by_key[{images[row].event_id, images[row].telescope_id}] = row;
        row_by_image_index[images[row].image_index] = row;
        if (waveform.time_reference == "image_first" &&
            std::isfinite(images[row].time_first_ns)) {
            reference_time_ns[row] = images[row].time_first_ns;
        } else if (waveform.time_reference == "image_mean" &&
                   std::isfinite(images[row].time_mean_ns)) {
            reference_time_ns[row] = images[row].time_mean_ns;
        }
    }

    std::map<std::size_t, SparseWaveformSample> samples_by_index;
    auto addSample = [&](std::size_t row,
                         std::size_t col,
                         std::size_t bin,
                         std::int32_t photon_count,
                         float pe,
                         float cherenkov_pe,
                         float nsb_pe) {
        if (row >= n_images || col >= n_pixels || bin >= n_bins) return;
        const std::size_t flat = (row * n_bins + bin) * n_pixels + col;
        auto& sample = samples_by_index[flat];
        sample.image_index = images[row].image_index;
        sample.pixel_id = pixel_id_axis[col];
        sample.time_bin = static_cast<std::int32_t>(bin);
        sample.photon_count += photon_count;
        sample.pe += pe;
        sample.cherenkov_pe += cherenkov_pe;
        sample.nsb_pe += nsb_pe;
    };

    for (const auto& item : waveforms) {
        const auto& value = item.second;
        const auto image = image_row_by_key.find(
            {value.event_id, value.telescope_id});
        const auto pixel = pixel_to_col.find(value.pixel_id);
        if (image == image_row_by_key.end() || pixel == pixel_to_col.end() ||
            value.time_bin < 0 ||
            static_cast<std::size_t>(value.time_bin) >= n_bins) {
            continue;
        }
        if (waveform.source == "photon_count") {
            addSample(image->second, pixel->second,
                      static_cast<std::size_t>(value.time_bin),
                      static_cast<std::int32_t>(value.photon_count),
                      0.0f, 0.0f, 0.0f);
        } else if (waveform.source == "pe") {
            const float pe = static_cast<float>(value.pe);
            addSample(image->second, pixel->second,
                      static_cast<std::size_t>(value.time_bin),
                      0, pe, pe, 0.0f);
        }
    }
    if (usesImageReference(waveform)) {
        for (const auto& hit : raw_hits) {
            const auto image = image_row_by_key.find(
                {hit.event_id, hit.telescope_id});
            const auto pixel = pixel_to_col.find(hit.pixel_id);
            if (image == image_row_by_key.end() || pixel == pixel_to_col.end()) {
                continue;
            }
            const int bin = binForTime(
                waveform, hit.time_ns - reference_time_ns[image->second]);
            if (bin < 0) continue;
            if (waveform.source == "photon_count") {
                addSample(image->second, pixel->second,
                          static_cast<std::size_t>(bin),
                          static_cast<std::int32_t>(hit.photon_count),
                          0.0f, 0.0f, 0.0f);
            } else if (waveform.source == "pe") {
                const float pe = static_cast<float>(hit.pe);
                addSample(image->second, pixel->second,
                          static_cast<std::size_t>(bin),
                          0, pe, pe, 0.0f);
            }
        }
    }
    if (waveform.source == "pe" && nsb.enabled &&
        nsb.rate_pe_per_ns_per_pixel > 0.0) {
        for (std::size_t row = 0; row < images.size(); ++row) {
            const auto realization = generateNsbRealization(
                nsb, images[row].event_id, images[row].telescope_id,
                nsb_pixel_ids, n_bins, waveform.time_bin_width_ns);
            for (const auto& sample : realization.samples) {
                addSample(row, sample.pixel_col, sample.time_bin,
                          0, sample.pe, 0.0f, sample.pe);
            }
        }
    }

    std::vector<SparseWaveformSample> samples;
    samples.reserve(samples_by_index.size());
    for (const auto& item : samples_by_index) {
        samples.push_back(item.second);
    }

    hid_t group = H5Gcreate2(
        file, "waveforms", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    writeStringAttribute(group, "source", waveform.source);
    writeStringAttribute(group, "time_reference", waveform.time_reference);
    writeStringAttribute(group, "storage",
                         output.hdf5_waveform_storage == "sparse"
                             ? "sparse_coo"
                             : "dense");
    writeStringAttribute(group, "shape", "image_index,time_bin,pixel_id_axis");
    writeStringAttribute(
        group, "note",
        output.hdf5_waveform_storage == "sparse"
            ? "sparse proxy waveform at camera/collector output; real electronics is not modeled"
            : "proxy waveform accumulated at camera/collector output; real electronics waveform is not modeled");
    writePlain1D(group, "pixel_id_axis", H5T_NATIVE_INT32, pixel_id_axis);
    writePlain1D(group, "time_edges_ns", H5T_NATIVE_DOUBLE, time_edges);
    writePlain1D(group, "time_centers_ns", H5T_NATIVE_DOUBLE, time_centers);
    writePlain1D(group, "reference_time_ns", H5T_NATIVE_DOUBLE, reference_time_ns);

    if (output.hdf5_waveform_storage == "sparse") {
        writeSparseSamples(group, samples);
    } else {
        const std::size_t size = n_images * n_bins * n_pixels;
        std::vector<std::int32_t> photon_count;
        std::vector<float> pe;
        std::vector<float> cherenkov_pe;
        std::vector<float> nsb_pe;
        if (waveform.source == "photon_count") {
            photon_count.assign(size, 0);
        } else if (waveform.source == "pe") {
            pe.assign(size, 0.0f);
            if (output.hdf5_write_components) {
                cherenkov_pe.assign(size, 0.0f);
                nsb_pe.assign(size, 0.0f);
            }
        }
        for (const auto& sample : samples) {
            const auto row = row_by_image_index.at(sample.image_index);
            const auto col = pixel_to_col.at(sample.pixel_id);
            const std::size_t flat =
                (row * n_bins + static_cast<std::size_t>(sample.time_bin)) *
                    n_pixels + col;
            if (waveform.source == "photon_count") {
                photon_count[flat] = sample.photon_count;
            } else if (waveform.source == "pe") {
                pe[flat] = sample.pe;
                if (output.hdf5_write_components) {
                    cherenkov_pe[flat] = sample.cherenkov_pe;
                    nsb_pe[flat] = sample.nsb_pe;
                }
            }
        }
        if (waveform.source == "photon_count") {
            writePlain3D(group, "photon_count", H5T_NATIVE_INT32, photon_count,
                         n_images, n_bins, n_pixels);
        } else if (waveform.source == "pe") {
            if (output.hdf5_write_components) {
                writePlain3D(group, "cherenkov_pe", H5T_NATIVE_FLOAT,
                             cherenkov_pe, n_images, n_bins, n_pixels);
                writePlain3D(group, "nsb_pe", H5T_NATIVE_FLOAT,
                             nsb_pe, n_images, n_bins, n_pixels);
            }
            writePlain3D(group, "pe", H5T_NATIVE_FLOAT, pe,
                         n_images, n_bins, n_pixels);
        }
    }
    H5Gclose(group);
}

} // namespace lact
