#include "io/PhotonCsvSource.hpp"
#include <algorithm>
#include <cctype>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <cstdint>
#include <cmath>
#include <vector>

namespace {

std::string trim(const std::string& s) {
    auto first = std::find_if_not(s.begin(), s.end(), [](unsigned char c) {
        return std::isspace(c);
    });
    auto last = std::find_if_not(s.rbegin(), s.rend(), [](unsigned char c) {
        return std::isspace(c);
    }).base();
    if (first >= last) return "";
    return std::string(first, last);
}

std::string lowerCopy(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return s;
}

std::vector<std::string> splitCsv(const std::string& line) {
    std::stringstream ss(line);
    std::string cell;
    std::vector<std::string> cells;
    while (std::getline(ss, cell, ',')) {
        cells.push_back(trim(cell));
    }
    return cells;
}

int headerIndex(const std::map<std::string, int>& header, const std::string& key) {
    auto it = header.find(lowerCopy(key));
    return it == header.end() ? -1 : it->second;
}

std::string optionalCell(const std::vector<std::string>& cells,
                         const std::map<std::string, int>& header,
                         const std::string& key) {
    int idx = headerIndex(header, key);
    if (idx < 0 || static_cast<std::size_t>(idx) >= cells.size()) return "";
    return cells[idx];
}

std::string requiredCell(const std::vector<std::string>& cells,
                         const std::map<std::string, int>& header,
                         const std::string& key,
                         int line_no) {
    std::string value = optionalCell(cells, header, key);
    if (value.empty()) {
        throw std::runtime_error("PhotonCsvSource line " + std::to_string(line_no) +
                                 " missing required column: " + key);
    }
    return value;
}

double optionalDouble(const std::vector<std::string>& cells,
                      const std::map<std::string, int>& header,
                      const std::string& key,
                      double fallback) {
    std::string value = optionalCell(cells, header, key);
    return value.empty() ? fallback : std::stod(value);
}

int optionalInt(const std::vector<std::string>& cells,
                const std::map<std::string, int>& header,
                const std::string& key,
                int fallback) {
    std::string value = optionalCell(cells, header, key);
    return value.empty() ? fallback : std::stoi(value);
}

std::uint64_t optionalUInt64(const std::vector<std::string>& cells,
                             const std::map<std::string, int>& header,
                             const std::string& key,
                             std::uint64_t fallback) {
    std::string value = optionalCell(cells, header, key);
    return value.empty() ? fallback : std::stoull(value);
}

bool optionalBool(const std::vector<std::string>& cells,
                  const std::map<std::string, int>& header,
                  const std::string& key,
                  bool fallback) {
    std::string value = lowerCopy(optionalCell(cells, header, key));
    if (value.empty()) return fallback;
    if (value == "1" || value == "true" || value == "yes" || value == "on") return true;
    if (value == "0" || value == "false" || value == "no" || value == "off") return false;
    throw std::runtime_error("PhotonCsvSource invalid boolean in column " + key +
                             ": " + value);
}

} // namespace

PhotonCsvSource::PhotonCsvSource(const PhotonCsvConfig& cfg)
    : cfg_(cfg) {
    if (cfg_.csv_path.empty()) {
        throw std::runtime_error("PhotonCsvSource: csv_path is required");
    }
    openAndPrime();
}

void PhotonCsvSource::reset() {
    openAndPrime();
}

bool PhotonCsvSource::next(PhotonBunch& out) {
    if (buffered_) {
        out = std::move(*buffered_);
        buffered_.reset();
        return true;
    }
    return readNext(out);
}

void PhotonCsvSource::openAndPrime() {
    buffered_.reset();
    input_.close();
    input_.clear();
    input_.open(cfg_.csv_path);
    if (!input_) {
        throw std::runtime_error("failed to open photon CSV: " + cfg_.csv_path);
    }

    std::string line;
    if (!std::getline(input_, line)) {
        throw std::runtime_error("photon CSV is empty: " + cfg_.csv_path);
    }
    auto header_cells = splitCsv(line);
    header_.clear();
    for (int i = 0; i < static_cast<int>(header_cells.size()); ++i) {
        header_[lowerCopy(header_cells[i])] = i;
    }
    line_no_ = 1;
    PhotonBunch first;
    if (!readNext(first)) {
        throw std::runtime_error("photon CSV has no data rows: " + cfg_.csv_path);
    }
    buffered_ = std::move(first);
}

bool PhotonCsvSource::readNext(PhotonBunch& out) {
    std::string line;
    while (std::getline(input_, line)) {
        ++line_no_;
        if (trim(line).empty()) continue;

        auto cells = splitCsv(line);
        PhotonBunch bunch;
        bunch.photon.pos = {
            std::stod(requiredCell(cells, header_, "x_m", line_no_)),
            std::stod(requiredCell(cells, header_, "y_m", line_no_)),
            std::stod(requiredCell(cells, header_, "z_m", line_no_))
        };
        bunch.photon.dir = {
            std::stod(requiredCell(cells, header_, "dir_x", line_no_)),
            std::stod(requiredCell(cells, header_, "dir_y", line_no_)),
            std::stod(requiredCell(cells, header_, "dir_z", line_no_))
        };
        const auto finiteVec = [](const Vec3& value) {
            return std::isfinite(value.x) &&
                   std::isfinite(value.y) &&
                   std::isfinite(value.z);
        };
        if (!finiteVec(bunch.photon.pos)) {
            throw std::runtime_error(
                "PhotonCsvSource line " + std::to_string(line_no_) +
                " has a non-finite position");
        }
        if (!finiteVec(bunch.photon.dir) ||
            !std::isfinite(bunch.photon.dir.norm2()) ||
            bunch.photon.dir.norm2() <= 1.0e-30) {
            throw std::runtime_error(
                "PhotonCsvSource line " + std::to_string(line_no_) +
                " has a non-finite or zero direction");
        }
        bunch.photon.normalizeDirection();
        bunch.photon.time_ns =
            optionalDouble(cells, header_, "time_ns", cfg_.default_time_ns);
        bunch.photon.wavelength_nm =
            optionalDouble(cells, header_, "wavelength_nm", cfg_.default_wavelength_nm);
        bunch.raw_wavelength_nm =
            optionalDouble(cells, header_, "raw_wavelength_nm",
                           bunch.photon.wavelength_nm);
        bunch.photon.weight =
            optionalDouble(cells, header_, "weight", cfg_.default_weight);
        bunch.photon.optical_efficiency_preapplied =
            optionalBool(cells, header_, "optical_efficiency_preapplied", false);
        bunch.multiplicity =
            optionalDouble(cells, header_, "multiplicity", cfg_.default_multiplicity);
        bunch.event_id =
            optionalInt(cells, header_, "event_id", cfg_.default_event_id);
        bunch.telescope_id =
            optionalInt(cells, header_, "telescope_id", cfg_.default_telescope_id);
        bunch.shower_event_id =
            optionalInt(cells, header_, "shower_event_id", bunch.event_id);
        bunch.array_id =
            optionalInt(cells, header_, "array_id", 0);
        bunch.source_bunch_index =
            optionalUInt64(cells, header_, "source_bunch_index",
                           static_cast<std::uint64_t>(line_no_ - 2));
        bunch.eventio_2d =
            optionalBool(cells, header_, "eventio_2d", cfg_.default_eventio_2d);
        bunch.emission_altitude_km =
            optionalDouble(cells, header_, "emission_altitude_km",
                           bunch.emission_altitude_km);
        const std::string origin = lowerCopy(optionalCell(cells, header_, "origin"));
        if (origin.empty() || origin == "cherenkov") {
            bunch.origin = PhotonOrigin::Cherenkov;
        } else if (origin == "nsb") {
            bunch.origin = PhotonOrigin::Nsb;
        } else if (origin == "dark") {
            bunch.origin = PhotonOrigin::Dark;
        } else {
            throw std::runtime_error(
                "PhotonCsvSource line " + std::to_string(line_no_) +
                " origin must be cherenkov, nsb, or dark");
        }

        if (!std::isfinite(bunch.photon.time_ns)) {
            throw std::runtime_error(
                "PhotonCsvSource line " + std::to_string(line_no_) +
                " has a non-finite time_ns");
        }
        if (!std::isfinite(bunch.photon.wavelength_nm) ||
            bunch.photon.wavelength_nm <= 0.0) {
            throw std::runtime_error(
                "PhotonCsvSource line " + std::to_string(line_no_) +
                " requires a finite, positive wavelength_nm");
        }
        if (!std::isfinite(bunch.raw_wavelength_nm)) {
            throw std::runtime_error(
                "PhotonCsvSource line " + std::to_string(line_no_) +
                " has a non-finite raw_wavelength_nm");
        }
        if (!std::isfinite(bunch.photon.weight) ||
            bunch.photon.weight < 0.0) {
            throw std::runtime_error(
                "PhotonCsvSource line " + std::to_string(line_no_) +
                " requires a finite, non-negative weight");
        }
        if (!std::isfinite(bunch.multiplicity) ||
            bunch.multiplicity < 0.0) {
            throw std::runtime_error(
                "PhotonCsvSource line " + std::to_string(line_no_) +
                " requires a finite, non-negative multiplicity");
        }
        const std::string emission_altitude_text =
            optionalCell(cells, header_, "emission_altitude_km");
        if (!emission_altitude_text.empty() &&
            !std::isfinite(bunch.emission_altitude_km)) {
            throw std::runtime_error(
                "PhotonCsvSource line " + std::to_string(line_no_) +
                " has a non-finite emission_altitude_km");
        }

        if (cfg_.filter_telescope_id &&
            bunch.telescope_id != cfg_.selected_telescope_id) {
            continue;
        }
        if (cfg_.filter_event_id &&
            bunch.event_id != cfg_.selected_event_id) {
            continue;
        }

        out = std::move(bunch);
        return true;
    }
    return false;
}
