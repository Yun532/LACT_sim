#include "io/PhotonCsvSource.hpp"
#include <algorithm>
#include <cctype>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>

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

} // namespace

PhotonCsvSource::PhotonCsvSource(const PhotonCsvConfig& cfg)
    : cfg_(cfg) {
    if (cfg_.csv_path.empty()) {
        throw std::runtime_error("PhotonCsvSource: csv_path is required");
    }
    load();
}

void PhotonCsvSource::reset() {
    index_ = 0;
}

bool PhotonCsvSource::next(PhotonBunch& out) {
    if (index_ >= rows_.size()) {
        return false;
    }
    out = rows_[index_++];
    return true;
}

void PhotonCsvSource::load() {
    std::ifstream ifs(cfg_.csv_path);
    if (!ifs) {
        throw std::runtime_error("failed to open photon CSV: " + cfg_.csv_path);
    }

    std::string line;
    if (!std::getline(ifs, line)) {
        throw std::runtime_error("photon CSV is empty: " + cfg_.csv_path);
    }
    auto header_cells = splitCsv(line);
    std::map<std::string, int> header;
    for (int i = 0; i < static_cast<int>(header_cells.size()); ++i) {
        header[lowerCopy(header_cells[i])] = i;
    }

    int line_no = 1;
    while (std::getline(ifs, line)) {
        ++line_no;
        if (trim(line).empty()) continue;

        auto cells = splitCsv(line);
        PhotonBunch bunch;
        bunch.photon.pos = {
            std::stod(requiredCell(cells, header, "x_m", line_no)),
            std::stod(requiredCell(cells, header, "y_m", line_no)),
            std::stod(requiredCell(cells, header, "z_m", line_no))
        };
        bunch.photon.dir = {
            std::stod(requiredCell(cells, header, "dir_x", line_no)),
            std::stod(requiredCell(cells, header, "dir_y", line_no)),
            std::stod(requiredCell(cells, header, "dir_z", line_no))
        };
        bunch.photon.normalizeDirection();
        bunch.photon.time_ns =
            optionalDouble(cells, header, "time_ns", cfg_.default_time_ns);
        bunch.photon.wavelength_nm =
            optionalDouble(cells, header, "wavelength_nm", cfg_.default_wavelength_nm);
        bunch.photon.weight =
            optionalDouble(cells, header, "weight", cfg_.default_weight);
        bunch.multiplicity =
            optionalDouble(cells, header, "multiplicity", cfg_.default_multiplicity);
        bunch.event_id =
            optionalInt(cells, header, "event_id", cfg_.default_event_id);
        bunch.telescope_id =
            optionalInt(cells, header, "telescope_id", cfg_.default_telescope_id);

        if (cfg_.filter_telescope_id &&
            bunch.telescope_id != cfg_.selected_telescope_id) {
            continue;
        }
        if (cfg_.filter_event_id &&
            bunch.event_id != cfg_.selected_event_id) {
            continue;
        }

        rows_.push_back(bunch);
    }

    if (rows_.empty()) {
        throw std::runtime_error("photon CSV has no data rows: " + cfg_.csv_path);
    }
}
