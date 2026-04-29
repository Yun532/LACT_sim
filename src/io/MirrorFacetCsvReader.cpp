#include "io/MirrorFacetCsvReader.hpp"
#include <algorithm>
#include <cctype>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include "geometry/MirrorFacetValidation.hpp"

namespace {

std::string trim(const std::string& s) {
    auto first = std::find_if_not(s.begin(), s.end(), [](unsigned char c) {
        return std::isspace(c);
    });
    auto last = std::find_if_not(s.rbegin(), s.rend(), [](unsigned char c) {
        return std::isspace(c);
    }).base();

    if (first >= last) {
        return "";
    }
    return std::string(first, last);
}

std::vector<std::string> splitCsvLine(const std::string& line) {
    std::vector<std::string> cells;
    std::stringstream ss(line);
    std::string cell;
    while (std::getline(ss, cell, ',')) {
        cells.push_back(trim(cell));
    }
    if (!line.empty() && line.back() == ',') {
        cells.push_back("");
    }
    return cells;
}

std::string lowerCopy(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return s;
}

bool parseSurfaceType(const std::string& text, SurfaceType& out) {
    std::string s = lowerCopy(trim(text));
    if (s == "spherical") {
        out = SurfaceType::Spherical;
        return true;
    }
    if (s == "parabolic") {
        out = SurfaceType::Parabolic;
        return true;
    }
    if (s == "planar") {
        out = SurfaceType::Planar;
        return true;
    }
    if (s == "polynomial") {
        out = SurfaceType::Polynomial;
        return true;
    }
    return false;
}

bool parseApertureShape(const std::string& text, ApertureShape& out) {
    std::string s = lowerCopy(trim(text));
    if (s == "circular") {
        out = ApertureShape::Circular;
        return true;
    }
    if (s == "hexagon") {
        out = ApertureShape::Hexagon;
        return true;
    }
    if (s == "square") {
        out = ApertureShape::Square;
        return true;
    }
    return false;
}

bool hasColumn(const std::map<std::string, std::size_t>& header, const std::string& name) {
    return header.find(name) != header.end();
}

std::string getCell(const std::vector<std::string>& cells,
                    const std::map<std::string, std::size_t>& header,
                    const std::string& name,
                    const std::string& fallback = "")
{
    auto it = header.find(name);
    if (it == header.end() || it->second >= cells.size()) {
        return fallback;
    }
    if (cells[it->second].empty()) {
        return fallback;
    }
    return cells[it->second];
}

double parseDouble(const std::string& text, const std::string& field) {
    try {
        std::size_t pos = 0;
        double value = std::stod(text, &pos);
        if (pos != text.size()) {
            throw std::invalid_argument("trailing characters");
        }
        return value;
    } catch (const std::exception&) {
        throw std::runtime_error("invalid numeric field '" + field + "': " + text);
    }
}

int parseInt(const std::string& text, const std::string& field) {
    try {
        std::size_t pos = 0;
        int value = std::stoi(text, &pos);
        if (pos != text.size()) {
            throw std::invalid_argument("trailing characters");
        }
        return value;
    } catch (const std::exception&) {
        throw std::runtime_error("invalid integer field '" + field + "': " + text);
    }
}

} // namespace

bool readMirrorFacetsCSV(const std::string& path,
                         std::vector<MirrorFacet>& facets,
                         std::string* error)
{
    facets.clear();

    std::ifstream ifs(path);
    if (!ifs) {
        if (error) {
            *error = "failed to open file: " + path;
        }
        return false;
    }

    std::string line;
    if (!std::getline(ifs, line)) {
        if (error) {
            *error = "empty CSV file";
        }
        return false;
    }

    std::vector<std::string> header_cells = splitCsvLine(line);
    std::map<std::string, std::size_t> header;
    for (std::size_t i = 0; i < header_cells.size(); ++i) {
        header[lowerCopy(header_cells[i])] = i;
    }

    const std::vector<std::string> required = {
        "id", "center_x", "center_y", "center_z",
        "normal_x", "normal_y", "normal_z",
        "surface_type", "radius_of_curvature",
        "aperture_shape", "size1"
    };

    for (const auto& name : required) {
        if (!hasColumn(header, name)) {
            if (error) {
                *error = "missing required column: " + name;
            }
            return false;
        }
    }

    int line_no = 1;
    try {
        while (std::getline(ifs, line)) {
            ++line_no;
            if (trim(line).empty()) {
                continue;
            }

            std::vector<std::string> cells = splitCsvLine(line);
            MirrorFacet facet;
            facet.id = parseInt(getCell(cells, header, "id"), "id");
            facet.center.x = parseDouble(getCell(cells, header, "center_x"), "center_x");
            facet.center.y = parseDouble(getCell(cells, header, "center_y"), "center_y");
            facet.center.z = parseDouble(getCell(cells, header, "center_z"), "center_z");
            facet.normal.x = parseDouble(getCell(cells, header, "normal_x"), "normal_x");
            facet.normal.y = parseDouble(getCell(cells, header, "normal_y"), "normal_y");
            facet.normal.z = parseDouble(getCell(cells, header, "normal_z"), "normal_z");
            facet.normal = facet.normal.normalized();

            if (!parseSurfaceType(getCell(cells, header, "surface_type"), facet.surface_type)) {
                throw std::runtime_error("invalid surface_type: " +
                                         getCell(cells, header, "surface_type"));
            }
            facet.radius_of_curvature = parseDouble(getCell(cells, header, "radius_of_curvature"),
                                                    "radius_of_curvature");

            if (!parseApertureShape(getCell(cells, header, "aperture_shape"), facet.aperture_shape)) {
                throw std::runtime_error("invalid aperture_shape: " +
                                         getCell(cells, header, "aperture_shape"));
            }

            facet.size1 = parseDouble(getCell(cells, header, "size1"), "size1");
            facet.size2 = parseDouble(getCell(cells, header, "size2", "0"), "size2");
            facet.aperture_rotation_rad = parseDouble(getCell(cells, header, "aperture_rotation_rad", "0"),
                                                      "aperture_rotation_rad");
            facet.reflectivity_scale = parseDouble(getCell(cells, header, "reflectivity_scale", "1"),
                                                   "reflectivity_scale");
            facet.roughness_sigma_rad = parseDouble(getCell(cells, header, "roughness_sigma_rad", "0"),
                                                    "roughness_sigma_rad");
            facet.misalign_sigma_rad = parseDouble(getCell(cells, header, "misalign_sigma_rad", "0"),
                                                   "misalign_sigma_rad");

            facets.push_back(facet);
        }
    } catch (const std::exception& ex) {
        if (error) {
            std::ostringstream oss;
            oss << "line " << line_no << ": " << ex.what();
            *error = oss.str();
        }
        facets.clear();
        return false;
    }

    std::string validation_error;
    if (!validateMirrorFacets(facets, &validation_error)) {
        if (error) {
            *error = validation_error;
        }
        facets.clear();
        return false;
    }

    if (error) {
        error->clear();
    }
    return true;
}
