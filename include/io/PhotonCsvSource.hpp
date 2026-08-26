#pragma once
#include <cstddef>
#include <fstream>
#include <map>
#include <optional>
#include <string>
#include "io/PhotonSource.hpp"

struct PhotonCsvConfig {
    std::string csv_path;
    bool local_telescope_frame = true;
    double default_wavelength_nm = 400.0;
    double default_time_ns = 0.0;
    double default_weight = 1.0;
    double default_multiplicity = 1.0;
    int default_event_id = 0;
    int default_telescope_id = 0;
    bool default_eventio_2d = false;
    bool filter_telescope_id = false;
    int selected_telescope_id = 0;
    bool filter_event_id = false;
    int selected_event_id = 0;
};

class PhotonCsvSource : public PhotonSource {
public:
    explicit PhotonCsvSource(const PhotonCsvConfig& cfg);

    bool next(PhotonBunch& out) override;
    void reset() override;

    bool localTelescopeFrame() const { return cfg_.local_telescope_frame; }

private:
    PhotonCsvConfig cfg_;
    std::ifstream input_;
    std::map<std::string, int> header_;
    std::optional<PhotonBunch> buffered_;
    int line_no_ = 0;

    void openAndPrime();
    bool readNext(PhotonBunch& out);
};
