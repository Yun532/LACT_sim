#pragma once
#include <cstddef>
#include <string>
#include <vector>
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

    std::size_t size() const { return rows_.size(); }
    bool localTelescopeFrame() const { return cfg_.local_telescope_frame; }

private:
    PhotonCsvConfig cfg_;
    std::size_t index_ = 0;
    std::vector<PhotonBunch> rows_;

    void load();
};
