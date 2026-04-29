#pragma once
#include "core/PhotonBunch.hpp"

class PhotonSource {
public:
    virtual ~PhotonSource() = default;

    // Return true if one bunch is produced, false if no more data.
    virtual bool next(PhotonBunch& out) = 0;

    // Reset to the beginning if supported.
    virtual void reset() = 0;
};
