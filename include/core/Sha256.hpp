#pragma once

#include <string>

namespace lact {

std::string sha256String(const std::string& value);
std::string sha256File(const std::string& path);

} // namespace lact
