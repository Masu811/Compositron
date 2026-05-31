#pragma once

#include "compositron/core/measurement.hpp"

#include <string>

namespace compositron::importers {

compositron::core::Measurement import_SLOPE_n42(const std::string& filename);

} // namespace compositron::importers
