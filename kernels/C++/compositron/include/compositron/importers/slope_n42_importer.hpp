#pragma once

#include <string>

#include "compositron/core/measurement.hpp"


namespace compositron::importers {


compositron::core::Measurement import_SLOPE_n42(const std::string& filename);


} // namespace compositron::importers
