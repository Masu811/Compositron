#pragma once

#include <string>

#include "compositron/core/measurement.hpp"


namespace compositron::importers {


compositron::core::Measurement import_MePS_dat(const std::string& filename);


} // namespace compositron::importers
