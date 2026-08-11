#pragma once

#include <vector>
#include <string>
#include <cstdint>


namespace compositron::importers {


struct Array2D {
    size_t width;
    size_t height;
    std::vector<uint32_t> data;
};

Array2D import_png(const std::string& filename);


} // namespace compositron::importers
