#pragma once

#include <vector>
#include <string>
#include <cstdint>

namespace compositron::importers {
    typedef struct {
        size_t width;
        size_t height;
        std::vector<uint32_t> data;
    } Array2D;

    Array2D import_png(const std::string& filename);
}
