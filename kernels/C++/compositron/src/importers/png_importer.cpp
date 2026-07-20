#include <compositron/importers/png_importer.hpp>

#define STB_IMAGE_IMPLEMENTATION
extern "C" {
    #include "compositron/importers/vendored/stb_image.h"
}

#include <vector>
#include <cstdint>
#include <stdexcept>
#include <format>
#include <utility>

namespace compositron::importers {

Array2D import_png(const std::string& filename) {
    int width, height, channels;

    unsigned char *img = stbi_load(
        filename.c_str(), &width, &height, &channels, 3
    );

    if (img == nullptr) {
        stbi_image_free(img);
        throw std::runtime_error(std::format(
            "Could not load PNG file {}", filename
        ));
    }

    if (!std::in_range<size_t>(height) || !std::in_range<size_t>(width)) {
        stbi_image_free(img);
        throw std::runtime_error(std::format(
            "PNG file {} has invalid format {} x {}", filename, width, height
        ));
    }

    if (!std::in_range<size_t>(3 * width * height)) {
        stbi_image_free(img);
        throw std::runtime_error(std::format(
            "PNG file {} is too large ({} x {})", filename, width, height
        ));
    }

    size_t length = width * height;

    std::vector<uint32_t> data(length);

    for (size_t i = 0; i < length; ++i) {
        size_t index = i * 3;
        data[i] = img[index] | (img[index + 1] << 8) | (img[index + 2] << 16);
    }

    stbi_image_free(img);

    return {
        static_cast<size_t>(width), static_cast<size_t>(height), data
    };
}

} // namespace compositron::importers
