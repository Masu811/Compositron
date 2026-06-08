#include "compositron/importers/png_importer.h"

#include <string.h>
#include <stdint.h>
#include <stddef.h>
#include <stdio.h>

#define STB_IMAGE_IMPLEMENTATION
#include "compositron/importers/vendored/stb_image.h"


Image import_png(const char *filename) {
    Image img = {0};

    size_t filename_size = strlen(filename);

    if (filename_size > 200) {
        fprintf(stderr, "File path too long");
        return img;
    }

    int width, height, channels;
    unsigned char *img_file = stbi_load(filename, &width, &height, &channels, 3);

    if (img_file == NULL) {
        char buf[225] = {0};
        sprintf(buf, "Could not open file %s", filename);
        perror(buf);
        return img;
    }

    img.data = (uint32_t *)malloc(sizeof(uint32_t) * width * height);

    if (img.data == NULL) {
        perror("Failed to allocate memory for coincidence spectrum");
        stbi_image_free(img_file);
        return img;
    }

    img.width = width;
    img.height = height;

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int index = (y * width + x) * 3;
            int value = (
                img_file[index]
                | (img_file[index + 1] << 8)
                | (img_file[index + 2] << 16)
            );
            img.data[y * width + x] = value;
        }
    }

    stbi_image_free(img_file);

    return img;
}
