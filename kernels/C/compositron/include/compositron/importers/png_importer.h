#ifndef COMPOSITRON_PNG_IMPORTER_H
#define COMPOSITRON_PNG_IMPORTER_H

#include <stddef.h>
#include <stdint.h>

typedef struct {
    size_t width;
    size_t height;
    uint32_t *data;
} Image;

// Import a PNG file, treating the channels as contiguous bins. The caller is
// responsible for freeing the allocated data in the returned Image struct.
Image import_png(const char *filename);

#endif
