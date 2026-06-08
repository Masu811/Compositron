#ifndef COMPOSITRON_UTILS_H
#define COMPOSITRON_UTILS_H

#include <stdlib.h>
#include <stdbool.h>

typedef struct {
    void *data;
    size_t element_size;
    size_t capacity;
    size_t size;
    void (*freeItem)(void*);
} Vector;

Vector *vector_new(size_t init_capacity, size_t element_size);

Vector *vector_new_with_deallocator(
    size_t init_capacity, size_t element_size, void (*freeItem)(void*)
);

bool vector_push(Vector *vec, void *item);

void *vector_get(Vector *vec, size_t idx);

void vector_set(Vector *vec, size_t idx, void *item);

void vector_free(Vector *vec);


typedef struct {
    Vector **table;
    size_t capacity;
    size_t size;
    Vector *keys;
    void (*freeValue)(void*);
} Hashmap;


typedef struct {
    char *key;
    void *value;
    size_t idx;
} HashmapItem;


Hashmap *hashmap_new(size_t capacity);

Hashmap *hashmap_new_with_deallocator(
    size_t capacity, void (*freeValue)(void*)
);

void *hashmap_get(Hashmap *map, const char *key);

HashmapItem hashmap_iter(Hashmap *map, size_t idx);

bool hashmap_insert(Hashmap *map, char *key, void *value);

void hashmap_free(Hashmap *map);

typedef enum {
    U8,
    U16,
    U32,
    U64,
} SpectrumType;

typedef struct {
    void *data;
    SpectrumType type;
    size_t size;
} Spectrum;

typedef struct {
    void *data;
    SpectrumType type;
    size_t width;
    size_t height;
} Spectrum2D;

// Convert data in spectrum to the smallest type to hold the data.
// spectrum.data may or may not point to the same data as before the operation.
// If not, the dropped data is freed automatically.
void toMinTypeSpectrum(Spectrum *spectrum);

// Convert data in spectrum to the smallest type to hold the data.
// spectrum.data may or may not point to the same data as before the operation.
// If not, the dropped data is freed automatically.
void toMinTypeSpectrum2D(Spectrum2D *spectrum);

char *trim(char *s);

#endif
