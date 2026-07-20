#include "compositron/core/utils.h"

#include <stdint.h>
#include <string.h>
#include <ctype.h>
#include <stdio.h>

Vector *vector_new(size_t init_capacity, size_t element_size) {
    Vector *vec = (Vector*)calloc(1, sizeof(Vector));

    if (vec == NULL) return NULL;

    if (init_capacity == 0) init_capacity = 25;

    void *data = calloc(init_capacity, element_size);

    if (data == NULL) {
        free(vec);
        return NULL;
    }

    vec->data = data;
    vec->element_size = element_size;
    vec->capacity = init_capacity;

    return vec;
}


Vector *vector_new_with_deallocator(
    size_t init_capacity, size_t element_size, void (*freeItem)(void*)
) {
   Vector *vec = vector_new(init_capacity, element_size);
   if (vec == NULL) return NULL;
   vec->freeItem = freeItem;
   return vec;
}


static bool vector_grow(Vector *vec) {
    if (vec->capacity > SIZE_MAX / (2 * vec->element_size)) return false;

    void *new_data = realloc(
        vec->data, 2 * vec->capacity * vec->element_size
    );

    if (new_data == NULL) return false;

    vec->data = new_data;
    vec->capacity *= 2;

    return true;
}


bool vector_push(Vector *vec, void *item) {
    if (vec->size == vec->capacity) {
        if (!vector_grow(vec)) return false;
    }

    memcpy(
        (unsigned char*)vec->data + vec->size++ * vec->element_size,
        item,
        vec->element_size
    );

    return true;
}


void *vector_get(Vector *vec, size_t idx) {
    if (idx >= vec->size) return NULL;
    return (unsigned char*)vec->data + vec->element_size * idx;
}


void vector_set(Vector *vec, size_t idx, void *item) {
    if (idx >= vec->size) return;

    memcpy(
        (unsigned char*)vec->data + idx * vec->element_size,
        item,
        vec->element_size
    );
}


void vector_free(Vector *vec) {
    if (vec == NULL) return;

    if (vec->data != NULL) {
        if (vec->freeItem != NULL) {
            for (size_t i = 0; i < vec->size; ++i) {
                vec->freeItem(vector_get(vec, i));
            }
        }

        free(vec->data);
    }

    free(vec);
}


// Hashmap


// Taken from
// Stroustrup, Bjarne (1997). The C++ Programming Language Third Edition.
// Reading Massachusetts: Addison-Wesley. p. 503. ISBN 0-201-88954-4.
static size_t hashmap_hash(const char *key) {
    size_t res = 0;

    size_t len = strlen(key);
    const unsigned char *p = (const unsigned char*)key;

    while (len--) res = (res << 1) ^ *p++;

    return res;
}


Hashmap *hashmap_new(size_t capacity) {
    if (capacity == 0) capacity = 25;

    Hashmap *map = (Hashmap*)calloc(1, sizeof(Hashmap));

    if (map == NULL) return NULL;

    void *table = calloc(capacity, sizeof(Vector*));

    if (table == NULL) {
        free(map);
        return NULL;
    }

    map->table = table;
    map->capacity = capacity;

    Vector *keys = vector_new(0, sizeof(char*));

    if (keys == NULL) {
        free(map);
        free(table);
        return NULL;
    }

    map->keys = keys;

    return map;
}


Hashmap *hashmap_new_with_deallocator(
    size_t capacity, void (*freeValue)(void*)
) {
    Hashmap *map = hashmap_new(capacity);
    if (map == NULL) return NULL;
    map->freeValue = freeValue;
    return map;
}


static Vector *create_bucket(Hashmap *map, size_t idx) {
    if (map->table[idx] == NULL) {
        map->table[idx] = vector_new(0, sizeof(HashmapItem));
        if (map->table[idx] == NULL) return NULL;
    }

    return map->table[idx];
}


void *hashmap_get(Hashmap *map, const char *key) {
    size_t idx = hashmap_hash(key) % map->capacity;
    Vector *bucket = map->table[idx];
    if (bucket == NULL) return NULL;

    for (size_t i = 0; i < bucket->size; ++i) {
        HashmapItem *item = vector_get(bucket, i);
        if (strcmp(item->key, key) == 0) {
            return item->value;
        }
    }

    return NULL;
}


HashmapItem hashmap_iter(Hashmap *map, size_t idx) {
    if (idx >= map->size) return (HashmapItem){0};

    char **key = NULL;

    do {
        key = vector_get(map->keys, idx++);
    } while (*key == NULL && idx++ < map->keys->size);

    if (*key == NULL) return (HashmapItem){0};

    void *value = hashmap_get(map, *key);

    return (HashmapItem){ .key = *key, .value = value, .idx = idx };
}


static bool hashmap_insert_without_grow(Hashmap *map, char *key, void *value) {
    size_t idx = hashmap_hash(key) % map->capacity;
    Vector *bucket = map->table[idx];
    if (bucket == NULL) bucket = create_bucket(map, idx);
    if (bucket == NULL) return false;

    // Check if key is aleady in map

    for (size_t i = 0; i < bucket->size; ++i) {
        HashmapItem *item = vector_get(bucket, i);
        if (strcmp(item->key, key) == 0) {
            if (map->freeValue) map->freeValue(item->value);
            item->value = value;
            // Index remains the same
            return true;
        }
    }

    // If not, create new entry

    bool inserted = vector_push(map->keys, &key);

    if (!inserted) return false;

    inserted = vector_push(bucket, &(HashmapItem){
        .key = key, .value = value, .idx = map->keys->size-1
    });

    if (!inserted) {
        char **key = vector_get(map->keys, map->size--);
        if (key != NULL) free(key);
    }

    map->size++;

    return true;
}


static bool hashmap_grow(Hashmap *map) {
    if (map->capacity > SIZE_MAX / 2) return false;

    Hashmap *new = hashmap_new(2 * map->capacity);

    if (new == NULL) return false;

    // Insert items into new table (insertion order remains the same)

    for (size_t i = 0; i < map->size; ++i) {
        HashmapItem item = hashmap_iter(map, i);
        hashmap_insert_without_grow(new, item.key, item.value);
    }

    // Clean up the old table

    for (size_t i = 0; i < map->capacity; ++i) {
        Vector *bucket = map->table[i];
        if (bucket != NULL) vector_free(bucket);
    }

    // This does not free the keys themselves since we still need them
    vector_free(map->keys);
    free(map->table);

    // Overwrite map with new_map

    map->table = new->table;
    map->keys = new->keys;
    map->capacity = new->capacity;

    free(new);

    return true;
}


bool hashmap_insert(Hashmap *map, char *key, void *value) {
    if (map == NULL) return false;

    if (map->size > 0.75 * map->capacity) {
        bool grown = hashmap_grow(map);
        if (!grown) return false;
    }

    return hashmap_insert_without_grow(map, key, value);
}


void hashmap_free(Hashmap *map) {
    for (size_t i = 0; i < map->capacity; ++i) {
        Vector *bucket = map->table[i];

        if (bucket == NULL) continue;

        for (size_t i = 0; i < bucket->size; ++i) {
            HashmapItem *item = vector_get(bucket, i);
            if (map->freeValue != NULL) map->freeValue(item->value);
            free(item->key);
        }

        vector_free(bucket);
    }

    free(map->table);
    vector_free(map->keys);
    free(map);
}


void hashmap_print(Hashmap *map) {
    for (size_t i = 0; i < map->capacity; i++) {
        printf("[%zu]:", i);
        Vector *bucket = map->table[i];

        if (bucket == NULL) {
            printf(" ---\n");
            continue;
        }

        for (size_t j = 0; j < bucket->size; j++) {
            HashmapItem *item = vector_get(bucket, j);
            printf(" [%zu]: {key: %s}", j, item->key);
        }
        printf("\n");
    }
}


static uint64_t generic_max(void *arr, SpectrumType type, size_t size) {
    uint64_t max;

    switch (type) {
        case U8:
            max = *(uint8_t*)arr;
            for (size_t i = 1; i < size; ++i) {
                if (((uint8_t*)arr)[i] > max) max = ((uint8_t*)arr)[i];
            }
            break;
        case U16:
            max = *(uint16_t*)arr;
            for (size_t i = 1; i < size; ++i) {
                if (((uint16_t*)arr)[i] > max) max = ((uint16_t*)arr)[i];
            }
            break;
        case U32:
            max = *(uint32_t*)arr;
            for (size_t i = 1; i < size; ++i) {
                if (((uint32_t*)arr)[i] > max) max = ((uint32_t*)arr)[i];
            }
            break;
        case U64:
            max = *(uint64_t*)arr;
            for (size_t i = 1; i < size; ++i) {
                if (((uint64_t*)arr)[i] > max) max = ((uint64_t*)arr)[i];
            }
            break;
    }

    return max;
}


#define CONVERT(src_t, dst_t) \
for (size_t i = 0; i < size; ++i) { \
    ((dst_t*)dst)[i] = ((src_t*)src)[i]; \
}


static void convert_arrays(
    void *src, SpectrumType src_t, void* dst, SpectrumType dst_t, size_t size
) {
    switch (src_t) {
        case U8:
            switch (dst_t) {
                case U8: CONVERT(uint8_t, uint8_t); break;
                case U16: CONVERT(uint8_t, uint16_t); break;
                case U32: CONVERT(uint8_t, uint32_t); break;
                case U64: CONVERT(uint8_t, uint64_t); break;
            }
            break;
        case U16:
            switch (dst_t) {
                case U8: CONVERT(uint16_t, uint8_t); break;
                case U16: CONVERT(uint16_t, uint16_t); break;
                case U32: CONVERT(uint16_t, uint32_t); break;
                case U64: CONVERT(uint16_t, uint64_t); break;
            }
            break;
        case U32:
            switch (dst_t) {
                case U8: CONVERT(uint32_t, uint8_t); break;
                case U16: CONVERT(uint32_t, uint16_t); break;
                case U32: CONVERT(uint32_t, uint32_t); break;
                case U64: CONVERT(uint32_t, uint64_t); break;
            }
            break;
        case U64:
            switch (dst_t) {
                case U8: CONVERT(uint64_t, uint8_t); break;
                case U16: CONVERT(uint64_t, uint16_t); break;
                case U32: CONVERT(uint64_t, uint32_t); break;
                case U64: CONVERT(uint64_t, uint64_t); break;
            }
            break;
    }
}


void toMinTypeSpectrum(Spectrum *spectrum) {
    if (spectrum == NULL) return;

    void *src_data = spectrum->data;
    SpectrumType src_type = spectrum->type;
    size_t size = spectrum->size;

    if (src_data == NULL || size == 0) return;

    size_t src_tsize;

    switch (src_type) {
        case U8:
            return;
        case U16:
            src_tsize = sizeof(uint16_t);
            break;
        case U32:
            src_tsize = sizeof(uint32_t);
            break;
        case U64:
            src_tsize = sizeof(uint64_t);
            break;
    }

    uint64_t max = generic_max(src_data, src_type, size);

    size_t dst_tsize;
    SpectrumType dst_type;

    if (max <= UINT8_MAX) {
        if (src_type == U8) return;
        dst_tsize = sizeof(uint8_t);
        dst_type = U8;
    } else if (max <= UINT16_MAX) {
        if (src_type == U16) return;
        dst_tsize = sizeof(uint16_t);
        dst_type = U16;
    } else if (max <= UINT32_MAX) {
        if (src_type == U32) return;
        dst_tsize = sizeof(uint32_t);
        dst_type = U32;
    } else {
        return;
    }

    // We can convert in-place, because dst is smaller than src
    convert_arrays(src_data, src_type, src_data, dst_type, size);

    spectrum->type = dst_type;

    // We still need realloc to shrink the mem-usage. realloc might shrink the
    // memory in place, but it's not guaranteed. Therefore new_data
    // is treated like a new pointer.
    void *new_data = realloc(src_data, size * dst_tsize);

    if (new_data == NULL) return;

    spectrum->data = new_data;
}

void toMinTypeSpectrum2D(Spectrum2D *spectrum) {
    if (
        spectrum == NULL || spectrum->data == NULL ||
        spectrum->width == 0 || spectrum->height == 0
    ) return;

    if (spectrum->height > SIZE_MAX / spectrum->width) return;

    Spectrum new_spectrum = {
        .data = spectrum->data,
        .type = spectrum-> type,
        .size = spectrum->width * spectrum->height
    };

    toMinTypeSpectrum(&new_spectrum);

    spectrum->data = new_spectrum.data;
    spectrum->type = new_spectrum.type;
}

static char *ltrim(char *s) {
    if (!s) return NULL;

    char *p = s;

    while (isspace(*p)) p++;

    if (p > s) {
        char *new = strdup(p);
        if (new == NULL) return s;
        free(s);
        return new;
    }

    return s;
}

static char *rtrim(char *s) {
    if (!s) return NULL;

    size_t size = strlen(s);

    if (size == 0) return s;

    char* back = s + size - 1;

    while (back > s && isspace(*back)) back--;

    *(back+1) = '\0';

    return s;
}

char *trim(char *s) {
    if (!s) return NULL;
    return ltrim(rtrim(s));
}
