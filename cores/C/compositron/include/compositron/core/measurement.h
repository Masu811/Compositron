#ifndef COMPOSITRON_MEASUREMENT_H
#define COMPOSITRON_MEASUREMENT_H

#include "compositron/core/utils.h"

#include <stddef.h>

typedef struct {
    size_t dbs;
    size_t cdbs;
} MeasurementShape;

typedef struct {
    char *filename;
    char *name;
    Hashmap *dbs;
    Hashmap *cdbs;
    Hashmap *metadata;
} Measurement;

Measurement *newMeasurement();

void freeMeasurement(Measurement *m);

MeasurementShape getMeasurementShape(Measurement *m);

#endif
