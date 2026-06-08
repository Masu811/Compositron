#include "compositron/core/measurement.h"

#include "compositron/core/utils.h"
#include "compositron/dbs/dbspectrum.h"
#include "compositron/cdbs/cdbspectrum.h"

#include <stdlib.h>
#include <stdio.h>


Measurement *newMeasurement() {
    Measurement *m = (Measurement*)calloc(1, sizeof(Measurement));

    m->dbs = hashmap_new_with_deallocator(0, freeDBSpectrum);
    if (m->dbs == NULL) goto failure;

    m->cdbs = hashmap_new_with_deallocator(0, freeCDBSpectrum);
    if (m->cdbs == NULL) goto failure;

    m->metadata = hashmap_new_with_deallocator(0, free);
    if (m->metadata == NULL) goto failure;

    return m;

failure:

    freeMeasurement(m);

    return NULL;
}


void freeMeasurement(Measurement *m) {
    if (m == NULL) return;

    if (m->filename != NULL) free(m->filename);
    if (m->name != NULL) free(m->name);
    if (m->dbs != NULL) hashmap_free(m->dbs);
    if (m->cdbs != NULL) hashmap_free(m->cdbs);
    if (m->metadata != NULL) hashmap_free(m->metadata);

    free(m);
}


MeasurementShape getMeasurementShape(Measurement *m) {
    MeasurementShape shape;

    shape.dbs = m->dbs->size;
    shape.cdbs = m->cdbs->size;

    return shape;
}
