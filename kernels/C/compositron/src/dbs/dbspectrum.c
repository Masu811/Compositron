#include "compositron/dbs/dbspectrum.h"

#include <stdio.h>
#include <math.h>

DBSpectrum *newDBSpectrum() {
    DBSpectrum *d = (DBSpectrum*)calloc(1, sizeof(DBSpectrum));

    if (d == NULL) {
        fprintf(stderr, "Could not allocate memory for DB Spectrum\n");
        return NULL;
    }

    d->s = NAN;
    d->ds = NAN;
    d->w = NAN;
    d->dw = NAN;
    d->v2p = NAN;
    d->dv2p = NAN;
    d->dcounts = NAN;
    d->peak_counts = NAN;
    d->dpeak_counts = NAN;

    return d;
}


void freeDBSpectrum(void *spectrum) {
    if (spectrum == NULL) return;

    DBSpectrum *d = (DBSpectrum*)spectrum;

    if (d->spectrum.data != NULL) free(d->spectrum.data);
    if (d->detname != NULL) free(d->detname);

    free(d);
}
