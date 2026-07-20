#include "compositron/cdbs/cdbspectrum.h"

#include <stdio.h>
#include <math.h>

CDBSpectrum *newCDBSpectrum() {
    CDBSpectrum *c = (CDBSpectrum*)calloc(1, sizeof(CDBSpectrum));

    if (c == NULL) {
        fprintf(stderr, "Could not allocate memory for CDB Spectrum\n");
        return NULL;
    }

    c->s = NAN;
    c->ds = NAN;
    c->w = NAN;
    c->dw = NAN;
    c->dcounts = NAN;
    c->roi_counts = NAN;
    c->droi_counts = NAN;

    return c;
}


void freeCDBSpectrum(void *spectrum) {
    if (spectrum == NULL) return;

    CDBSpectrum *c = spectrum;

    if (c->spectrum.data != NULL) free(c->spectrum.data);
    if (c->detpair != NULL) free(c->detpair);

    free(c);
}
