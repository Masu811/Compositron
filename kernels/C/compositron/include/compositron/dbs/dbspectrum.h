#ifndef COMPOSITRON_DBSPECTRUM_H
#define COMPOSITRON_DBSPECTRUM_H

#include "compositron/core/utils.h"

#include <stddef.h>
#include <stdint.h>

typedef struct {
    double c1;
    double c2;
} d_ecal_t;

typedef struct {
    Spectrum spectrum;
    char* detname;
    d_ecal_t ecal;
    double s;
    double ds;
    double w;
    double dw;
    double v2p;
    double dv2p;
    uint64_t counts;
    double dcounts;
    double peak_counts;
    double dpeak_counts;
} DBSpectrum;

DBSpectrum *newDBSpectrum();

void freeDBSpectrum(void *d);

#endif
