#ifndef COMPOSITRON_CDBSPECTRUM_H
#define COMPOSITRON_CDBSPECTRUM_H

#include "compositron/core/utils.h"
#include "compositron/dbs/dbspectrum.h"
#include <stddef.h>
#include <stdint.h>

typedef struct {
    d_ecal_t ecal_1;
    d_ecal_t ecal_2;
} c_ecal_t;

typedef struct {
    Spectrum2D spectrum;
    char* detpair;
    c_ecal_t ecal;
    double s;
    double ds;
    double w;
    double dw;
    uint64_t counts;
    double dcounts;
    double roi_counts;
    double droi_counts;
} CDBSpectrum;

CDBSpectrum *newCDBSpectrum();

void freeCDBSpectrum(void *c);

#endif
