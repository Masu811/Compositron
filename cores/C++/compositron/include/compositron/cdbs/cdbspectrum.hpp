#pragma once

#include "compositron/core/utils.hpp"

#include <cmath>
#include <cstdint>
#include <tuple>
#include <vector>
#include <string>

namespace compositron::cdbs {

typedef std::tuple<
    std::tuple<double, double>, std::tuple<double, double>
> ecal_t;

class CDBSpectrum {
public:
    compositron::core::utils::Spectrum2DPtr spectrum;
    std::string detpair;
    ecal_t ecal;
    double s = NAN;
    double ds = NAN;
    double w = NAN;
    double dw = NAN;
    uint64_t counts = 0;
    double dcounts = NAN;
    double roi_counts = NAN;
    double droi_counts = NAN;

    CDBSpectrum() = delete;

    CDBSpectrum(
        core::utils::Spectrum2DPtr spectrum,
        std::string detname,
        ecal_t ecal = {}
    );

    template<typename U>
    CDBSpectrum(
        std::vector<U>&& spectrum,
        std::tuple<size_t, size_t> shape,
        std::string detname,
        ecal_t ecal = {}
    );
};

} // namespace compositron::cdbs
