#pragma once

#include "compositron/core/utils.hpp"
#include <cmath>
#include <cstdint>
#include <tuple>
#include <vector>
#include <string>

namespace compositron::dbs {

typedef std::tuple<double, double> ecal_t;

class DBSpectrum {
public:
    core::utils::SpectrumPtr spectrum;
    std::string detname;
    ecal_t ecal;
    double s = NAN;
    double ds = NAN;
    double w = NAN;
    double dw = NAN;
    double v2p = NAN;
    double dv2p = NAN;
    uint64_t counts = 0;
    double dcounts = NAN;
    double peak_counts = NAN;
    double dpeak_counts = NAN;

    DBSpectrum() = delete;

    DBSpectrum(
        core::utils::SpectrumPtr spectrum,
        std::string detname,
        ecal_t ecal = {}
    );

    template<typename T>
    DBSpectrum(
        std::vector<T>&& spectrum, std::string detname, ecal_t ecal = {}
    );
};

} // namespace compositron::dbs
