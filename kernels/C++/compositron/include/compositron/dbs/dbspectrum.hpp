#pragma once

#include <cstdint>
#include <vector>
#include <array>
#include <map>

#include "compositron/core/utils.hpp"
#include "compositron/core/fitting.hpp"


namespace compositron::dbs {


struct PeakArea {
    core::utils::Unit width;
};

struct PeakAreaWithOffset {
    core::utils::Unit width;
    core::utils::Unit offset;
};

struct SpectrumArea {
    double left_bnd;
    double right_bnd;
};

struct LineshapeParam {
    double val;
    double err;
    std::vector<SpectrumArea> num;
    std::vector<SpectrumArea> denom;
};


class DBSpectrum {
public:
    core::utils::Spectrum spectrum;
    core::utils::EnergyDetector detector;
    uint64_t counts;
    double dcounts;
    std::optional<Eigen::VectorXd> peak;
    std::optional<std::array<size_t, 2>> peak_bnds;
    double peak_counts;
    double dpeak_counts;
    std::map<std::string, core::fitting::SimpleFitParam> peak_params;
    std::map<std::string, LineshapeParam> lineshape_params;

    DBSpectrum() = delete;

    DBSpectrum(DBSpectrum& other) = delete;

    DBSpectrum(DBSpectrum&& other) = default;

    DBSpectrum(
        core::utils::Spectrum&& spectrum,
        core::utils::EnergyDetector detector
    );

    template<typename T>
    DBSpectrum(
        const std::vector<T>& spectrum, core::utils::EnergyDetector detector
    ) : DBSpectrum(core::utils::Spectrum(spectrum), detector) {}
};


} // namespace compositron::dbs
