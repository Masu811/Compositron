#pragma once

#include <cmath>
#include <cstdint>
#include <vector>
#include <string>
#include <map>

#include "compositron/core/utils.hpp"
#include "compositron/core/fitting.hpp"
#include "compositron/cdbs/anti_aliasing.hpp"


namespace compositron::cdbs {


using anti_aliasing::Polygon;


struct DiagonalArea {
    core::utils::Unit width_cel;
    core::utils::Unit width_cml;
    core::utils::Unit offset_cel;
    core::utils::Unit offset_cml;
};

struct AxisAlignedArea {
    std::array<core::utils::Unit, 2> first_det_bnds;
    std::array<core::utils::Unit, 2> second_det_bnds;
};

struct EllipseArea {
    core::utils::Unit radius_cel;
    core::utils::Unit radius_cml;
};

struct LineshapeParam {
    double val;
    double err;
    std::vector<Polygon> num;
    std::vector<Polygon> denom;
};


class CDBSpectrum {
public:
    compositron::core::utils::Spectrum2D spectrum;
    compositron::core::utils::EnergyDetectorPair detpair;
    uint64_t counts = 0;
    double dcounts = NAN;
    std::optional<Eigen::MatrixXd> peak;
    std::optional<std::array<std::array<size_t, 2>, 2>> peak_bnds;
    double peak_counts;
    double dpeak_counts;
    std::map<std::string, core::fitting::SimpleFitParam> peak_params;
    std::map<std::string, LineshapeParam> lineshape_params;

    CDBSpectrum() = delete;

    CDBSpectrum(CDBSpectrum& other) = delete;

    CDBSpectrum(CDBSpectrum&& other) = default;

    CDBSpectrum(
        core::utils::Spectrum2D&& spectrum,
        core::utils::EnergyDetectorPair detpair
    ) : spectrum(spectrum), detpair(detpair) {}

    template<typename U>
    CDBSpectrum(
        const std::vector<U>& spectrum,
        size_t nrows,
        size_t ncols,
        core::utils::StorageOrder order,
        core::utils::EnergyDetectorPair detpair
    ) : CDBSpectrum(
            core::utils::Spectrum2D(spectrum, nrows, ncols, order),
            detpair
    ) {}
};


} // namespace compositron::cdbs
