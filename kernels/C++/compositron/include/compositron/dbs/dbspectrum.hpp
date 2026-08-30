#pragma once

#include <cmath>
#include <cstdint>
#include <variant>
#include <vector>
#include <array>

#include "compositron/core/utils.hpp"
#include "compositron/core/fitting.hpp"


namespace compositron::dbs {


enum PeakModel {
    GAUSS,
    ERFLINEAR1GAUSS,
    ERFLINEAR2GAUSS,
    ERFLINEAR3GAUSS,
};


class PeakArea {
public:
    core::utils::Unit width;

    std::array<double, 2> to_bnds_kev(double eres) const;
};


class PeakAreaWithOffset {
public:
    core::utils::Unit width;
    core::utils::Unit offset;

    std::array<double, 2> to_bnds_kev(double eres) const;
};


class SpectrumArea {
public:
    double left_bnd_kev;
    double right_bnd_kev;

    std::array<double, 2> to_bnds_kev(double eres) const;
};


class Area {
private:
    std::variant<PeakArea, PeakAreaWithOffset, SpectrumArea> data;

public:
    Area() = delete;

    Area(core::utils::Unit width);
    Area(core::utils::Unit width, core::utils::Unit offset);
    Area(double left_bnd_kev, double right_bnd_kev);

    std::array<double, 2> to_bnds_kev(double eres) const;
};


enum BinIntegrationScheme {
    CONST,
    LINEAR,
};


struct LineshapeParam {
    double val;
    double err;
    std::vector<Area> num;
    std::vector<Area> denom;
    bool in_corrected_peak;
};


struct LineshapeParamDefinition {
    std::string name;
    std::vector<Area> num;
    std::vector<Area> denom;
    bool in_corrected_peak;
};


const std::array<LineshapeParamDefinition, 4> STD_LINESHAPE_PARAMS{{
    {
        "S",
        {Area(core::utils::Unit(2., core::utils::KEV))},
        {Area(core::utils::Unit(20., core::utils::KEV))},
        true,
    },
    {
        "W",
        {
            Area(
                core::utils::Unit(1., core::utils::KEV),
                core::utils::Unit(3., core::utils::KEV)
            ),
            Area(
                core::utils::Unit(1., core::utils::KEV),
                core::utils::Unit(-3., core::utils::KEV)
            )
        },
        {Area(core::utils::Unit(20., core::utils::KEV))},
        true,
    },
    {
        "V/P", {Area(400, 500)}, {Area(501, 521)}, false,
    },
    {
        "P/T",
        {Area(501, 521)},
        {Area(-INFINITY, INFINITY)},
        false,
    },
}};


class DBSpectrum {
public:
    core::utils::Spectrum spectrum;
    core::utils::EnergyDetector detector;
    uint64_t counts;
    std::optional<Eigen::VectorXd> peak;
    std::optional<std::array<size_t, 2>> peak_bnd_idcs;
    std::optional<double> peak_counts;
    std::unordered_map<std::string, core::fitting::SimpleFitParam> peak_params;
    std::unordered_map<std::string, LineshapeParam> lineshape_params;

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

private:
    Eigen::VectorXd get_peak_energies();

    void fit_peak(PeakModel peak_model);

    void subtract_bg(PeakModel peak_model);

public:
    void extract_peak(double peak_width, PeakModel peak_model);

    void correct_ecal(core::utils::EcalCorrectionOrder order);

    double integrate(
        Area area,
        bool in_corrected_peak,
        BinIntegrationScheme bin_integration_scheme
    );

    LineshapeParam& calc_lineshape_param(
        LineshapeParamDefinition definition,
        BinIntegrationScheme bin_integration_scheme
    );

    void default_analyze();
};


} // namespace compositron::dbs
