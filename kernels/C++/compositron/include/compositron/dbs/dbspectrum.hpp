#pragma once

#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>
#include <array>
#include <map>

#include "compositron/core/utils.hpp"
#include "compositron/core/fitting.hpp"


namespace compositron::dbs {


enum PeakModel {
    GAUSS,
    ERFLINEAR1GAUSS,
    ERFLINEAR2GAUSS,
    ERFLINEAR3GAUSS,
};


class Area {
private:
    enum AreaType {
        PEAK,
        PEAKOFFSET,
        SPECTRUM
    };

    AreaType type;

    double left_bnd;
    double right_bnd;
    core::utils::Unit width;
    core::utils::Unit offset;

public:
    Area() = delete;

    Area(core::utils::Unit width);
    Area(core::utils::Unit width, core::utils::Unit offset);
    Area(double left_bnd, double right_bnd);

    std::array<double, 2> to_bnds_kev(double eres);
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
        {Area(core::utils::Unit(1., core::utils::KEV))},
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

private:
    Eigen::VectorXd get_peak_energies();

    bool fit_peak(PeakModel peak_model);

    void subtract_bg(PeakModel peak_model);

public:
    void extract_peak(double peak_width, bool bg_corr, PeakModel peak_model);

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
