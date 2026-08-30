#include "compositron/dbs/dbspectrum.hpp"

#include <iostream>

#include <cmath>
#include <cstddef>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <variant>
#include <ranges>

#include <eigen3/Eigen/Core>

#include "compositron/constants.hpp"
#include "compositron/core/utils.hpp"
#include "compositron/dbs/fitting.hpp"


namespace compositron::dbs {


namespace {


template <typename T>
double integrate_const(
    Eigen::VectorX<T> spectrum,
    double lower_ch_bnd,
    double upper_ch_bnd
) {
    double counts = 0;

    size_t lowest_ch = std::round(lower_ch_bnd);
    size_t highest_ch = std::round(upper_ch_bnd);

    counts += spectrum(
        Eigen::seq(lowest_ch, highest_ch)
    ).template cast<double>().sum();

    // Left edge counts

    double lowest_ch_counts = spectrum[lowest_ch];
    double x1 = lowest_ch - 0.5;
    counts -= lowest_ch_counts * (lower_ch_bnd - x1);

    // Right edge counts

    double highest_ch_counts = spectrum[highest_ch];
    double x2 = highest_ch + 0.5;
    counts -= highest_ch_counts * (x2 - upper_ch_bnd);

    return counts;
}


template <typename T>
double integrate_linear(
    Eigen::VectorX<T> spectrum,
    double lower_ch_bnd,
    double upper_ch_bnd
) {
    double counts = 0;

    // index of lowest channel still inside ROI
    size_t lowest_ch = std::round(lower_ch_bnd);
    // index of highest channel still inside ROI
    size_t highest_ch = std::round(upper_ch_bnd);

    counts += spectrum(
        Eigen::seq(lowest_ch + 1, highest_ch - 1)
    ).template cast<double>().sum();

    counts += 0.5 * spectrum[lowest_ch];
    counts += 0.5 * spectrum[highest_ch];

    // Left edge counts

    double x1 = lowest_ch;

    double y1 = spectrum[lowest_ch];
    double y2 = spectrum[lowest_ch + 1];

    double y_l = (y2 - y1) * (lower_ch_bnd - x1) + y1;

    counts -= 0.5 * (y1 + y_l) * (lower_ch_bnd - x1);

    // Right edge counts

    x1 = highest_ch - 1;
    double x2 = highest_ch;

    y1 = spectrum[highest_ch - 1];
    y2 = spectrum[highest_ch];

    double y_u = (y2 - y1) * (upper_ch_bnd - x1) + y1;

    counts -= 0.5 * (y_u + y2) * (x2 - upper_ch_bnd);

    return counts;
}


template <typename T>
double integrate_generic(
    Eigen::VectorX<T> spectrum,
    double lower_ch_bnd,
    double upper_ch_bnd,
    BinIntegrationScheme bin_integration_scheme
) {
    if (std::isinf(lower_ch_bnd)) {
        lower_ch_bnd = 0;
    }
    if (std::isinf(upper_ch_bnd)) {
        upper_ch_bnd = spectrum.size() - 1;
    }

    if (
        lower_ch_bnd < 0 || lower_ch_bnd > spectrum.size() - 1 ||
        upper_ch_bnd < 0 || upper_ch_bnd > spectrum.size() - 1
    ) {
        throw std::runtime_error(
            "Attempting to integrate an area larger than the spectrum or "
            "extracted peak"
        );
    }

    switch (bin_integration_scheme) {
        case BinIntegrationScheme::CONST: {
            return integrate_const(spectrum, lower_ch_bnd, upper_ch_bnd);
        }
        case BinIntegrationScheme::LINEAR: {
            return integrate_linear(spectrum, lower_ch_bnd, upper_ch_bnd);
        }
        default:
            return 0;
    }

}


std::vector<std::array<double, 2>> get_area_bnds(
    std::vector<Area>& areas, double eres
) {
    auto bnds = areas
        | std::views::transform([eres](auto x) {
            return x.to_bnds_kev(eres);
        })
        | std::ranges::to<std::vector<std::array<double, 2>>>();

    for (size_t i = 0; i < bnds.size() - 1; ++i) {
        for (size_t j = i + 1; j < bnds.size(); ++j) {
            if (!(bnds[i][1] < bnds[j][0] || bnds[j][1] < bnds[i][0])) {
                throw std::runtime_error(
                    "Lineshape numerator or denominator areas overlap"
                );
            }
        }
    }

    return bnds;
}


} // anonymous namespace


// class Area

std::array<double, 2> PeakArea::to_bnds_kev(double eres) const {
    double width_kev = width.to_kev(eres);

    double left_bnd_kev = M_E_KEV - width_kev / 2;
    double right_bnd_kev = M_E_KEV + width_kev / 2;

    return { left_bnd_kev, right_bnd_kev };
}


std::array<double, 2> PeakAreaWithOffset::to_bnds_kev(double eres) const {
    double width_kev = width.to_kev(eres);
    double offset_kev = offset.to_kev(eres);

    double left_bnd_kev = M_E_KEV - width_kev / 2 + offset_kev;
    double right_bnd_kev = M_E_KEV + width_kev / 2 + offset_kev;

    return { left_bnd_kev, right_bnd_kev };
}


std::array<double, 2> SpectrumArea::to_bnds_kev(double eres) const {
    return { left_bnd_kev, right_bnd_kev };
}


Area::Area(core::utils::Unit width) : data(PeakArea(width)) {}


Area::Area(
    core::utils::Unit width, core::utils::Unit offset
) : data(PeakAreaWithOffset(width, offset)) {}


Area::Area(
    double left_bnd_kev, double right_bnd_kev
) : data(SpectrumArea(left_bnd_kev, right_bnd_kev)) {}


std::array<double, 2> Area::to_bnds_kev(double eres) const {
    return std::visit(
        [this, eres](const auto& area) {
            return area.to_bnds_kev(eres);
        },
        data
    );
}


// class DBSpectrum

DBSpectrum::DBSpectrum(
    core::utils::Spectrum&& spectrum,
    core::utils::EnergyDetector detector
) : spectrum(spectrum), detector(detector) {
    counts = std::visit(
        [](const auto& spectrum) {
            return std::accumulate(
                spectrum.begin(),
                spectrum.end(),
                0,
                [](std::uint64_t acc, auto& x) { return acc + x; }
            );
        },
        this->spectrum.data
    );
}


Eigen::VectorXd DBSpectrum::get_peak_energies() {
    auto ecal = detector.corrected_ecal.value_or(detector.ecal);

    auto [left_peak_idx, right_peak_idx] = peak_bnd_idcs.value();

    size_t len = right_peak_idx - left_peak_idx + 1;

    auto linspace = Eigen::VectorXd::LinSpaced(len, left_peak_idx, right_peak_idx);

    return linspace.array() * ecal.scale + ecal.offset;
}


void DBSpectrum::fit_peak(PeakModel peak_model) {
    auto x = get_peak_energies();
    auto& y = peak.value();

    auto max = y.maxCoeff();

    if (peak_model == GAUSS) {
        peak_params = fitting::fit_gauss(x, y, {max, 511., 3.});
        return;
    }

    auto min = y.minCoeff();

    double first = *y.begin();
    double last = *std::prev(y.end());
    double erf_amp = 0.5 * (first - last);

    switch (peak_model) {
        case ERFLINEAR1GAUSS: {
            peak_params = fitting::fit_erf_linear_1_gauss(
                x, y, {max, 511., 1., erf_amp, 0., min}
            );
            break;
        }
        case ERFLINEAR2GAUSS: {
            peak_params = fitting::fit_erf_linear_2_gauss(
                x, y, {max / 2, 511., 1., max / 2, 511., 1.5, erf_amp, 0., min}
            );
            break;
        }
        case ERFLINEAR3GAUSS: {
            peak_params = fitting::fit_erf_linear_3_gauss(
                x, y, {
                    max / 3, 511., 1., max / 3, 511., 1.5, max / 3,
                    511., 2., erf_amp, 0., min
                }
            );
            break;
        }
        default: {
            throw std::runtime_error("Unreachable code reached");
        }
    }
}


void DBSpectrum::subtract_bg(PeakModel peak_model) {
    if (peak_model == PeakModel::GAUSS) return;

    fit_peak(peak_model);

    auto x = get_peak_energies();
    auto& y = peak.value();

    double amp = peak_params["erf_amp"].val;
    double x0 = peak_params["x0_1"].val;
    double sig = peak_params["sig_1"].val;
    double lin = peak_params["lin"].val;
    double off = peak_params["const"].val;

    auto corr = x.unaryExpr([amp, x0, sig, lin, off](double x_i) {
        return fitting::erf_linear_background(x_i, amp, x0, sig, lin, off);
    });

    for (auto&& [y_i, c_i] : std::views::zip(y, corr)) {
        y_i = std::max(y_i - c_i, 0.);
    }
}


void DBSpectrum::extract_peak(double peak_width, PeakModel peak_model)  {
    auto ecal = detector.corrected_ecal.value_or(detector.ecal);

    auto left_peak_idx = ecal.to_index_rounded(M_E_KEV - peak_width / 2);
    auto right_peak_idx = ecal.to_index_rounded(M_E_KEV + peak_width / 2);

    if (right_peak_idx >= spectrum.size()) {
        throw std::runtime_error(
            "Attempting to extract peak larger than spectrum"
        );
    }

    peak = std::visit(
        [left_peak_idx, right_peak_idx](const auto& spectrum) {
            auto peak = spectrum(Eigen::seq(left_peak_idx, right_peak_idx));
            return Eigen::VectorX<double>(peak.template cast<double>());
        },
        spectrum.data
    );

    peak_bnd_idcs = {left_peak_idx, right_peak_idx};

    subtract_bg(peak_model);
}


void DBSpectrum::correct_ecal(core::utils::EcalCorrectionOrder order) {
    if (order == core::utils::EcalCorrectionOrder::NONE) return;

    if (!peak.has_value()) {
        extract_peak(20., PeakModel::GAUSS);
    }

    fit_peak(PeakModel::GAUSS);

    double c = peak_params["x0_1"].val;

    auto ecal = detector.ecal;

    switch (order) {
        case core::utils::EcalCorrectionOrder::ZEROTH:
            ecal.offset += M_E_KEV - c;
            break;
        case core::utils::EcalCorrectionOrder::FIRST:
            ecal.scale *= (M_E_KEV - ecal.offset) / (c - ecal.offset);
            break;
        case core::utils::EcalCorrectionOrder::NONE:
            break;
    }

    detector.corrected_ecal = ecal;
}


double DBSpectrum::integrate(
    Area area,
    bool in_corrected_peak,
    BinIntegrationScheme bin_integration_scheme
) {
    auto [lower_bnd, upper_bnd] = area.to_bnds_kev(detector.eres);

    if (lower_bnd >= upper_bnd) {
        throw std::runtime_error(
            "Input parameter violates constaint lower bound < upper bound"
        );
    }

    auto ecal = detector.corrected_ecal.value_or(detector.ecal);

    double lower_ch_bnd = ecal.to_index_double(lower_bnd);
    double upper_ch_bnd = ecal.to_index_double(upper_bnd);

    if (in_corrected_peak) {
        if (!peak.has_value()) {
            throw std::runtime_error(
                "Attempting analysis on background subtracted peak, but no "
                "subtraction has been performed yet"
            );
        }

        auto [lower_peak_bnd, _] = peak_bnd_idcs.value();

        lower_ch_bnd -= lower_peak_bnd;
        upper_ch_bnd -= lower_peak_bnd;

        return integrate_generic(
            peak.value(), lower_ch_bnd, upper_ch_bnd, bin_integration_scheme
        );
    } else {
        return std::visit(
            [lower_ch_bnd, upper_ch_bnd, bin_integration_scheme](
                const auto& spectrum
            ){
                return integrate_generic(
                    spectrum, lower_ch_bnd, upper_ch_bnd, bin_integration_scheme
                );
            },
            spectrum.data
        );
    }
}


LineshapeParam& DBSpectrum::calc_lineshape_param(
    LineshapeParamDefinition definition,
    BinIntegrationScheme bin_integration_scheme
) {
    auto num_bnds = get_area_bnds(definition.num, detector.eres);
    auto denom_bnds = get_area_bnds(definition.denom, detector.eres);

    auto num_areas = definition.num
        | std::views::transform(
            [this, &definition, bin_integration_scheme](auto x) {
                return integrate(
                    x, definition.in_corrected_peak, bin_integration_scheme
                );
            }
        );

    double num = std::accumulate(num_areas.begin(), num_areas.end(), 0);

    auto denom_areas = definition.denom
        | std::views::transform(
            [this, &definition, bin_integration_scheme](auto x) {
                return integrate(
                    x, definition.in_corrected_peak, bin_integration_scheme
                );
            }
        );

    double denom = std::accumulate(denom_areas.begin(), denom_areas.end(), 0);

    double common = 0;

    for (const auto num_bnd : num_bnds) {
        for (const auto denom_bnd : denom_bnds) {
            double left_bnd_kev = std::max(num_bnd[0], denom_bnd[0]);
            double right_bnd_kev = std::min(num_bnd[1], denom_bnd[1]);
            if (left_bnd_kev < right_bnd_kev) {
                common += integrate(
                    Area(left_bnd_kev, right_bnd_kev),
                    definition.in_corrected_peak,
                    bin_integration_scheme
                );
            }
        }
    }

    double val = num / denom;
    double err = val * std::sqrt(
        (num - common) / std::pow(num, 2) +
        (denom - common) / std::pow(denom, 2) +
        common * std::pow((denom - num) / (denom * num), 2)
    );

    LineshapeParam ls_param {
        val, err, definition.num, definition.denom, definition.in_corrected_peak
    };

    lineshape_params[definition.name] = ls_param;

    return lineshape_params.at(definition.name);
}


void DBSpectrum::default_analyze() {
    correct_ecal(core::utils::EcalCorrectionOrder::FIRST);

    extract_peak(30, PeakModel::ERFLINEAR2GAUSS);

    for (const auto& param : STD_LINESHAPE_PARAMS) {
        calc_lineshape_param(param, BinIntegrationScheme::CONST);
    }

    peak = std::nullopt;
    peak_bnd_idcs = std::nullopt;
}



} // namespace compositron::dbs
