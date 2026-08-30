#include "compositron/cdbs/cdbspectrum.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <optional>
#include <ranges>
#include <stdexcept>
#include <fstream>
#include <string>

#include <ceres/types.h>

#include "compositron/cdbs/anti_aliasing.hpp"
#include "compositron/cdbs/fitting.hpp"
#include "compositron/constants.hpp"
#include "compositron/core/utils.hpp"
#include "compositron/dbs/dbspectrum.hpp"


namespace compositron::cdbs {


namespace {


inline double ij_to_xy(double i, core::utils::LinearCalibration ecal) {
    return ecal.from_index(i);
}


inline double xy_to_ij(double x, core::utils::LinearCalibration ecal) {
    return ecal.to_index_double(x);
}


inline std::array<double, 2> xy_to_uv(double x, double y) {
    return {0.5 * (x - y), 0.5 * (x + y) - M_E_KEV};
}


inline std::array<double, 2> uv_to_xy(double u, double v) {
    return {M_E_KEV + u + v, M_E_KEV - u + v};
}


inline std::array<double, 2> ij_to_uv(
    double i, double j, std::array<core::utils::LinearCalibration, 2> ecal
) {
    double x = ij_to_xy(j, ecal[1]);
    double y = ij_to_xy(i, ecal[0]);
    return xy_to_uv(x, y);
}


inline std::array<double, 2> uv_to_ij(
    double u, double v, std::array<core::utils::LinearCalibration, 2> ecal
) {
    auto [x, y] = uv_to_xy(u, v);
    return {xy_to_ij(y, ecal[0]), xy_to_ij(x, ecal[1])};
}


anti_aliasing::Rectangle convert_rectangle(
    std::array<double, 2> first_det_bnds,
    std::array<double, 2> second_det_bnds,
    std::array<core::utils::LinearCalibration, 2> ecal,
    std::array<std::array<size_t, 2>, 2> peak_bnds
) {
    return {
        ecal[0].to_index_double(first_det_bnds[1]) - peak_bnds[0][0],
        ecal[0].to_index_double(first_det_bnds[0]) - peak_bnds[0][0],
        ecal[1].to_index_double(second_det_bnds[0]) - peak_bnds[1][0],
        ecal[1].to_index_double(second_det_bnds[1]) - peak_bnds[1][0]
    };
}


anti_aliasing::Parallelogram convert_parallelogram(
    double width_cel,
    double width_cml,
    double offset_cel,
    double offset_cml,
    std::array<core::utils::LinearCalibration, 2> ecal,
    std::array<std::array<size_t, 2>, 2> peak_bnds
) {
    auto [i, j] = uv_to_ij(offset_cel, offset_cml, ecal);

    double w1 = std::sqrt(0.5 * (
        std::pow(width_cel / ecal[1].scale, 2)
        + std::pow(width_cel / ecal[0].scale, 2)
    ));
    double w2 = std::sqrt(0.5 * (
        std::pow(width_cml / ecal[1].scale, 2)
        + std::pow(-width_cml / ecal[0].scale, 2)
    ));

    return {
        i - peak_bnds[0][0],
        j - peak_bnds[1][0],
        -ecal[1].scale / ecal[0].scale,
        ecal[1].scale / ecal[0].scale,
        w1,
        w2
    };
}


anti_aliasing::Ellipse convert_ellipse(
    double radius_cel,
    double radius_cml,
    std::array<core::utils::LinearCalibration, 2> ecal,
    std::array<std::array<size_t, 2>, 2> peak_bnds
) {
    double w1 = std::sqrt(0.5 * (
        std::pow(radius_cel / ecal[1].scale, 2) + std::pow(radius_cel / ecal[0].scale, 2)
    ));
    double w2 = std::sqrt(0.5 * (
        std::pow(radius_cml / ecal[1].scale, 2) + std::pow(-radius_cml / ecal[0].scale, 2)
    ));

    return {
        ecal[0].to_index_double(M_E_KEV) - peak_bnds[0][0],
        ecal[1].to_index_double(M_E_KEV) - peak_bnds[1][0],
        w1,
        w2,
        -std::atan(ecal[1].scale / ecal[0].scale)
    };
}


} // anonymous namespace


// class Projection


Eigen::VectorXd Projection::get_energies() {
    if (ecal.has_value()) {
        size_t len = spectrum.size();
        auto e = ecal.value();
        return Eigen::VectorXd::NullaryExpr(len, [e](auto i) {
            return e.from_index((size_t)i);
        });
    } else if (bins.has_value()) {
        auto b = bins.value();
        size_t len = b.size();
        return Eigen::VectorXd::NullaryExpr(len, [b](auto i) {
            return 0.5 * (b[i] + b[i + 1]);
        });
    }

    throw std::runtime_error(
        "Could not fold projection as bins are not symmetric around 0"
    );
}


FoldedProjection Projection::fold() {
    if (ecal.has_value()) {
        auto e = ecal.value();
        double center_idx = e.to_index_double(0);

        if (std::abs(center_idx - (std::floor(center_idx) + 0.5)) > 1e-2) {
            throw std::runtime_error(
                "Could not fold projection as bins are not symmetric around 0"
            );
        }

        size_t last_left_idx = e.to_index_rounded(0);
        size_t first_right_idx = last_left_idx + 1;

        size_t len_left = last_left_idx + 1;
        size_t len_right = spectrum.size() - len_left;

        size_t len_new = std::min(len_left, len_right);

        auto folded_spectrum = Eigen::VectorXd::NullaryExpr(
            len_new, [this, last_left_idx, first_right_idx](auto i) {
                return spectrum[last_left_idx - i] + spectrum[first_right_idx + i];
            }
        );

        auto energies = Eigen::VectorXd::NullaryExpr(
            len_new, [e, first_right_idx](auto i) {
                return e.from_index(first_right_idx + i);
            }
        );

        return FoldedProjection { energies, folded_spectrum };
    } else if (bins.has_value()) {
        auto b = bins.value();
        auto center_idx = b.cwiseAbs().minCoeff();

        size_t len_left = center_idx;
        size_t len_right = b.size() - len_left - 1;

        size_t len_new = std::min(len_left, len_right);

        std::vector<double> folded_spectrum, energies;

        for (size_t i = 0; i < len_new; ++i) {
            auto left_bin_left_edge = b[center_idx - i - 1];
            auto left_bin_right_edge = b[center_idx - i];

            auto right_bin_left_edge = b[center_idx + i];
            auto right_bin_right_edge = b[center_idx + i + 1];

            if (
                left_bin_left_edge + right_bin_right_edge > 1e-2
                || left_bin_right_edge + right_bin_left_edge > 1e-2
            ) {
                throw std::runtime_error(
                    "Could not fold projection as bins are not symmetric around 0"
                );
            }

            folded_spectrum.push_back(
                spectrum[center_idx - i - 1] + spectrum[center_idx + i]
            );
            energies.push_back(
                0.5 * (right_bin_left_edge + right_bin_right_edge)
            );
        }

        auto e = energies;

        return FoldedProjection {
            Eigen::Map<Eigen::VectorXd>(e.data(), e.size()),
            Eigen::Map<Eigen::VectorXd>(folded_spectrum.data(), folded_spectrum.size())
        };
    }

    throw std::runtime_error(
        "Could not fold projection as bins are not symmetric around 0"
    );
}


dbs::DBSpectrum Projection::to_dbspectrum() {
    if (!ecal.has_value()) {
        throw std::runtime_error(
            "Could not convert projection to DBSpectrum due to missing "
            "energy calibration"
        );
    }

    auto s = spectrum.cast<uint64_t>();

    auto ss = std::vector<uint64_t>(s.begin(), s.end());

    return dbs::DBSpectrum(
        ss,
        core::utils::EnergyDetector {
            parent_detector_name,
            ecal.value(),
            std::nullopt,
            NAN
        }
    );
}


// class CDBSpectrum

void CDBSpectrum::default_analyze() {
    correct_ecal(core::utils::EcalCorrectionOrder::FIRST);

    extract_peak({10, 10}, BackgroundModel::NONE);

    for (const auto& param : STD_LINESHAPE_PARAMS) {
        calc_lineshape_param(param);
    }
}


Projection CDBSpectrum::project_axes(Axis axis) {
    switch (axis) {
        case FIRST_DET: {
            Eigen::VectorXd projected_spectrum = std::visit([](const auto& arr){
                return Eigen::VectorXd(arr.template cast<double>().rowwise().sum());
            }, spectrum.data);

            return Projection {
                projected_spectrum,
                detpair.first_det.name,
                detpair.first_det.corrected_ecal.value_or(detpair.first_det.ecal),
                std::nullopt
            };
        }
        case SECOND_DET: {
            auto projected_spectrum = std::visit([](const auto& arr){
                return Eigen::VectorXd(arr.template cast<double>().colwise().sum());
            }, spectrum.data);

            return Projection {
                projected_spectrum,
                detpair.second_det.name,
                detpair.second_det.corrected_ecal.value_or(detpair.second_det.ecal),
                std::nullopt
            };
        }
        default:
            throw std::runtime_error("Unreachable code reached");
    }
}


void CDBSpectrum::fit_2d_peak(BackgroundModel bg_model) {
    if (!peak.has_value() || !peak_bnds.has_value()) {
        throw std::runtime_error(
            "Attempting analysis on background subtracted peak, but no "
            "subtraction has been performed yet"
        );
    }

    auto p = peak.value();
    auto p_bnds = peak_bnds.value();

    auto ecal_1 = detpair.first_det.corrected_ecal.value_or(detpair.first_det.ecal);
    auto ecal_2 = detpair.second_det.corrected_ecal.value_or(detpair.second_det.ecal);

    size_t nrows = p_bnds[0][1] - p_bnds[0][0];
    size_t ncols = p_bnds[1][1] - p_bnds[1][0];

    auto x = Eigen::VectorXd::NullaryExpr(ncols, [ecal_2, p_bnds](auto i) {
        return ecal_2.from_index(p_bnds[1][0] + i);
    });
    auto y = Eigen::VectorXd::NullaryExpr(nrows, [ecal_1, p_bnds](auto i) {
        return ecal_1.from_index(p_bnds[0][0] + i);
    });

    switch (bg_model) {
        case NONE: {
            Eigen::Index x0, y0;
            double max = p.maxCoeff(&y0, &x0);

            size_t x0_idx = x0 + p_bnds[1][0];
            size_t y0_idx = y0 + p_bnds[0][0];

            peak_params = fitting::fit_gauss2d(
                x, y, p, {
                    max, ecal_2.from_index(x0_idx), ecal_1.from_index(y0_idx),
                    0.7, 1.8, -0.78
                }
            );

            break;
        }
        default:
            throw std::runtime_error("Unreachable code reached");
    }
}


void CDBSpectrum::subtract_bg(BackgroundModel bg_model) {
    if (bg_model == BackgroundModel::NONE) return;

    switch (bg_model) {
        default:
            throw std::runtime_error("Unreachable code reached");
    }
}


void CDBSpectrum::extract_peak(
    std::array<double, 2> peak_window_kev, BackgroundModel bg_model
) {
    auto ecal_1 = detpair.first_det.corrected_ecal.value_or(detpair.first_det.ecal);
    auto ecal_2 = detpair.second_det.corrected_ecal.value_or(detpair.second_det.ecal);

    double first_row_e = ecal_1.to_index_double(M_E_KEV - peak_window_kev[0] / 2);
    double last_row_e = ecal_1.to_index_double(M_E_KEV + peak_window_kev[0] / 2);

    double first_col_e = ecal_2.to_index_double(M_E_KEV - peak_window_kev[1] / 2);
    double last_col_e = ecal_2.to_index_double(M_E_KEV + peak_window_kev[1] / 2);

    auto [spectrum_nrows, spectrum_ncols] = std::visit([](const auto& spectrum) {
        return std::make_pair(spectrum.rows(), spectrum.cols());
    }, spectrum.data);

    size_t first_row = std::max(std::min(first_row_e, spectrum_nrows - 1.), 0.);
    size_t last_row = std::max(std::min(last_row_e + 1, spectrum_nrows - 1.), 0.);

    size_t first_col = std::max(std::min(first_col_e, spectrum_ncols - 1.), 0.);
    size_t last_col = std::max(std::min(last_col_e + 1, spectrum_ncols - 1.), 0.);

    size_t nrows = last_row - first_row;
    size_t ncols = last_col - first_col;

    peak = std::visit([first_row, first_col, nrows, ncols](const auto& spectrum) {
        return Eigen::MatrixXd(
            spectrum.block(first_row, first_col, nrows, ncols).template cast<double>()
        );
    }, spectrum.data);

    peak_bnds = {{{first_row, last_row}, {first_col, last_col}}};

    subtract_bg(bg_model);

    peak_counts = peak.value().sum();
}


void CDBSpectrum::correct_ecal(core::utils::EcalCorrectionOrder order) {
    if (order == core::utils::EcalCorrectionOrder::NONE) return;

    if (!peak.has_value()) {
        extract_peak({4, 4}, BackgroundModel::NONE);
    }

    if (!peak_params.contains("x0") || !peak_params.contains("y0")) {
        fit_2d_peak(BackgroundModel::NONE);
    }

    double x0 = peak_params.at("x0").val;
    double y0 = peak_params.at("y0").val;

    auto ecal_1 = detpair.first_det.ecal;
    auto ecal_2 = detpair.second_det.ecal;

    switch (order) {
        case core::utils::ZEROTH: {
            ecal_1.offset += M_E_KEV - y0;
            ecal_2.offset += M_E_KEV - x0;
            break;
        }
        case core::utils::FIRST: {
            ecal_1.scale *= (M_E_KEV - ecal_1.offset) / (y0 - ecal_1.offset);
            ecal_2.scale *= (M_E_KEV - ecal_2.offset) / (x0 - ecal_2.offset);
            break;
        }
        default:
            throw std::runtime_error("Unreachable code reached");
    }

    detpair.first_det.corrected_ecal = ecal_1;
    detpair.second_det.corrected_ecal = ecal_2;
}


double CDBSpectrum::get_max_projection_length(
    double v,
    std::array<core::utils::LinearCalibration, 2> ecal,
    std::array<std::array<size_t, 2>, 2> peak_bnds
) {
    auto [i1, j1] = uv_to_ij(0, v, ecal);
    auto [i2, j2] = uv_to_ij(0, -v, ecal);
    double m = -ecal[1].scale / ecal[0].scale;
    double t1 = i1 - m * j1;
    double t2 = i2 - m * j2;

    auto [peak_row_bnds, peak_col_bnds] = peak_bnds;
    auto [i_min_idx, i_max_idx] = peak_row_bnds;
    auto [j_min_idx, j_max_idx] = peak_col_bnds;

    double i_min = i_min_idx - 0.49;
    double i_max = i_max_idx + 0.49;
    double j_min = j_min_idx - 0.49;
    double j_max = j_max_idx + 0.49;

    double x1 = std::abs(ij_to_uv(i_min, (i_min - t1) / m, ecal)[0]);
    double x2 = std::abs(ij_to_uv(i_min, (i_min - t2) / m, ecal)[0]);
    double x3 = std::abs(ij_to_uv(m * j_min + t1, j_min, ecal)[0]);
    double x4 = std::abs(ij_to_uv(m * j_min + t2, j_min, ecal)[0]);
    double x5 = std::abs(ij_to_uv(i_max, (i_max - t1) / m, ecal)[0]);
    double x6 = std::abs(ij_to_uv(i_max, (i_max - t2) / m, ecal)[0]);
    double x7 = std::abs(ij_to_uv(m * j_max + t1, j_max, ecal)[0]);
    double x8 = std::abs(ij_to_uv(m * j_max + t2, j_max, ecal)[0]);

    return std::min({x1, x2, x3, x4, x5, x6, x7, x8});
}


anti_aliasing::Polygon CDBSpectrum::calculate_polygon_boundaries(
    const DiagonalArea& area,
    std::array<core::utils::LinearCalibration, 2> ecal,
    std::array<std::array<size_t, 2>, 2> peak_bnds
) {
    double width_cel = area.width_cel.to_kev(detpair.eres);
    double width_cml = area.width_cml.to_kev(detpair.eres);

    if (std::isinf(width_cel)) {
        width_cel = get_max_projection_length(0.5 * width_cml, ecal, peak_bnds);
    }

    if (!std::isfinite(width_cel) || !std::isfinite(width_cml)) {
        throw std::runtime_error("Invalid boundaries for area");
    }

    return anti_aliasing::Polygon { convert_parallelogram(
        width_cel, width_cml, 0, 0, ecal, peak_bnds
    ).to_vertices() };
}


anti_aliasing::Polygon CDBSpectrum::calculate_polygon_boundaries(
    const DiagonalAreaWithOffset& area,
    std::array<core::utils::LinearCalibration, 2> ecal,
    std::array<std::array<size_t, 2>, 2> peak_bnds
) {
    double width_cel = area.width_cel.to_kev(detpair.eres);
    double width_cml = area.width_cml.to_kev(detpair.eres);
    double offset_cel = area.offset_cel.to_kev(detpair.eres);
    double offset_cml = area.offset_cml.to_kev(detpair.eres);

    if (std::isinf(width_cel)) {
        width_cel = get_max_projection_length(0.5 * width_cml, ecal, peak_bnds);
    }

    if (
        !std::isfinite(width_cel) || !std::isfinite(width_cml)
        || !std::isfinite(offset_cel) || !std::isfinite(offset_cml)
    ) {
        throw std::runtime_error("Invalid boundaries for area");
    }

    return anti_aliasing::Polygon { convert_parallelogram(
        width_cel, width_cml, offset_cel, offset_cml, ecal, peak_bnds
    ).to_vertices() };
}


anti_aliasing::Polygon CDBSpectrum::calculate_polygon_boundaries(
    const AxisAlignedArea& area,
    std::array<core::utils::LinearCalibration, 2> ecal,
    std::array<std::array<size_t, 2>, 2> peak_bnds
) {
    auto [y_min, y_max] = area.first_det_bnds;
    auto [x_min, x_max] = area.second_det_bnds;

    if (std::isinf(y_min)) {
        y_min = ecal[0].from_index(peak_bnds[0][0] - 0.49);
    }
    if (std::isinf(y_max)) {
        y_max = ecal[0].from_index(peak_bnds[0][1] + 0.49);
    }
    if (std::isinf(x_min)) {
        x_min = ecal[1].from_index(peak_bnds[1][0] - 0.49);
    }
    if (std::isinf(x_max)) {
        x_max = ecal[1].from_index(peak_bnds[1][1] + 0.49);
    }

    if (
        !std::isfinite(y_min) || !std::isfinite(y_max)
        || ! std::isfinite(x_min) || !std::isfinite(x_max)
    ) {
        throw std::runtime_error("Invalid boundaries for area");
    }

    return anti_aliasing::Polygon { convert_rectangle(
        {y_min, y_max}, {x_min, x_max}, ecal, peak_bnds
    ).to_vertices() };
}


anti_aliasing::Polygon CDBSpectrum::calculate_polygon_boundaries(
    const EllipseArea& area,
    std::array<core::utils::LinearCalibration, 2> ecal,
    std::array<std::array<size_t, 2>, 2> peak_bnds
) {
    double radius_cel = area.radius_cel.to_kev(detpair.eres);
    double radius_cml = area.radius_cml.to_kev(detpair.eres);

    if (!std::isfinite(radius_cel) || !std::isfinite(radius_cml)) {
        throw std::runtime_error("Invalid boundaries for area");
    }

    return anti_aliasing::Polygon { convert_ellipse(
        radius_cel, radius_cml, ecal, peak_bnds
    ).to_vertices() };
}


double CDBSpectrum::blend_and_sum(
    const Eigen::MatrixX<uint8_t>& weights,
    size_t upper_row,
    size_t left_col,
    size_t nrows,
    size_t ncols,
    bool in_corrected_peak
) {
    if (in_corrected_peak) {
        auto p = peak.value();

        auto view = p.block(upper_row, left_col, nrows, ncols);

        return view.cwiseProduct(weights.cast<double>()).sum() / 255.;
    } else {
        return std::visit(
            [upper_row, left_col, nrows, ncols, weights](const auto& spectrum) {
                auto view = spectrum.block(upper_row, left_col, nrows, ncols);

                return view.template cast<uint64_t>()
                    .cwiseProduct(weights.cast<uint64_t>())
                    .sum() / 255.;
            }, spectrum.data
        );
    }
}


double CDBSpectrum::integrate(const DiagonalArea& area, bool in_corrected_peak) {
    return CDBSpectrum::integrate_generic(area, in_corrected_peak);
}


double CDBSpectrum::integrate(const DiagonalAreaWithOffset& area, bool in_corrected_peak) {
    return CDBSpectrum::integrate_generic(area, in_corrected_peak);
}


double CDBSpectrum::integrate(const AxisAlignedArea& area, bool in_corrected_peak) {
    return CDBSpectrum::integrate_generic(area, in_corrected_peak);
}


double CDBSpectrum::integrate(const EllipseArea& area, bool in_corrected_peak) {
    return CDBSpectrum::integrate_generic(area, in_corrected_peak);
}


LineshapeParam& CDBSpectrum::calc_lineshape_param(
    const LineshapeParamDefinition& definition
) {
    auto ecal_1 = detpair.first_det.corrected_ecal.value_or(detpair.first_det.ecal);
    auto ecal_2 = detpair.second_det.corrected_ecal.value_or(detpair.second_det.ecal);
    std::array<core::utils::LinearCalibration, 2> ecal = {ecal_1, ecal_2};

    auto [spectrum_nrows, spectrum_ncols] = std::visit([](const auto& spectrum) {
        return std::make_pair<size_t, size_t>(spectrum.rows(), spectrum.cols());
    }, spectrum.data);

    std::array<std::array<size_t, 2>, 2> integral_bnds;

    if (definition.in_corrected_peak) {
        if (peak_bnds.has_value()) {
            integral_bnds = peak_bnds.value();
        } else {
            throw std::runtime_error(
                "Attempting analysis on background subtracted peak, but no "
                "subtraction has been performed yet"
            );
        }
    } else {
        integral_bnds = {{{0, spectrum_nrows - 1}, {0, spectrum_ncols - 1}}};
    }

    size_t nrows = integral_bnds[0][1] - integral_bnds[0][0];
    size_t ncols = integral_bnds[1][1] - integral_bnds[1][0];

    auto num_polygons = definition.num
        | std::views::transform([this, ecal, integral_bnds](const Area& area) {
            return std::visit([this, ecal, integral_bnds](const auto& area) {
                return calculate_polygon_boundaries(area, ecal, integral_bnds);
            }, area);
        })
        | std::ranges::to<std::vector<anti_aliasing::Polygon>>();

    auto denom_polygons = definition.denom
        | std::views::transform([this, ecal, integral_bnds](const Area& area) {
            return std::visit([this, ecal, integral_bnds](const auto& area) {
                return calculate_polygon_boundaries(area, ecal, integral_bnds);
            }, area);
        })
        | std::ranges::to<std::vector<anti_aliasing::Polygon>>();

    for (const auto& field : {num_polygons, denom_polygons}) {
        int len = field.size();
        for (auto i : std::views::iota(0, len - 1)) {
            for (auto j : std::views::iota(i + 1, len)) {
                auto intersection = anti_aliasing::intersect_convex_polygons(
                    field[i], field[j]
                );

                if (intersection.vertices.size() > 0) {
                    throw std::runtime_error(
                        "Lineshape numerator or denominator areas overlap"
                    );
                }
            }
        }
    }

    auto num_areas = definition.num
        | std::views::transform([this, definition](const Area& area) {
            return std::visit([this, definition](const auto& area) {
                return integrate(area, definition.in_corrected_peak);
            }, area);
        })
        | std::ranges::to<std::vector<double>>();

    auto denom_areas = definition.denom
        | std::views::transform([this, definition](const Area& area) {
            return std::visit([this, definition](const auto& area) {
                return integrate(area, definition.in_corrected_peak);
            }, area);
        })
        | std::ranges::to<std::vector<double>>();

    std::vector<double> common_areas;

    for (const auto& num_area : num_polygons) {
        for (const auto& denom_area : denom_polygons) {
            auto intersection = anti_aliasing::intersect_convex_polygons(
                num_area, denom_area
            );

            if (intersection.vertices.size() == 0) continue;

            auto vertices = intersection.to_vertices();

            auto [left_col, right_col, lower_row, upper_row] =
                anti_aliasing::get_bounding_box(vertices);

            if (right_col > spectrum_ncols - 1 || lower_row > nrows - 1) {
                throw std::runtime_error("Polygon vertices are out of matrix bounds");
            }

            std::transform(
                vertices.begin(), vertices.end(), intersection.vertices.begin(),
                [left_col, upper_row](auto& v) {
                    return anti_aliasing::Vertex { v[0] - left_col, v[1] - upper_row };
                }
            );

            size_t nrows_view = lower_row - upper_row + 1;
            size_t ncols_view = right_col - left_col + 1;

            auto weights = anti_aliasing::agg_aa(nrows_view, ncols_view, intersection);

            common_areas.push_back(blend_and_sum(
                weights, upper_row, left_col,
                nrows_view, ncols_view, definition.in_corrected_peak
            ));
        }
    }

    double num = std::ranges::fold_left(num_areas, 0., std::plus{});
    double denom = std::ranges::fold_left(denom_areas, 0., std::plus{});
    double common = std::ranges::fold_left(common_areas, 0., std::plus{});

    double val = num / denom;
    double err = val * std::sqrt(
        (num - common) / std::pow(num, 2) +
        (denom - common) / std::pow(denom, 2) +
        common * std::pow((denom - num) / (denom * num), 2)
    );

    LineshapeParam param {
        val, err, definition.num, definition.denom, definition.in_corrected_peak
    };

    lineshape_params[definition.name] = param;

    return lineshape_params.at(definition.name);
}


Projection CDBSpectrum::project_diagonal(
    core::utils::Unit bin_width,
    core::utils::Unit width,
    core::utils::Unit max_length,
    bool in_corrected_peak
) {
    auto ecal_1 = detpair.first_det.corrected_ecal.value_or(detpair.first_det.ecal);
    auto ecal_2 = detpair.second_det.corrected_ecal.value_or(detpair.second_det.ecal);
    std::array<core::utils::LinearCalibration, 2> ecal = {ecal_1, ecal_2};

    auto [spectrum_nrows, spectrum_ncols] = std::visit([](const auto& spectrum) {
        return std::make_pair<size_t, size_t>(spectrum.rows(), spectrum.cols());
    }, spectrum.data);

    std::array<std::array<size_t, 2>, 2> integral_bnds;

    if (in_corrected_peak) {
        if (peak_bnds.has_value()) {
            integral_bnds = peak_bnds.value();
        } else {
            throw std::runtime_error(
                "Attempting analysis on background subtracted peak, but no "
                "subtraction has been performed yet"
            );
        }
    } else {
        integral_bnds = {{{0, spectrum_nrows - 1}, {0, spectrum_ncols - 1}}};
    }

    double width_kev = width.to_kev(detpair.eres);
    double max_length_kev = max_length.to_kev(detpair.eres);

    double u_bin = bin_width.to_kev(detpair.eres);

    double u_max = get_max_projection_length(width_kev / 2, ecal, integral_bnds);
    u_max = std::min(u_max, max_length_kev);

    size_t num_bins = std::max(2 * (u_max / u_bin), 0.);

    if (num_bins < 2) {
        throw std::runtime_error("Invalid bins for projection");
    }

    double max_offset = (num_bins / 2 - 0.5) * u_bin;

    core::utils::LinearCalibration proj_ecal{-max_offset, u_bin};

    auto projection = Eigen::VectorXd::NullaryExpr(
        num_bins, [this, u_bin, width, proj_ecal, in_corrected_peak](auto i) {
            return integrate(
                DiagonalAreaWithOffset(
                    core::utils::Unit(u_bin, core::utils::KEV),
                    width,
                    core::utils::Unit(proj_ecal.from_index((size_t)i), core::utils::KEV),
                    core::utils::Unit(u_bin, core::utils::KEV)
                ),
                in_corrected_peak
            );
        }
    );

    return Projection {projection, detpair.name, proj_ecal, std::nullopt};
}


Projection CDBSpectrum::project_diagonal(
    std::vector<core::utils::Unit> bins,
    core::utils::Unit width,
    core::utils::Unit max_length,
    bool in_corrected_peak
) {
    auto ecal_1 = detpair.first_det.corrected_ecal.value_or(detpair.first_det.ecal);
    auto ecal_2 = detpair.second_det.corrected_ecal.value_or(detpair.second_det.ecal);
    std::array<core::utils::LinearCalibration, 2> ecal = {ecal_1, ecal_2};

    auto [spectrum_nrows, spectrum_ncols] = std::visit([](const auto& spectrum) {
        return std::make_pair<size_t, size_t>(spectrum.rows(), spectrum.cols());
    }, spectrum.data);

    std::array<std::array<size_t, 2>, 2> integral_bnds;

    if (in_corrected_peak) {
        if (peak_bnds.has_value()) {
            integral_bnds = peak_bnds.value();
        } else {
            throw std::runtime_error(
                "Attempting analysis on background subtracted peak, but no "
                "subtraction has been performed yet"
            );
        }
    } else {
        integral_bnds = {{{0, spectrum_nrows - 1}, {0, spectrum_ncols - 1}}};
    }

    double width_kev = width.to_kev(detpair.eres);
    double max_length_kev = max_length.to_kev(detpair.eres);

    if (bins.size() < 2) {
        throw std::runtime_error("Invalid bins for projection");
    }

    auto bins_kev = bins
        | std::views::transform([this](const auto& bin) {
            return bin.to_kev(detpair.eres);
        })
        | std::ranges::to<std::vector<double>>();

    auto projection = Eigen::VectorXd::NullaryExpr(
        bins.size(), [this, bins_kev, in_corrected_peak](auto i) {
            double left_edge = bins_kev[i];
            double right_edge = bins_kev[i + 1];

            double width_cel = right_edge - left_edge;
            double offset = 0.5 * (left_edge + right_edge);

            return integrate(
                DiagonalAreaWithOffset(
                    core::utils::Unit(width_cel, core::utils::KEV),
                    core::utils::Unit(0, core::utils::KEV),
                    core::utils::Unit(offset, core::utils::KEV),
                    core::utils::Unit(0, core::utils::KEV)
                ),
                in_corrected_peak
            );
        }
    );

    return Projection {
        projection,
        detpair.name,
        std::nullopt,
        Eigen::Map<Eigen::VectorXd>(bins_kev.data(), bins_kev.size())
    };
}



} // namespace compositron::cdbs
