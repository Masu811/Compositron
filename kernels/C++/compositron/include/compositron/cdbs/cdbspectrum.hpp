#pragma once

#include <cmath>
#include <cstdint>
#include <optional>
#include <variant>
#include <vector>
#include <string>

#include "compositron/core/utils.hpp"
#include "compositron/core/fitting.hpp"
#include "compositron/cdbs/anti_aliasing.hpp"
#include "compositron/dbs/dbspectrum.hpp"


namespace compositron::cdbs {


enum BackgroundModel {
    NONE,
};


enum Axis {
    FIRST_DET,
    SECOND_DET
};


class DiagonalArea {
public:
    core::utils::Unit width_cel;
    core::utils::Unit width_cml;

    DiagonalArea(
        core::utils::Unit width_cel,
        core::utils::Unit width_cml
    ) : width_cel(width_cel),
        width_cml(width_cml) {}
};


class DiagonalAreaWithOffset {
public:
    core::utils::Unit width_cel;
    core::utils::Unit width_cml;
    core::utils::Unit offset_cel;
    core::utils::Unit offset_cml;

    DiagonalAreaWithOffset(
        core::utils::Unit width_cel,
        core::utils::Unit width_cml,
        core::utils::Unit offset_cel,
        core::utils::Unit offset_cml
    ) : width_cel(width_cel),
        width_cml(width_cml),
        offset_cel(offset_cel),
        offset_cml(offset_cml) {}
};

class AxisAlignedArea {
public:
    std::array<double, 2> first_det_bnds;
    std::array<double, 2> second_det_bnds;

    AxisAlignedArea(
        std::array<double, 2> first_det_bnds,
        std::array<double, 2> second_det_bnds
    ) : first_det_bnds(first_det_bnds), second_det_bnds(second_det_bnds) {}
};


class EllipseArea {
public:
    core::utils::Unit radius_cel;
    core::utils::Unit radius_cml;
};


using Area = std::variant<
    DiagonalArea,
    DiagonalAreaWithOffset,
    AxisAlignedArea,
    EllipseArea
>;


class LineshapeParam {
public:
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


const std::array<LineshapeParamDefinition, 3> STD_LINESHAPE_PARAMS{{
    {
        "S",
        {
            DiagonalArea(
                core::utils::Unit(2, core::utils::KEV),
                core::utils::Unit(1, core::utils::ERES)
            )
        },
        {
            DiagonalArea(
                core::utils::Unit(INFINITY, core::utils::KEV),
                core::utils::Unit(1, core::utils::ERES)
            )
        },
        true
    },
    {
        "W",
        {
            DiagonalAreaWithOffset(
                core::utils::Unit(1, core::utils::KEV),
                core::utils::Unit(1, core::utils::ERES),
                core::utils::Unit(3, core::utils::KEV),
                core::utils::Unit(0, core::utils::KEV)
            ),
            DiagonalAreaWithOffset(
                core::utils::Unit(1, core::utils::KEV),
                core::utils::Unit(1, core::utils::ERES),
                core::utils::Unit(-3, core::utils::KEV),
                core::utils::Unit(0, core::utils::KEV)
            )
        },
        {
            DiagonalArea(
                core::utils::Unit(INFINITY, core::utils::KEV),
                core::utils::Unit(1, core::utils::ERES)
            )
        },
        true
    },
    {
        "P/T",
        {
            DiagonalArea(
                core::utils::Unit(INFINITY, core::utils::KEV),
                core::utils::Unit(1, core::utils::ERES)
            )
        },
        {
            AxisAlignedArea({-INFINITY, INFINITY}, {-INFINITY, INFINITY})
        },
        false
    },
}};


struct FoldedProjection {
    Eigen::VectorXd energies;
    Eigen::VectorXd spectrum;
};


class Projection {
public:
    Eigen::VectorXd spectrum;
    std::string parent_detector_name;
    std::optional<core::utils::LinearCalibration> ecal;
    std::optional<Eigen::VectorXd> bins;
    double counts;

    Projection(
        Eigen::VectorXd spectrum,
        std::string parent_detector_name,
        std::optional<core::utils::LinearCalibration> ecal,
        std::optional<Eigen::VectorXd> bins
    ) : spectrum(spectrum),
        parent_detector_name(parent_detector_name),
        ecal(ecal),
        bins(bins) {}

    Eigen::VectorXd get_energies();

    FoldedProjection fold();

    dbs::DBSpectrum to_dbspectrum();
};


class CDBSpectrum {
public:
    compositron::core::utils::Spectrum2D spectrum;
    compositron::core::utils::EnergyDetectorPair detpair;
    uint64_t counts = 0;
    std::optional<Eigen::MatrixXd> peak;
    std::optional<std::array<std::array<size_t, 2>, 2>> peak_bnds;
    std::optional<double> peak_counts;
    std::unordered_map<std::string, core::fitting::SimpleFitParam> peak_params;
    std::unordered_map<std::string, LineshapeParam> lineshape_params;

    CDBSpectrum() = delete;

    CDBSpectrum(CDBSpectrum& other) = delete;

    CDBSpectrum(CDBSpectrum&& other) = default;

    CDBSpectrum(
        core::utils::Spectrum2D&& spectrum,
        core::utils::EnergyDetectorPair detpair
    ) : spectrum(spectrum), detpair(detpair) {}

    template<typename T>
    CDBSpectrum(
        const std::vector<T>& spectrum,
        size_t nrows,
        size_t ncols,
        core::utils::StorageOrder order,
        core::utils::EnergyDetectorPair detpair
    ) : CDBSpectrum(
            core::utils::Spectrum2D(spectrum, nrows, ncols, order),
            detpair
    ) {}

    Projection project_axes(Axis axis);

private:
    void fit_2d_peak(BackgroundModel bg_model);

    void subtract_bg(BackgroundModel bg_model);

public:
    void extract_peak(
        std::array<double, 2> peak_window_kev, BackgroundModel bg_model
    );

    void correct_ecal(core::utils::EcalCorrectionOrder order);

private:
    double get_max_projection_length(
        double v,
        std::array<core::utils::LinearCalibration, 2> ecal,
        std::array<std::array<size_t, 2>, 2> peak_bnds
    );

    anti_aliasing::Polygon calculate_polygon_boundaries(
        const DiagonalArea& area,
        std::array<core::utils::LinearCalibration, 2> ecal,
        std::array<std::array<size_t, 2>, 2> peak_bnds
    );

    anti_aliasing::Polygon calculate_polygon_boundaries(
        const DiagonalAreaWithOffset& area,
        std::array<core::utils::LinearCalibration, 2> ecal,
        std::array<std::array<size_t, 2>, 2> peak_bnds
    );

    anti_aliasing::Polygon calculate_polygon_boundaries(
        const AxisAlignedArea& area,
        std::array<core::utils::LinearCalibration, 2> ecal,
        std::array<std::array<size_t, 2>, 2> peak_bnds
    );

    anti_aliasing::Polygon calculate_polygon_boundaries(
        const EllipseArea& area,
        std::array<core::utils::LinearCalibration, 2> ecal,
        std::array<std::array<size_t, 2>, 2> peak_bnds
    );

    double blend_and_sum(
        const Eigen::MatrixX<uint8_t>& weights,
        size_t upper_row,
        size_t left_col,
        size_t nrows,
        size_t ncols,
        bool in_corrected_peak
    );

    template <typename T>
    double integrate_generic(const T& area, bool in_corrected_peak) {
        auto ecal_1 = detpair.first_det.corrected_ecal.value_or(detpair.first_det.ecal);
        auto ecal_2 = detpair.second_det.corrected_ecal.value_or(detpair.second_det.ecal);
        std::array<core::utils::LinearCalibration, 2> ecal = {ecal_1, ecal_2};

        size_t spectrum_nrows = spectrum.rows();
        size_t spectrum_ncols = spectrum.cols();

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

        size_t nrows = integral_bnds[0][1] - integral_bnds[0][0] + 1;
        size_t ncols = integral_bnds[1][1] - integral_bnds[1][0] + 1;

        auto polygon = calculate_polygon_boundaries(area, ecal, integral_bnds);

        auto vertices = polygon.to_vertices();

        auto [left_col, right_col, lower_row, upper_row] =
            anti_aliasing::get_bounding_box(vertices);

        if (right_col > spectrum_ncols - 1 || lower_row > spectrum_nrows - 1) {
            throw std::runtime_error("Polygon vertices are out of matrix bounds");
        }

        std::transform(
            polygon.vertices.begin(), polygon.vertices.end(), polygon.vertices.begin(),
            [left_col, upper_row](auto& v) {
                return anti_aliasing::Vertex{v[0] - left_col, v[1] - upper_row};
            }
        );

        size_t nrows_view = lower_row - upper_row + 1;
        size_t ncols_view = right_col - left_col + 1;

        auto weights = anti_aliasing::agg_aa(nrows_view, ncols_view, polygon);

        return blend_and_sum(
            weights, upper_row, left_col, nrows_view, ncols_view, in_corrected_peak
        );
    }

public:
    double integrate(const DiagonalArea& area, bool in_corrected_peak);

    double integrate(const DiagonalAreaWithOffset& area, bool in_corrected_peak);

    double integrate(const AxisAlignedArea& area, bool in_corrected_peak);

    double integrate(const EllipseArea& area, bool in_corrected_peak);

    LineshapeParam& calc_lineshape_param(
        const LineshapeParamDefinition& definition
    );

    void default_analyze();

    Projection project_diagonal(
        core::utils::Unit bin_width,
        core::utils::Unit width,
        core::utils::Unit max_length,
        bool in_corrected_peak
    );

    Projection project_diagonal(
        std::vector<core::utils::Unit> bins,
        core::utils::Unit width,
        core::utils::Unit max_length,
        bool in_corrected_peak
    );
};


} // namespace compositron::cdbs
