#pragma once

#include <ceres/jet.h>

#include "compositron/constants.hpp"
#include "compositron/pals/model.hpp"


namespace compositron::pals::fitting {


struct FitResult {
    model::LifetimeModel model;
    ceres::Solver::Summary summary;
    ceres::TerminationType fit_status;
    int n_eval;
    std::optional<double> red_chi_2;
    std::optional<Eigen::MatrixXd> cov;
};


class LifetimeResidual {
public:
    double x;
    double y;
    double w;
    size_t n_l;
    size_t n_r;
    int last_l_vary_idx;
    int last_r_vary_idx;
    const std::vector<uint8_t> l_vary;
    const std::vector<uint8_t> r_vary;


    // TODO: outsource into a callback to only compute once per evaluation cycle
    // (requires implementing an analytic jacobian)
    template <typename T>
    std::array<std::vector<T>, 2> transform_intensities(
        T const* const* params
    ) const {
        T available_l_intensity(1);
        T available_r_intensity(1);

        std::vector<T> abs_l_intensities;
        std::vector<T> abs_r_intensities;

        for (int i = 0; i < n_l; ++i) {
            T l_int;
            if (l_vary[i]) {
                l_int = T(available_l_intensity * params[3 + 2 * i + 1][0]);
            } else if (i == last_l_vary_idx) {
                l_int = T(available_l_intensity);
            } else {
                l_int = T(params[3 + 2 * i + 1][0]);
            }
            abs_l_intensities.push_back(l_int);
            available_l_intensity -= l_int;
        }

        for (int i = 0; i < n_r; ++i) {
            T r_int;
            if (r_vary[i]) {
                r_int = T(available_r_intensity * params[3 + 2 * n_l + 3 * i + 1][0]);
            } else if (i == last_r_vary_idx) {
                r_int = T(available_r_intensity);
            } else {
                r_int = T(params[3 + 2 * n_l + 3 * i + 1][0]);
            }
            abs_r_intensities.push_back(r_int);
            available_r_intensity -= r_int;
        }

        return {abs_l_intensities, abs_r_intensities};
    }


    template<typename T>
    bool operator()(T const* const* params, T* residuals) const {
        T f(0);

        T n = params[0][0];
        T bg = params[1][0];
        T t0 = params[2][0];

        auto [abs_l_intensities, abs_r_intensities] = transform_intensities(params);

        for (size_t i = 0; i < n_l; ++i) {
            T tau = params[3 + 2 * i][0];
            T l_int = abs_l_intensities[i];

            for (size_t j = 0; j < n_r; ++j) {
                T sig = params[3 + 2 * n_l + 3 * j][0];
                T r_int = abs_r_intensities[j];
                T r_t0 = params[3 + 2 * n_l + 3 * j + 2][0];

                T u = sig / (tau * SQRT_2);

                T x_i = x - t0 - r_t0;

                T e = ceres::exp(-x_i / tau + u * u);
                T c = ceres::erfc(u - x_i / (sig * SQRT_2));

                f += 0.5 * l_int * r_int / tau * e * c;
            }
        }

        if (!ceres::isfinite(f)) {
            return false;
        }

        residuals[0] = w * (y - (n * f + bg));

        return true;
    }
};


FitResult fit_lifetime_model(
    const Eigen::VectorXd& x,
    const Eigen::VectorXd& y,
    const model::LifetimeModel& model,
    uint64_t counts,
    double peak_center
);


} // namespace compositron::pals::fitting
