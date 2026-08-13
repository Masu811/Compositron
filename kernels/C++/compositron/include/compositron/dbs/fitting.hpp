#pragma once

#include <ceres/autodiff_cost_function.h>
#include <ceres/cost_function.h>
#include <ceres/problem.h>
#include <ceres/solver.h>

#include "compositron/core/fitting.hpp"


namespace compositron::dbs::fitting {


double gauss(double x, double amp, double x0, double sig);


struct GaussResidual {
    const double x;
    const double y;
    const double w;

    GaussResidual(double x, double y, double w) : x(x), y(y), w(w) {}

    template <typename T>
    bool operator()(
        const T* const amp,
        const T* const x0,
        const T* const sig,
        T* residual
    ) const {
        if (sig[0] == 0) {
            return false;
        }

        T u = (x - x0[0]) / sig[0];

        residual[0] = w * (y - amp[0] * ceres::exp(-0.5 * u * u));

        return true;
    }
};

core::fitting::FitResult fit_gaussian(
    const Eigen::VectorXd& x,
    const Eigen::VectorXd& y,
    const std::array<double, 3>& init
);


double erf_linear_background(
    double x, double amp, double x0, double sig, double lin, double off
);


struct ErfLinear1GaussResidual {
    const double x;
    const double y;
    const double w;

    ErfLinear1GaussResidual(double x, double y, double w) : x(x), y(y), w(w) {}

    template <typename T>
    bool operator()(
        const T* const amp_gauss_1,
        const T* const x0_1,
        const T* const sig_1,
        const T* const amp_erf,
        const T* const lin,
        const T* const off,
        T* residual
    ) const {
        T u_1 = (x - x0_1[0]) / sig_1[0];

        T gauss_1 = amp_gauss_1[0] * ceres::exp(-0.5 * u_1 * u_1);

        T erf = amp_erf[0] * ceres::erf(u_1);

        residual[0] = w * (y - (gauss_1 + erf + lin[0] * x + off[0]));

        return true;
    }
};


core::fitting::FitResult fit_erf_linear_1_gauss(
    const Eigen::VectorXd& x,
    const Eigen::VectorXd& y,
    const std::array<double, 6>& init
);


struct ErfLinear2GaussResidual {
    const double x;
    const double y;
    const double w;

    ErfLinear2GaussResidual(double x, double y, double w) : x(x), y(y), w(w) {}

    template <typename T>
    bool operator()(
        const T* const amp_gauss_1,
        const T* const x0_1,
        const T* const sig_1,
        const T* const amp_gauss_2,
        const T* const x0_2,
        const T* const sig_2,
        const T* const amp_erf,
        const T* const lin,
        const T* const off,
        T* residual
    ) const {
        T u_1 = (x - x0_1[0]) / sig_1[0];
        T u_2 = (x - x0_2[0]) / sig_2[0];

        T gauss_1 = amp_gauss_1[0] * ceres::exp(-0.5 * u_1 * u_1);
        T gauss_2 = amp_gauss_2[0] * ceres::exp(-0.5 * u_2 * u_2);

        T erf = amp_erf[0] * ceres::erf(u_1);

        residual[0] = w * (y - (
            gauss_1 + gauss_2 + erf + lin[0] * x + off[0]
        ));

        return true;
    }
};


core::fitting::FitResult fit_erf_linear_2_gauss(
    const Eigen::VectorXd& x,
    const Eigen::VectorXd& y,
    const std::array<double, 9>& init
);


struct ErfLinear3GaussResidual {
    const double x;
    const double y;
    const double w;

    ErfLinear3GaussResidual(double x, double y, double w) : x(x), y(y), w(w) {}

    template <typename T>
    bool operator()(
        const T* const amp_gauss_1,
        const T* const x0_1,
        const T* const sig_1,
        const T* const amp_gauss_2,
        const T* const x0_2,
        const T* const sig_2,
        const T* const amp_gauss_3,
        const T* const x0_3,
        const T* const sig_3,
        const T* const amp_erf,
        const T* const lin,
        const T* const off,
        T* residual
    ) const {
        T u_1 = (x - x0_1[0]) / sig_1[0];
        T u_2 = (x - x0_2[0]) / sig_2[0];
        T u_3 = (x - x0_3[0]) / sig_3[0];

        T gauss_1 = amp_gauss_1[0] * ceres::exp(-0.5 * u_1 * u_1);
        T gauss_2 = amp_gauss_2[0] * ceres::exp(-0.5 * u_2 * u_2);
        T gauss_3 = amp_gauss_3[0] * ceres::exp(-0.5 * u_3 * u_3);

        T erf = amp_erf[0] * ceres::erf(u_1);

        residual[0] = w * (y - (
            gauss_1 + gauss_2 + gauss_3 + erf + lin[0] * x + off[0]
        ));

        return true;
    }
};


core::fitting::FitResult fit_erf_linear_3_gauss(
    const Eigen::VectorXd& x,
    const Eigen::VectorXd& y,
    const std::array<double, 12>& init
);


} // namespace compositron::dbs::fitting
