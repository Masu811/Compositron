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
    bool operator()(const T* const params, T* residual) const {
        T amp = params[0];
        T x0 = params[1];
        T sig = params[2];

        T u = (x - x0) / sig;

        residual[0] = w * (y - amp * ceres::exp(-0.5 * u * u));

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
    bool operator()(const T* const params, T* residual) const {
        T amp_gauss_1 = params[0];
        T x0_1 = params[1];
        T sig_1 = params[2];
        T amp_erf = params[3];
        T lin = params[4];
        T off = params[5];

        T u_1 = (x - x0_1) / sig_1;

        T gauss_1 = amp_gauss_1 * ceres::exp(-0.5 * u_1 * u_1);

        T erf = amp_erf * ceres::erf(u_1);

        residual[0] = w * (y - (gauss_1 + erf + lin * x + off));

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
    bool operator()(const T* const params, T* residual) const {
        T amp_gauss_1 = params[0];
        T x0_1 = params[1];
        T sig_1 = params[2];
        T amp_gauss_2 = params[3];
        T x0_2 = params[4];
        T sig_2 = params[5];
        T amp_erf = params[6];
        T lin = params[7];
        T off = params[8];

        T u_1 = (x - x0_1) / sig_1;
        T u_2 = (x - x0_2) / sig_2;

        T gauss_1 = amp_gauss_1 * ceres::exp(-0.5 * u_1 * u_1);
        T gauss_2 = amp_gauss_2 * ceres::exp(-0.5 * u_2 * u_2);

        T erf = amp_erf * ceres::erf(u_1);

        residual[0] = w * (y - (gauss_1 + gauss_2 + erf + lin * x + off));

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
    bool operator()(const T* const params, T* residual) const {
        T amp_gauss_1 = params[0];
        T x0_1 = params[1];
        T sig_1 = params[2];
        T amp_gauss_2 = params[3];
        T x0_2 = params[4];
        T sig_2 = params[5];
        T amp_gauss_3 = params[6];
        T x0_3 = params[7];
        T sig_3 = params[8];
        T amp_erf = params[9];
        T lin = params[10];
        T off = params[11];

        T u_1 = (x - x0_1) / sig_1;
        T u_2 = (x - x0_2) / sig_2;
        T u_3 = (x - x0_3) / sig_3;

        T gauss_1 = amp_gauss_1 * ceres::exp(-0.5 * u_1 * u_1);
        T gauss_2 = amp_gauss_2 * ceres::exp(-0.5 * u_2 * u_2);
        T gauss_3 = amp_gauss_3 * ceres::exp(-0.5 * u_3 * u_3);

        T erf = amp_erf * ceres::erf(u_1);

        residual[0] = w * (y - (
            gauss_1 + gauss_2 + gauss_3 + erf + lin * x + off
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
