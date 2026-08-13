#include "compositron/dbs/fitting.hpp"
#include "compositron/core/fitting.hpp"

#include <cmath>


namespace compositron::dbs::fitting {


double gauss(double x, double amp, double x0, double sig) {
    double u = (x - x0) / sig;
    return amp * std::exp(-0.5 * u * u);
}


core::fitting::FitResult fit_gaussian(
    const Eigen::VectorXd& x,
    const Eigen::VectorXd& y,
    const std::array<double, 3>& init
) {
    ceres::Problem problem;

    std::vector<double> opt(init.begin(), init.end());

    for (size_t i = 0; i < x.size(); ++i) {
        auto* cost_function =
            new ceres::AutoDiffCostFunction<GaussResidual, 1, 1, 1, 1>(
                x[i], y[i], 1 / std::sqrt(std::max(y[i], 1.))
            );

        problem.AddResidualBlock(cost_function, nullptr, &opt[0], &opt[1], &opt[2]);
    }

    return core::fitting::fit(problem, opt);
}


double erf_linear_background(
    double x, double amp, double x0, double sig, double lin, double off
) {
    return amp * std::erf((x - x0) / sig) + lin * x + off;
}


core::fitting::FitResult fit_erf_linear_1_gauss(
    const Eigen::VectorXd& x,
    const Eigen::VectorXd& y,
    const std::array<double, 6>& init
) {
    ceres::Problem problem;

    std::vector<double> opt(init.begin(), init.end());

    for (size_t i = 0; i < x.size(); ++i) {
        auto* cost_function =
            new ceres::AutoDiffCostFunction<ErfLinear1GaussResidual, 1, 1, 1, 1, 1, 1, 1>(
                x[i], y[i], 1 / std::sqrt(std::max(y[i], 1.))
            );

        problem.AddResidualBlock(
            cost_function, nullptr, &opt[0], &opt[1], &opt[2], &opt[3], &opt[4], &opt[5]
        );
    }

    return core::fitting::fit(problem, opt);
}


core::fitting::FitResult fit_erf_linear_2_gauss(
    const Eigen::VectorXd& x,
    const Eigen::VectorXd& y,
    const std::array<double, 9>& init
) {
    ceres::Problem problem;

    std::vector<double> opt(init.begin(), init.end());

    for (size_t i = 0; i < x.size(); ++i) {
        auto* cost_function =
            new ceres::AutoDiffCostFunction<ErfLinear2GaussResidual, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1>(
                x[i], y[i], 1 / std::sqrt(std::max(y[i], 1.))
            );

        problem.AddResidualBlock(
            cost_function, nullptr, &opt[0], &opt[1], &opt[2], &opt[3], &opt[4], &opt[5], &opt[6], &opt[7], &opt[8]
        );
    }

    return core::fitting::fit(problem, opt);
}


core::fitting::FitResult fit_erf_linear_3_gauss(
    const Eigen::VectorXd& x,
    const Eigen::VectorXd& y,
    const std::array<double, 12>& init
) {
    ceres::Problem problem;

    std::vector<double> opt(init.begin(), init.end());

    for (size_t i = 0; i < x.size(); ++i) {
        auto* cost_function =
            new ceres::AutoDiffCostFunction<ErfLinear3GaussResidual, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1>(
                x[i], y[i], 1 / std::sqrt(std::max(y[i], 1.))
            );

        problem.AddResidualBlock(
            cost_function, nullptr, &opt[0], &opt[1], &opt[2], &opt[3], &opt[4], &opt[5], &opt[6], &opt[7], &opt[8], &opt[9], &opt[10], &opt[11]
        );
    }

    return core::fitting::fit(problem, opt);
}


} // namespace compositron::dbs::fitting
