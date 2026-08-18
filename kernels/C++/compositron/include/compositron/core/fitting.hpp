#pragma once

#include <ceres/loss_function.h>
#include <ceres/problem.h>
#include <cstddef>
#include <optional>
#include <vector>

#include <ceres/ceres.h>
#include <ceres/solver.h>
#include <ceres/types.h>


namespace compositron::core::fitting {


struct SimpleFitParam {
    double val;
    double err;
};


struct FitParam {
    double val;
    double err = NAN;
    double min = -std::numeric_limits<double>::infinity();
    double max = std::numeric_limits<double>::infinity();
    bool vary = true;
    bool default_initialized = false;

    FitParam() = default;
    FitParam(double val) : val(val) {}
    FitParam(double val, bool vary) : val(val), vary(vary) {}
    FitParam(double val, double min) : val(val), min(min) {}
    FitParam(double val, bool vary, double min) : val(val), vary(vary), min(min) {}
    FitParam(double val, double min, double max) : val(val), min(min), max(max) {}
    FitParam(
        double val, bool vary, double min, double max
    ) : val(val), vary(vary), min(min), max(max) {}
};


struct FitResult {
    std::vector<SimpleFitParam> opt;
    ceres::Solver::Summary summary;
    ceres::TerminationType fit_status;
    size_t n_eval;
    std::optional<double> red_chi_2;
    std::optional<Eigen::MatrixXd> cov;
};


double compute_red_chi2(const ceres::Solver::Summary& summary);


std::optional<Eigen::MatrixXd> compute_covariance_single_block(
    ceres::Problem& problem, const std::vector<double>& params
);


std::optional<Eigen::MatrixXd> compute_covariance_multi_blocks(
    ceres::Problem& problem, const std::vector<double>& params
);


FitResult fit(ceres::Problem& problem, std::vector<double>& opt);


template <typename CostFunctor, size_t N>
core::fitting::FitResult fit_generic(
    const Eigen::VectorXd& x,
    const Eigen::VectorXd& y,
    const std::array<double, N>& init
) {
    ceres::Problem problem;

    std::vector<double> opt(init.begin(), init.end());

    for (size_t i = 0; i < x.size(); ++i) {
        auto* cost_function =
            new ceres::AutoDiffCostFunction<CostFunctor, 1, N>(
                x[i], y[i], 1 / std::sqrt(std::max(y[i], 1.))
            );

        problem.AddResidualBlock(cost_function, nullptr, opt.data());
    }

    return core::fitting::fit(problem, opt);
}


template <typename CostFunctor, size_t N>
core::fitting::FitResult fit_generic_2d(
    const Eigen::VectorXd& x,
    const Eigen::VectorXd& y,
    const Eigen::MatrixXd& z,
    const std::array<double, N>& init
) {
    ceres::Problem problem;

    std::vector<double> opt(init.begin(), init.end());

    for (size_t i = 0; i < y.size(); ++i) {
        for (size_t j = 0; j < x.size(); ++j) {
            auto* cost_function =
                new ceres::AutoDiffCostFunction<CostFunctor, 1, N>(
                    x[j], y[i], z(i, j), 1 / std::sqrt(std::max(z(i, j), 1.))
                );

            problem.AddResidualBlock(cost_function, nullptr, opt.data());
        }
    }

    return core::fitting::fit(problem, opt);
}


} // namespace compositron::core::fitting
