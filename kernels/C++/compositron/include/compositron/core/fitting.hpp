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


struct FitResult {
    std::vector<SimpleFitParam> opt;
    ceres::Solver::Summary summary;
    ceres::TerminationType fit_status;
    size_t n_eval;
    std::optional<double> red_chi_2;
    std::optional<Eigen::MatrixXd> cov;
};


double compute_red_chi2(const ceres::Solver::Summary& summary);


std::optional<Eigen::MatrixXd> compute_covariance(
    ceres::Problem& problem, const std::vector<double>& params
);


FitResult fit(ceres::Problem& problem, std::vector<double>& opt);


} // namespace compositron::core::fitting
