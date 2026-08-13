#include "compositron/core/fitting.hpp"
#include <cmath>


namespace compositron::core::fitting {


double compute_red_chi2(const ceres::Solver::Summary& summary) {
    double chi_square = 2.0 * summary.final_cost;
    int num_residuals = summary.num_residuals;
    int num_parameters = summary.num_effective_parameters;
    return chi_square / (num_residuals - num_parameters);
}


std::optional<Eigen::MatrixXd> compute_covariance(
    ceres::Problem& problem, const std::vector<double>& params
) {
    ceres::Covariance::Options cov_options;
    ceres::Covariance covariance(cov_options);

    std::vector<std::pair<const double*, const double*>> covariance_blocks;

    size_t n = params.size();

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i; j < n; ++j) {
            covariance_blocks.push_back(std::make_pair(&params[i], &params[j]));
        }
    }

    bool success = covariance.Compute(covariance_blocks, &problem);

    if (!success) {
        return std::nullopt;
    }

    Eigen::MatrixXd cov(n, n);

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i; j < n; ++j) {
            covariance.GetCovarianceBlock(&params[i], &params[j], &cov(i, j));
            covariance.GetCovarianceBlock(&params[i], &params[j], &cov(j, i));
        }
    }

    return cov;
}


FitResult fit(ceres::Problem& problem, std::vector<double>& opt) {
    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    ceres::Solver::Summary summary;
    Solve(options, &problem, &summary);

    core::fitting::FitResult fit_result;

    fit_result.fit_status = summary.termination_type;

    if (fit_result.fit_status != ceres::CONVERGENCE) {
        return fit_result;
    }

    fit_result.red_chi_2 = core::fitting::compute_red_chi2(summary);

    auto cov = core::fitting::compute_covariance(problem, opt);

    if (cov.has_value()) {
        auto c = cov.value();
        for (size_t i = 0; i < opt.size(); ++i) {
            double val = opt[i];
            double err = std::sqrt(c(i, i));
            fit_result.opt.push_back(core::fitting::SimpleFitParam{ val, err });
        }

        fit_result.cov = cov;
    } else {
        for (size_t i = 0; i < opt.size(); ++i) {
            double val = opt[i];
            fit_result.opt.push_back(core::fitting::SimpleFitParam{ val, NAN });
        }
    }

    return fit_result;
}


} // namespace compositron::core::fitting
