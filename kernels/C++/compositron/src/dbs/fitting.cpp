#include "compositron/dbs/fitting.hpp"

#include <cmath>
#include <stdexcept>
#include <unordered_map>

#include <ceres/types.h>

#include "compositron/core/fitting.hpp"


namespace compositron::dbs::fitting {


double gauss(double x, double amp, double x0, double sig) {
    double u = (x - x0) / sig;
    return amp * std::exp(-0.5 * u * u);
}


std::unordered_map<std::string, core::fitting::SimpleFitParam> fit_gauss(
    const Eigen::VectorXd& x,
    const Eigen::VectorXd& y,
    const std::array<double, 3>& init
) {
    auto result = core::fitting::fit_generic<GaussResidual, 3>(x, y, init);

    if (result.fit_status != ceres::CONVERGENCE) {
        throw std::runtime_error("Gauss fit failed");
    }

    std::unordered_map<std::string, core::fitting::SimpleFitParam> map{{
        {"amp_1", result.opt[0]},
        {"x0_1", result.opt[1]},
        {"sig_1", result.opt[2]},
    }};

    return map;
}


double erf_linear_background(
    double x, double amp, double x0, double sig, double lin, double off
) {
    return amp * std::erfc((x - x0) / sig) + lin * x + off;
}


std::unordered_map<std::string, core::fitting::SimpleFitParam> fit_erf_linear_1_gauss(
    const Eigen::VectorXd& x,
    const Eigen::VectorXd& y,
    const std::array<double, 6>& init
) {
    auto result = core::fitting::fit_generic<ErfLinear1GaussResidual, 6>(x, y, init);

    if (result.fit_status != ceres::CONVERGENCE) {
        throw std::runtime_error("Gauss fit failed");
    }

    std::unordered_map<std::string, core::fitting::SimpleFitParam> map{{
        {"amp_1", result.opt[0]},
        {"x0_1", result.opt[1]},
        {"sig_1", result.opt[2]},
        {"erf_amp", result.opt[3]},
        {"lin", result.opt[4]},
        {"const", result.opt[5]},
    }};

    return map;
}


std::unordered_map<std::string, core::fitting::SimpleFitParam> fit_erf_linear_2_gauss(
    const Eigen::VectorXd& x,
    const Eigen::VectorXd& y,
    const std::array<double, 9>& init
) {
    auto result = core::fitting::fit_generic<ErfLinear2GaussResidual, 9>(x, y, init);

    if (result.fit_status != ceres::CONVERGENCE) {
        throw std::runtime_error("Gauss fit failed");
    }

    std::unordered_map<std::string, core::fitting::SimpleFitParam> map{{
        {"amp_1", result.opt[0]},
        {"x0_1", result.opt[1]},
        {"sig_1", result.opt[2]},
        {"amp_2", result.opt[3]},
        {"x0_2", result.opt[4]},
        {"sig_2", result.opt[5]},
        {"erf_amp", result.opt[6]},
        {"lin", result.opt[7]},
        {"const", result.opt[8]},
    }};

    return map;
}


std::unordered_map<std::string, core::fitting::SimpleFitParam> fit_erf_linear_3_gauss(
    const Eigen::VectorXd& x,
    const Eigen::VectorXd& y,
    const std::array<double, 12>& init
) {
    auto result = core::fitting::fit_generic<ErfLinear3GaussResidual, 12>(x, y, init);

    if (result.fit_status != ceres::CONVERGENCE) {
        throw std::runtime_error("Gauss fit failed");
    }

    std::unordered_map<std::string, core::fitting::SimpleFitParam> map{{
        {"amp_1", result.opt[0]},
        {"x0_1", result.opt[1]},
        {"sig_1", result.opt[2]},
        {"amp_2", result.opt[3]},
        {"x0_2", result.opt[4]},
        {"sig_2", result.opt[5]},
        {"amp_3", result.opt[6]},
        {"x0_3", result.opt[7]},
        {"sig_3", result.opt[8]},
        {"erf_amp", result.opt[9]},
        {"lin", result.opt[10]},
        {"const", result.opt[11]},
    }};

    return map;
}


} // namespace compositron::dbs::fitting
