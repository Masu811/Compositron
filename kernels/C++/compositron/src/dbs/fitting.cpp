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
    return core::fitting::fit_generic<GaussResidual, 3>(x, y, init);
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
    return core::fitting::fit_generic<ErfLinear1GaussResidual, 6>(x, y, init);
}


core::fitting::FitResult fit_erf_linear_2_gauss(
    const Eigen::VectorXd& x,
    const Eigen::VectorXd& y,
    const std::array<double, 9>& init
) {
    return core::fitting::fit_generic<ErfLinear2GaussResidual, 9>(x, y, init);
}


core::fitting::FitResult fit_erf_linear_3_gauss(
    const Eigen::VectorXd& x,
    const Eigen::VectorXd& y,
    const std::array<double, 12>& init
) {
    return core::fitting::fit_generic<ErfLinear3GaussResidual, 12>(x, y, init);
}


} // namespace compositron::dbs::fitting
