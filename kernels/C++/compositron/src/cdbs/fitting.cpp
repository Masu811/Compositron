#include "compositron/cdbs/fitting.hpp"
#include "compositron/core/fitting.hpp"

#include <cmath>


namespace compositron::cdbs::fitting {


double gauss2d(
    double x, double y, double amp, double x0, double y0,
    double sig_x, double sig_y, double phi
) {
    double cos_phi = std::cos(phi);
    double sin_phi = std::sin(phi);

    double x_rotated = (x - x0) * cos_phi - (y - y0) * sin_phi;
    double y_rotated = (x - x0) * sin_phi + (y - y0) * cos_phi;

    double u = x_rotated / sig_x;
    double v = y_rotated / sig_y;

    return amp * std::exp(-0.5 * u * u - 0.5 * v * v);
}


core::fitting::FitResult fit_gauss2d(
    const Eigen::VectorXd& x,
    const Eigen::VectorXd& y,
    const Eigen::MatrixXd& z,
    const std::array<double, 6>& init
) {
    return core::fitting::fit_generic_2d<Gauss2DResidual, 6>(x, y, z, init);
}


} // namespace compositron::cdbs::fitting
