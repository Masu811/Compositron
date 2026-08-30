#include "compositron/cdbs/fitting.hpp"
#include "compositron/core/fitting.hpp"

#include <cmath>
#include <unordered_map>


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


std::unordered_map<std::string, core::fitting::SimpleFitParam> fit_gauss2d(
    const Eigen::VectorXd& x,
    const Eigen::VectorXd& y,
    const Eigen::MatrixXd& z,
    const std::array<double, 6>& init
) {
    auto result = core::fitting::fit_generic_2d<Gauss2DResidual, 6>(x, y, z, init);

    if (result.fit_status != ceres::CONVERGENCE) {
        throw std::runtime_error("Fit failed");
    }

    std::unordered_map<std::string, core::fitting::SimpleFitParam> map;

    map["amp"] = result.opt[0];
    map["x0"] = result.opt[1];
    map["y0"] = result.opt[2];
    map["sig_x"] = result.opt[3];
    map["sig_y"] = result.opt[4];
    map["phi"] = result.opt[5];

    return map;
}


} // namespace compositron::cdbs::fitting
