#pragma once

#include "compositron/core/fitting.hpp"


namespace compositron::cdbs::fitting {


double gauss2d(
    double x, double y, double amp, double x0, double y0,
    double sig_x, double sig_y, double phi
);


struct Gauss2DResidual {
    const double x;
    const double y;
    const double z;
    const double w;

    Gauss2DResidual(
        double x, double y, double z, double w
    ) : x(x), y(y), z(z), w(w) {}

    template <typename T>
    bool operator()(const T* const params, T* residual) const {
        T amp = params[0];
        T x0 = params[1];
        T y0 = params[2];
        T sig_x = params[3];
        T sig_y = params[4];
        T phi = params[5];

        T cos_phi = ceres::cos(phi);
        T sin_phi = ceres::sin(phi);

        T x_rotated = (x - x0) * cos_phi - (y - y0) * sin_phi;
        T y_rotated = (x - x0) * sin_phi + (y - y0) * cos_phi;

        T u = x_rotated / sig_x;
        T v = y_rotated / sig_y;

        // Not weighted on purpose!
        residual[0] = z - amp * ceres::exp(-0.5 * u * u - 0.5 * v * v);

        return true;
    }
};


std::unordered_map<std::string, core::fitting::SimpleFitParam> fit_gauss2d(
    const Eigen::VectorXd& x,
    const Eigen::VectorXd& y,
    const Eigen::MatrixXd& z,
    const std::array<double, 6>& init
);


} // namespace compositron::cdbs::fitting
