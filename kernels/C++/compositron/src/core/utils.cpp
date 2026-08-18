#include "compositron/core/utils.hpp"

#include "compositron/constants.hpp"
#include <limits>
#include <stdexcept>


namespace compositron::core::utils {


using Eigen::Dynamic;
using Eigen::RowMajor;
using Eigen::ColMajor;


// class LinearCalibration

double LinearCalibration::from_index(size_t index) const {
    return scale * index + offset;
}

double LinearCalibration::from_index(double index) const {
    return scale * index + offset;
}

size_t LinearCalibration::to_index(double value) const {
    return std::max((value - offset) / scale, 0.);
}

size_t LinearCalibration::to_index_rounded(double value) const {
    return std::max(std::round((value - offset) / scale), 0.);
}

double LinearCalibration::to_index_double(double value) const {
    return std::max((value - offset) / scale, 0.);
}


// class Unit

double Unit::to_kev(double eres) const {
    switch (type) {
        case UnitType::KEV:
            return value;
        case UnitType::M0C:
            return value * KEV_PER_M0C;
        case UnitType::ERES:
            if (std::isnan(eres)) {
                throw std::runtime_error(
                    "Could not convert units of eres to keV due to missing eres"
                );
            }
            return value * eres;
        default:
            return std::numeric_limits<double>::quiet_NaN();
    }
}


// class Spectrum

std::size_t Spectrum::size() const {
    return std::visit(
        [](const auto& data) {
            return data.size();
        },
        data
    );
}

// class Spectrum 2D

std::size_t Spectrum2D::size() const {
    return std::visit(
        [](const auto& data) {
            return data.size();
        },
        data
    );
}

std::size_t Spectrum2D::rows() const {
    return std::visit(
        [](const auto& data) {
            return data.rows();
        },
        data
    );
}

std::size_t Spectrum2D::cols() const {
    return std::visit(
        [](const auto& data) {
            return data.cols();
        },
        data
    );
}


} // namespace compositron::core::utils
