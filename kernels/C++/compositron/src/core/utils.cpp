#include "compositron/core/utils.hpp"

#include "compositron/constants.hpp"


namespace compositron::core::utils {


using Eigen::Dynamic;
using Eigen::RowMajor;
using Eigen::ColMajor;


// class LinearCalibration

double LinearCalibration::from_index(size_t index) {
    return scale * index + offset;
}

double LinearCalibration::from_index(double index) {
    return scale * index + offset;
}

size_t LinearCalibration::to_index(double value) {
    return (value - offset) / scale;
}

size_t LinearCalibration::to_index_rounded(double value) {
    return std::round((value - offset) / scale);
}

double LinearCalibration::to_float_index(double value) {
    return (value - offset) / scale;
}


// class Unit

std::optional<double> Unit::to_kev(std::optional<double> eres) {
    switch (type) {
        case UnitType::KEV:
            return value;
            break;
        case UnitType::M0C:
            return value * KEV_PER_M0C;
            break;
        case UnitType::ERES:
            if (eres.has_value()) {
                return value * eres.value();
            } else {
                return std::nullopt;
            }
            break;
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


using Spectrum2DData = std::variant<
    Eigen::MatrixX<std::uint8_t>,
    Eigen::MatrixX<std::uint16_t>,
    Eigen::MatrixX<std::uint32_t>,
    Eigen::MatrixX<std::uint64_t>
>;


std::size_t Spectrum2D::size() const {
    return std::visit(
        [](const auto& data) {
            return data.size();
        },
        data
    );
}


} // namespace compositron::core::utils
