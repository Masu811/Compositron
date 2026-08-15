#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <variant>
#include <vector>
#include <cassert>

#include <Eigen/Core>
#include <Eigen/Dense>


namespace compositron::core::utils {


// String trimming (https://stackoverflow.com/questions/216823/how-can-i-trim-a-stdstring)

// Trim from the start (in place)
inline void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
        return !std::isspace(ch);
    }));
}

// Trim from the end (in place)
inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

// Trim from both ends (in place)
inline void trim(std::string &s) {
    rtrim(s);
    ltrim(s);
}


class LinearCalibration {
public:
    double offset;
    double scale;

    double from_index(size_t index);
    double from_index(double index);
    size_t to_index(double value);
    size_t to_index_rounded(double value);
    double to_index_double(double value);
};


// DBS detector
struct EnergyDetector {
    std::string name;
    LinearCalibration ecal;
    std::optional<LinearCalibration> corrected_ecal;
    double eres;
};


// CDBS detector
struct EnergyDetectorPair {
    std::string name;
    EnergyDetector first_det;
    EnergyDetector second_det;
    double eres;
};


// PALS detector
struct TimingDetectorPair {
    std::string name;
    LinearCalibration tcal;
    std::optional<LinearCalibration> corrected_tcal;
    double tres;
};


enum UnitType {
    KEV,
    M0C,
    ERES,
};


class Unit {
public:
    UnitType type;
    double value;

    Unit() = default;

    Unit(double value, UnitType type) : type(type), value(value) {}

    double to_kev(double eres);
};


enum EcalCorrectionOrder {
    ZEROTH,
    FIRST,
    NONE,
};


using SpectrumData = std::variant<
    Eigen::VectorX<std::uint8_t>,
    Eigen::VectorX<std::uint16_t>,
    Eigen::VectorX<std::uint32_t>,
    Eigen::VectorX<std::uint64_t>
>;


class Spectrum {
public:
    SpectrumData data;

    Spectrum() = delete;

    template <typename T>
    Spectrum(std::vector<T>& vec_data) {
        auto max = std::max_element(vec_data.begin(), vec_data.end());

        if (*max <= UINT8_MAX) {
            data = Eigen::VectorX<uint8_t>(
                Eigen::Map<Eigen::VectorX<std::uint8_t>>(
                    std::vector<uint8_t>(vec_data.begin(), vec_data.end()).data(),
                    vec_data.size()
                )
            );
        } else if (*max <= UINT16_MAX) {
            data = Eigen::VectorX<uint16_t>(
                Eigen::Map<Eigen::VectorX<std::uint16_t>>(
                    std::vector<uint16_t>(vec_data.begin(), vec_data.end()).data(),
                    vec_data.size()
                )
            );
        } else if (*max <= UINT32_MAX) {
            data = Eigen::VectorX<uint32_t>(
                Eigen::Map<Eigen::VectorX<std::uint32_t>>(
                    std::vector<uint32_t>(vec_data.begin(), vec_data.end()).data(),
                    vec_data.size()
                )
            );
        } else if (*max <= UINT64_MAX) {
            data = Eigen::VectorX<uint64_t>(
                Eigen::Map<Eigen::VectorX<std::uint64_t>>(
                    std::vector<uint64_t>(vec_data.begin(), vec_data.end()).data(),
                    vec_data.size()
                )
            );
        } else {
            throw std::range_error(
                "Greater than 64 bit uint values are not supported"
            );
        }
    }

    std::size_t size() const;
};


enum StorageOrder {
    COLMAJ,
    ROWMAJ
};


using Spectrum2DData = std::variant<
    Eigen::MatrixX<std::uint8_t>,
    Eigen::MatrixX<std::uint16_t>,
    Eigen::MatrixX<std::uint32_t>,
    Eigen::MatrixX<std::uint64_t>
>;


class Spectrum2D {
private:
    Spectrum2DData data;

private:
    template <typename T>
    Eigen::MatrixX<T> convertSpectrum(
        std::vector<T>&& vec_data,
        size_t nrows,
        size_t ncols,
        StorageOrder order
    ) {
        using Eigen::Dynamic;
        using Eigen::ColMajor;
        using Eigen::RowMajor;

        if (order == StorageOrder::ROWMAJ) {
            return Eigen::Matrix<T, Dynamic, Dynamic, ColMajor>(
                Eigen::Map<Eigen::Matrix<T, Dynamic, Dynamic, RowMajor>>(
                    vec_data.data(), nrows, ncols
                )
            );
        }

        return Eigen::Matrix<T, Dynamic, Dynamic, ColMajor>(
            Eigen::Map<Eigen::Matrix<T, Dynamic, Dynamic, ColMajor>>(
                vec_data.data(), nrows, ncols
            )
        );
    }

public:
    Spectrum2D() = delete;

    template <typename T>
    Spectrum2D(
        std::vector<T>& vec_data,
        size_t nrows,
        size_t ncols,
        StorageOrder order
    ) {
        auto max = std::max_element(vec_data.begin(), vec_data.end());

        if (*max <= UINT8_MAX) {
            data = convertSpectrum<std::uint8_t>(
                std::vector<uint8_t>(vec_data.begin(), vec_data.end()),
                nrows, ncols, order
            );
        } else if (*max <= UINT16_MAX) {
            data = convertSpectrum<std::uint16_t>(
                std::vector<uint16_t>(vec_data.begin(), vec_data.end()),
                nrows, ncols, order
            );
        } else if (*max <= UINT32_MAX) {
            data = convertSpectrum<std::uint32_t>(
                std::vector<uint32_t>(vec_data.begin(), vec_data.end()),
                nrows, ncols, order
            );
        } else if (*max <= UINT64_MAX) {
            data = convertSpectrum<std::uint64_t>(
                std::vector<uint64_t>(vec_data.begin(), vec_data.end()),
                nrows, ncols, order
            );
        } else {
            throw std::range_error(
                "Greater than 64 bit uint values are not supported"
            );
        }
    }

    std::size_t size() const;
};


} // namespace compositron::core::utils
