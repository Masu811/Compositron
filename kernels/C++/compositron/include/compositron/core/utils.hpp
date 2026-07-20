#pragma once

#include <algorithm>
#include <cstdint>
#include <memory>
#include <vector>
#include <tuple>
#include <cassert>

namespace compositron::core::utils {

// class Spectrum

class SpectrumBase { };

template<typename T>
class Spectrum : public SpectrumBase {
private:
    std::vector<T> _data;

public:
    Spectrum() = delete;

    Spectrum(std::vector<T>&& data) {
        _data = data;
    }

    const std::vector<T>& data() const {
        return _data;
    }

    size_t size() const{
        return _data.size();
    }

    T& operator[](size_t i) {
        return _data[i];
    }
};

typedef std::unique_ptr<SpectrumBase> SpectrumPtr;

template<typename T>
SpectrumPtr toMinTypeSpectrum(const std::vector<T>& data) {
    auto max = std::max_element(data.begin(), data.end());

    if (*max <= UINT8_MAX) {
        return std::make_unique<Spectrum<uint8_t>>(
            std::vector<uint8_t>(data.begin(), data.end())
        );
    } else if (*max <= UINT16_MAX) {
        return std::make_unique<Spectrum<uint16_t>>(
            std::vector<uint16_t>(data.begin(), data.end())
        );
    } else if (*max <= UINT32_MAX) {
        return std::make_unique<Spectrum<uint32_t>>(
            std::vector<uint32_t>(data.begin(), data.end())
        );
    } else if (*max <= UINT64_MAX) {
        return std::make_unique<Spectrum<uint64_t>>(
            std::vector<uint64_t>(data.begin(), data.end())
        );
    }

    throw std::range_error("Greater than 64 bit uint values are not supported");
}

// class Spectrum 2D

class Spectrum2DBase { };

template<typename T>
class Spectrum2D : public Spectrum2DBase {
private:
    std::tuple<size_t, size_t> _shape;
    std::vector<T> _data;

public:
    Spectrum2D() = delete;

    Spectrum2D(std::vector<T>&& data, std::tuple<size_t, size_t> shape) {
        auto [n_rows, n_cols] = shape;
        assert(data.size() == n_rows * n_cols);
        _shape = shape;
        _data = data;
    }

    const std::vector<T>& data() const {
        return _data;
    }

    std::tuple<size_t, size_t> shape() const {
        return _shape;
    }

    T& operator()(size_t i, size_t j) {
        auto [n_rows, n_cols] = _shape;
        return _data[i * n_cols + j];
    }
};

typedef std::unique_ptr<Spectrum2DBase> Spectrum2DPtr;

template<typename T>
Spectrum2DPtr toMinTypeSpectrum2D(
    const std::vector<T>& data, std::tuple<size_t, size_t> shape
) {
    auto max = std::max_element(data.begin(), data.end());

    if (*max <= UINT8_MAX) {
        return std::make_unique<Spectrum2D<uint8_t>>(
            std::vector<uint8_t>(data.begin(), data.end()), shape
        );
    } else if (*max <= UINT16_MAX) {
        return std::make_unique<Spectrum2D<uint16_t>>(
            std::vector<uint16_t>(data.begin(), data.end()), shape
        );
    } else if (*max <= UINT32_MAX) {
        return std::make_unique<Spectrum2D<uint32_t>>(
            std::vector<uint32_t>(data.begin(), data.end()), shape
        );
    } else if (*max <= UINT64_MAX) {
        return std::make_unique<Spectrum2D<uint64_t>>(
            std::vector<uint64_t>(data.begin(), data.end()), shape
        );
    }

    throw std::range_error("Greater than 64 bit uint values are not supported");
}

} // namespace compositron::core::utils
