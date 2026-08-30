#pragma once

#include <Eigen/Core>
#include <optional>
#include <string>
#include <unordered_map>

#include "compositron/core/fitting.hpp"
#include "compositron/core/utils.hpp"
#include "compositron/pals/fitting.hpp"
#include "compositron/pals/model.hpp"


namespace compositron::pals {


class PALSpectrum {
public:
    core::utils::Spectrum spectrum;
    core::utils::TimingDetectorPair detpair;
    uint64_t counts;
    std::unordered_map<std::string, core::fitting::SimpleFitParam> peak_params;
    std::optional<fitting::FitResult> fit_result;

    PALSpectrum() = delete;

    PALSpectrum(PALSpectrum& other) = delete;

    PALSpectrum(PALSpectrum&& other) = default;

    PALSpectrum(
        core::utils::Spectrum&& spectrum,
        core::utils::TimingDetectorPair detpair
    );

    template<typename T>
    PALSpectrum(
        const std::vector<T>& spectrum, core::utils::TimingDetectorPair detpair
    ) : PALSpectrum(core::utils::Spectrum(spectrum), detpair) {}

    Eigen::VectorXd get_peak_energies(size_t fit_start_idx, size_t fit_end_idx);

    Eigen::VectorXd extract_peak(size_t fit_start_idx, size_t fit_end_idx);

    double get_peak_center();

    void fit(
        const model::LifetimeModel& model, size_t fit_start_idx, size_t fit_end_idx
    );

    std::string fit_report();
};


} // namespace compositron::pals
