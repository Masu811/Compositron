#pragma once

#include <optional>
#include <string>

#include "compositron/core/utils.hpp"
#include "compositron/pals/fitting.hpp"
#include "compositron/pals/model.hpp"


namespace compositron::pals {


class PALSpectrum {
public:
    core::utils::Spectrum spectrum;
    core::utils::TimingDetectorPair detpair;
    uint64_t counts;
    double dcounts;
    std::optional<std::string> name;
    std::optional<Eigen::VectorXd> peak;
    std::optional<std::array<size_t, 2>> peak_bnds;
    std::optional<double> peak_counts;
    std::optional<double> dpeak_counts;
    std::optional<double> peak_center;
    std::optional<double> dpeak_center;
    std::optional<double> peak_fwhm;
    std::optional<double> dpeak_fwhm;
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

    Eigen::VectorXd get_peak_energies();

    void extract_peak(size_t fit_start_idx, size_t fit_end_idx);

    void get_peak_center();

    void fit(
        const model::LifetimeModel& model, size_t fit_start_idx, size_t fit_end_idx
    );

    std::string fit_report();
};


} // namespace compositron::pals
