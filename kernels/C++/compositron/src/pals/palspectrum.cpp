#include "compositron/pals/palspectrum.hpp"

#include <Eigen/Core>
#include <algorithm>
#include <ranges>
#include <stdexcept>
#include <variant>

#include <ceres/types.h>

#include "compositron/constants.hpp"
#include "compositron/core/fitting.hpp"
#include "compositron/dbs/fitting.hpp"
#include "compositron/pals/fitting.hpp"


namespace compositron::pals {


// class PALSpectrum

PALSpectrum::PALSpectrum(
    core::utils::Spectrum&& spectrum,
    core::utils::TimingDetectorPair detpair
) : spectrum(spectrum), detpair(detpair) {
    counts = std::visit(
        [](const auto& spectrum) {
            return std::accumulate(
                spectrum.begin(),
                spectrum.end(),
                0,
                [](std::uint64_t acc, auto& x) { return acc + x; }
            );
        },
        this->spectrum.data
    );
}


Eigen::VectorXd PALSpectrum::extract_peak(size_t fit_start_idx, size_t fit_end_idx) {
    size_t max_size = spectrum.size();

    if (fit_start_idx >= max_size || fit_end_idx >= max_size) {
        throw std::runtime_error("Fit range lies outside of spectrum");
    }

    return std::visit(
        [fit_start_idx, fit_end_idx](const auto& spectrum) {
            auto peak = spectrum(Eigen::seq(fit_start_idx, fit_end_idx));
            return Eigen::VectorX<double>(peak.template cast<double>());
        },
        spectrum.data
    );
}


Eigen::VectorXd PALSpectrum::get_peak_energies(
    size_t fit_start_idx, size_t fit_end_idx
) {
    auto tcal = detpair.corrected_tcal.value_or(detpair.tcal);

    size_t len = fit_end_idx - fit_start_idx + 1;

    auto linspace = Eigen::VectorXd::LinSpaced(len, fit_start_idx, fit_end_idx);

    return linspace.array() * tcal.scale + tcal.offset;
}


double PALSpectrum::get_peak_center() {
    auto result = std::visit(
        [](const auto& spectrum) {
            Eigen::Index argmax;
            double max = spectrum.maxCoeff(&argmax);

            int peak_center_window = 15;

            int left_idx = argmax - peak_center_window;
            int right_idx = argmax + peak_center_window;

            auto x = Eigen::VectorXd::LinSpaced(
                2 * peak_center_window, left_idx, right_idx
            );

            auto peak = spectrum(Eigen::seq(left_idx, right_idx));
            auto y = Eigen::VectorX<double>(peak.template cast<double>());

            return dbs::fitting::fit_gauss(
                x, y, {max, static_cast<double>(argmax), 100.}
            );
        },
        spectrum.data
    );

    auto peak_height = result["amp_1"];
    auto peak_center = result["x0_1"];
    auto peak_width = result["sig_1"];

    peak_params["center_idx"] = peak_center;
    peak_params["height"] = peak_height;
    peak_params["fwhm"] = core::fitting::SimpleFitParam{
        peak_width.val * FWHM_OVER_SIGMA,
        peak_width.err * FWHM_OVER_SIGMA,
    };

    return peak_center.val;
}


void PALSpectrum::fit(
    const model::LifetimeModel& model, size_t fit_start_idx, size_t fit_end_idx
) {
    auto peak_center = get_peak_center();

    auto x = get_peak_energies(fit_start_idx, fit_end_idx);
    auto y = extract_peak(fit_start_idx, fit_end_idx);

    fit_result = fitting::fit_lifetime_model(
        x, y, model, counts, detpair.tcal.from_index(peak_center)
    );
}


std::string print_params(
    std::vector<core::fitting::FitParam> params,
    std::vector<std::string> names
) {
    std::string text;
    int nrows = names.size();

    auto values = std::views::zip(params, names)
        | std::views::transform([](const auto pair) {
            auto [param, name] = pair;
            if (name.starts_with("int")) {
                return std::format("{:.2f}", 100 * param.val);
            } else if (std::abs(param.val) >= 1e6) {
                return std::format("{:.3e}", param.val);
            } else {
                return std::format("{:.2f}", param.val);
            }
        })
        | std::ranges::to<std::vector<std::string>>();

    auto abs_errors = std::views::zip(params, names)
        | std::views::transform([](const auto pair) {
            auto [param, name] = pair;
            if (name.starts_with("int")) {
                return std::format("{:.2f}", 100 * param.err);
            } else if (std::abs(param.err) >= 1e6) {
                return std::format("{:.3e}", param.err);
            } else {
                return std::format("{:.2f}", param.err);
            }
        })
        | std::ranges::to<std::vector<std::string>>();

    auto rel_errors = params
        | std::views::transform([](const auto& param) {
            return std::format("{:.2f} %", std::abs(param.err / param.val * 100));
        })
        | std::ranges::to<std::vector<std::string>>();

    auto mins = std::views::zip(params, names)
        | std::views::transform([](const auto pair) {
            auto [param, name] = pair;
            if (name.starts_with("int")) {
                return std::format("{:.1f}", 100 * param.min);
            } else {
                return std::format("{:.1f}", param.min);
            }
        })
        | std::ranges::to<std::vector<std::string>>();

    auto maxs = std::views::zip(params, names)
        | std::views::transform([](const auto pair) {
            auto [param, name] = pair;
            if (name.starts_with("int")) {
                return std::format("{:.1f}", 100 * param.max);
            } else {
                return std::format("{:.1f}", param.max);
            }
        })
        | std::ranges::to<std::vector<std::string>>();

    auto fixed = params
        | std::views::transform([](const auto& param) {
            return param.vary ? "" : "Fixed";
        })
        | std::ranges::to<std::vector<std::string>>();

    auto at_bnd = params
        | std::views::transform([](const auto& param) {
            if (std::abs(param.val - param.min) < 1e-3 || std::abs(param.val - param.max) < 1e-3) {
                return "At Boundary";
            } else {
                return "";
            }
        })
        | std::ranges::to<std::vector<std::string>>();

    std::vector<std::vector<std::string>> columns{
        names, values, abs_errors, rel_errors, mins, maxs, fixed, at_bnd
    };

    auto col_widths = columns
        | std::views::transform([](const auto& column) {
            return std::ranges::max(
                column
                    | std::views::transform([](const auto& row) {
                        return row.size();
                    })
            );
        })
        | std::ranges::to<std::vector<size_t>>();

    std::vector<std::string> headers{
        "Parameter", "Value", "Stderr (abs)", "Stderr (rel)", "Min", "Max",
        "Fixed", ""
    };

    auto head_widths = headers
        | std::views::transform([](const auto& header) {
            return header.size();
        });

    auto widths = std::views::zip(head_widths, col_widths)
        | std::views::transform([](const auto& pair) {
            auto [head_width, col_width] = pair;
            return std::max(head_width, col_width);
        });

    std::vector<std::string> formats{"<", ">", ">", ">", ">", ">", "<", "<"};

    for (const auto [i, head] : std::views::enumerate(headers)) {
        text += std::format(" {:^{}} ", head, widths[i]);
    }

    text += "\n";

    for (const auto wid : widths) {
        text += std::format(" {:-^{}} ", "", wid);
    }

    text += "\n";

    for (const auto row : std::views::iota(0, nrows)) {
        for (const auto& [col, column] : std::views::enumerate(columns)) {
            if (formats[col] == "<") {
                text += std::format(" {:<{}} ", column[row], widths[col]);
            } else {
                text += std::format(" {:>{}} ", column[row], widths[col]);
            }
        }
        text += "\n";
    }

    return text;
}


std::string PALSpectrum::fit_report() {
    if (!fit_result.has_value()) {
        throw std::runtime_error("No fit has been performed yet");
    }

    std::string text;

    auto& result = fit_result.value();
    auto& model = result.model;

    int n_l = model.lifetime_components.size();
    int n_r = model.resolution_components.size();

    std::string thick_hline(100, '#');
    std::string thin_hline(100, '-');

    thick_hline += "\n";
    thin_hline += "\n";

    text += thick_hline;
    text += std::format("{:^100}\n", "Fit report");
    text += thick_hline;
    text += "\n";
    text += std::format("Statistics: {}\n", counts);
    text += std::format("Number of Data Points: {}\n", result.n_dpoints);
    text += "\n";
    text += std::format("Lifetime Components:   {}\n", n_l);
    text += std::format("Resolution Components: {}\n", n_r);
    if (model.background_component.default_initialized) {
        text += "Background Components: 0\n";
    } else {
        text += "Background Components: 1\n";
    }
    text += "\n";
    text += "Fit Method: Trust-Region Levenberg-Marquardt\n";
    text += "Fit Status: ";
    switch (result.fit_status) {
        case ceres::CONVERGENCE:
            text += "Convergence\n";
            break;
        case ceres::NO_CONVERGENCE:
            text += "No Convergence\n";
            break;
        case ceres::FAILURE:
            text += "Failure\n";
            break;
        case ceres::USER_SUCCESS:
            text += "User Success\n";
            break;
        case ceres::USER_FAILURE:
            text += "User Failure\n";
            break;
    }
    text += std::format("Reason for Termination: {}\n", result.summary.message);
    text += std::format("Number of Function Evaluations: {}\n", result.n_eval);
    if (result.red_chi_2.has_value()) {
        text += std::format("Reduced Chi Squared: {}\n", result.red_chi_2.value());
    } else {
        text += "Reduced Chi Squared: NaN\n";
    }
    text += "\n";

    text += thin_hline;
    text += std::format("{:^100}\n", " Parameters ");
    text += thin_hline;

    text += "\n";
    text += std::format("{:=^100}\n", " General Parameters ");
    text += "\n";

    text += print_params(
        {model.scale_component, model.shift_component, model.background_component},
        {"N", "t0", "background"}
    );

    text += "\n";
    text += std::format("{:=^100}\n", " Lifetime Components ");
    text += "\n";

    double lt_norm = std::accumulate(
        model.lifetime_components.begin(),
        model.lifetime_components.end(),
        0.,
        [](double acc, const auto& param) { return acc + param.intensity.val; }
    );

    if (std::abs(lt_norm - 1) > 1e-3) {
        text += "!!! intensities don't add up to 100% !!!\n\n";
    }

    std::vector<core::fitting::FitParam> lt_params;
    std::vector<std::string> lt_names;

    for (const auto i : std::views::iota(0, n_l)) {
        lt_params.push_back(model.lifetime_components[i].lifetime);
        lt_names.push_back(std::format("lifetime_{}", i+1));
    }
    for (const auto i : std::views::iota(0, n_l)) {
        lt_params.push_back(model.lifetime_components[i].intensity);
        lt_names.push_back(std::format("intensity_{}", i+1));
    }

    text += print_params(lt_params, lt_names);

    text += "\n";
    text += std::format("{:=^100}\n", " Resolution Components ");
    text += "\n";

    double res_norm = std::accumulate(
        model.resolution_components.begin(),
        model.resolution_components.end(),
        0.,
        [](double acc, const auto& param) { return acc + param.intensity.val; }
    );

    if (std::abs(res_norm - 1) > 1e-3) {
        text += "!!! intensities don't add up to 100% !!!\n\n";
    }

    std::vector<core::fitting::FitParam> res_params;
    std::vector<std::string> res_names;

    for (const auto i : std::views::iota(0, n_r)) {
        res_params.push_back(model.resolution_components[i].fwhm);
        res_names.push_back(std::format("fwhm_{}", i+1));
    }
    for (const auto i : std::views::iota(0, n_r)) {
        res_params.push_back(model.resolution_components[i].intensity);
        res_names.push_back(std::format("intensity_{}", i+1));
    }
    for (const auto i : std::views::iota(0, n_r)) {
        res_params.push_back(model.resolution_components[i].t0);
        res_names.push_back(std::format("t0_{}", i+1));
    }

    text += print_params(res_params, res_names);

    text += "\n";
    text += thin_hline;

    return text;
}


}
