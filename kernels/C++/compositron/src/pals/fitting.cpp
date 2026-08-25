#include "compositron/pals/fitting.hpp"

#include <algorithm>
#include <optional>
#include <ranges>

#include <ceres/dynamic_autodiff_cost_function.h>
#include <ceres/problem.h>
#include <ceres/solver.h>
#include <ceres/types.h>

#include "compositron/constants.hpp"
#include "compositron/core/fitting.hpp"
#include "compositron/pals/model.hpp"


namespace compositron::pals::fitting {


namespace {


class LifetimeFit {
private:
    const Eigen::VectorXd& x;
    const Eigen::VectorXd& y;
    const model::LifetimeModel& input_model;
    std::optional<model::LifetimeModel> output_model;
    int n_l;
    int n_r;
    std::vector<uint8_t> l_vary;
    std::vector<uint8_t> r_vary;
    int last_l_vary_idx = -1;
    int last_r_vary_idx = -1;
    std::vector<double> opt;
    std::vector<double> values;
    std::vector<double> errors;
    ceres::Problem problem;
    ceres::Solver::Summary summary;


    std::vector<core::fitting::FitParam> flatten_params() {
        auto& model = output_model.value();

        std::vector<core::fitting::FitParam> params;

        params.push_back(model.scale_component);
        params.push_back(model.background_component);
        params.push_back(model.shift_component);

        for (auto& lt_comp : model.lifetime_components) {
            params.push_back(lt_comp.lifetime);
            params.push_back(lt_comp.intensity);
        }

        for (auto& res_comp : model.resolution_components) {
            params.push_back(res_comp.fwhm);
            params.push_back(res_comp.intensity);
            params.push_back(res_comp.t0);
        }

        return params;
    }


    void set_model_defaults(uint64_t counts, double peak_center) {
        auto& model = output_model.value();

        if (model.scale_component.default_initialized) {
            model.scale_component.val = 2 * counts;
            model.scale_component.min = 0;
        }

        if (model.shift_component.default_initialized) {
            double l = x[0];
            double r = x[x.size() - 1];
            double d = r - l;

            model.shift_component.min = l - d;
            model.shift_component.max = r + d;
            model.shift_component.val = peak_center;
        }

        for (auto& comp : model.lifetime_components) {
            comp.lifetime.min = std::max(0., comp.lifetime.min);
            comp.intensity.min = std::max(0., comp.intensity.min);
            comp.intensity.max = std::min(1., comp.intensity.max);
        }

        for (auto& comp : model.resolution_components) {
            comp.fwhm.min = std::max(0., comp.fwhm.min);
            comp.intensity.min = std::max(0., comp.intensity.min);
            comp.intensity.max = std::min(1., comp.intensity.max);
        }
    }


    void normalize_and_fix_params() {
        auto& model = output_model.value();

        double l_int_sum = std::accumulate(
            model.lifetime_components.begin(),
            model.lifetime_components.end(),
            0.,
            [](double acc, const auto& comp) { return acc + comp.intensity.val; }
        );

        for (auto& comp : model.lifetime_components) {
            comp.intensity.val /= l_int_sum;
        }

        double r_int_sum = std::accumulate(
            model.resolution_components.begin(),
            model.resolution_components.end(),
            0.,
            [](double acc, const auto& comp) { return acc + comp.intensity.val; }
        );

        for (auto& comp : model.resolution_components) {
            comp.intensity.val /= r_int_sum;
        }

        for (const int i : std::views::iota(1, n_l + 1)) {
            int l_int_idx = n_l - i;
            core::fitting::FitParam& param = model.lifetime_components[l_int_idx].intensity;
            if (param.vary) {
                last_l_vary_idx = l_int_idx;
                param.val = 1;
                param.vary = false;
                break;
            }
        }

        for (const int i : std::views::iota(1, n_r + 1)) {
            int r_int_idx = n_r - i;
            core::fitting::FitParam& param = model.resolution_components[r_int_idx].intensity;
            if (param.vary) {
                last_r_vary_idx = r_int_idx;
                param.val = 1;
                param.vary = false;
                break;
            }
        }

        if (std::ranges::all_of(
            model.resolution_components, [](const auto& comp) { return comp.t0.vary; }
        )) {
            model.resolution_components[0].t0.vary = false;
        }
    }


    void get_transformed_variables() {
        auto& model = output_model.value();

        opt.push_back(model.scale_component.val);
        opt.push_back(model.background_component.val);
        opt.push_back(model.shift_component.val);

        double remaining_l_int = 1;

        for (auto& lt_comp : model.lifetime_components) {
            opt.push_back(lt_comp.lifetime.val);
            opt.push_back(lt_comp.intensity.val / remaining_l_int);
            remaining_l_int -= lt_comp.intensity.val;
        }

        double remaining_r_int = 1;

        for (auto& res_comp : model.resolution_components) {
            opt.push_back(res_comp.fwhm.val / FWHM_OVER_SIGMA);
            opt.push_back(res_comp.intensity.val / remaining_r_int);
            remaining_r_int -= res_comp.intensity.val;
            opt.push_back(res_comp.t0.val);
        }
    }


    void build_problem() {
        auto& model = output_model.value();

        for (const auto& comp : model.lifetime_components) {
            l_vary.push_back(comp.intensity.vary);
        }

        for (const auto& comp : model.resolution_components) {
            r_vary.push_back(comp.intensity.vary);
        }

        auto params = flatten_params();
        size_t n_params = params.size();

        std::vector<double*> opt_ptrs;

        for (double& val : opt) {
            opt_ptrs.push_back(&val);
        }

        for (size_t i = 0; i < x.size(); ++i) {
            auto* cost_function =
                new ceres::DynamicAutoDiffCostFunction<LifetimeResidual, 1>(
                    x[i], y[i], 1 / std::sqrt(std::max(y[i], 1.)),
                    n_l, n_r, last_l_vary_idx, last_r_vary_idx, l_vary, r_vary
                );

            for (size_t i = 0; i < n_params; ++i) {
                cost_function->AddParameterBlock(1);
            }
            cost_function->SetNumResiduals(1);

            problem.AddResidualBlock(cost_function, nullptr, opt_ptrs);
        }

        for (size_t i = 0; i < n_params; ++i) {
            const auto& param = params[i];

            if (!param.vary) {
                problem.SetParameterBlockConstant(opt_ptrs[i]);
                continue;
            }
            if (std::isfinite(param.min)) {
                problem.SetParameterLowerBound(opt_ptrs[i], 0, param.min);
            }
            if (std::isfinite(param.max)) {
                problem.SetParameterUpperBound(opt_ptrs[i], 0, param.max);
            }
        }
    }


    std::vector<double> get_backtransformed_variables() {
        values.push_back(opt[0]);  // scale_component
        values.push_back(opt[1]);  // background_component
        values.push_back(opt[2]);  // shift_component

        double remaining_l_int = 1;

        for (auto i : std::views::iota(0, n_l)) {
            values.push_back(opt[3 + 2 * i]);  // lifetime
            if (l_vary[i]) {
                double l_int = remaining_l_int * opt[3 + 2 * i + 1];
                values.push_back(l_int);  // intensity
                remaining_l_int -= l_int;
            } else if (i == last_l_vary_idx) {
                double l_int = remaining_l_int;
                values.push_back(l_int);  // intensity
            } else {
                double l_int = opt[3 + 2 * i + 1];
                values.push_back(l_int);  // intensity
            }
        }

        double remaining_r_int = 1;

        for (auto i : std::views::iota(0, n_r)) {
            values.push_back(opt[3 + 2 * n_l + 3 * i] * FWHM_OVER_SIGMA);  // fwhm
            if (r_vary[i]) {
                double r_int = remaining_r_int * opt[3 + 2 * n_l + 3 * i + 1];
                values.push_back(r_int);  // intensity
                remaining_r_int -= r_int;
            } else if (i == last_r_vary_idx) {
                double r_int = remaining_r_int;
                values.push_back(r_int);  // intensity
            } else {
                double r_int = opt[3 + 2 * n_l + 3 * i + 1];
                values.push_back(r_int);  // intensity
            }
            values.push_back(opt[3 + 2 * n_l + 3 * i + 2]);  // t0
        }

        return values;
    }


    void paste_fit_result_values() {
        auto& model = output_model.value();

        model.scale_component.val = values[0];
        model.background_component.val = values[1];
        model.shift_component.val = values[2];

        for (auto i : std::views::iota(0, n_l)) {
            model.lifetime_components[i].lifetime.val = values[3 + 2 * i];
            model.lifetime_components[i].intensity.val = values[3 + 2 * i + 1];
        }

        for (auto i : std::views::iota(0, n_r)) {
            model.resolution_components[i].fwhm.val = values[3 + 2 * n_l + 3 * i];
            model.resolution_components[i].intensity.val = values[3 + 2 * n_l + 3 * i + 1];
            model.resolution_components[i].t0.val = values[3 + 2 * n_l + 3 * i + 2];
        }
    }


    void get_errors(Eigen::MatrixXd& cov) {
        errors.push_back(std::sqrt(cov(0, 0)));  // scale_component
        errors.push_back(std::sqrt(cov(1, 1)));  // background_component
        errors.push_back(std::sqrt(cov(2, 2)));  // shift_component

        double l_cumul_sq_err = 0;

        int cov_idx = 3;

        for (auto i : std::views::iota(0, n_l)) {
            errors.push_back(std::sqrt(cov(cov_idx, cov_idx)));  // lifetime
            cov_idx++;

            size_t i_idx = 3 + 2 * i + 1;

            if (l_vary[i]) {
                // intensity
                double err = std::sqrt(cov(cov_idx, cov_idx));
                cov_idx++;

                errors.push_back(
                    values[i_idx] * std::sqrt(
                        l_cumul_sq_err + std::pow(err / (opt[i_idx] + 1e-12), 2)
                    )
                );

                l_cumul_sq_err += std::pow(
                    err / (1 - opt[i_idx] + 1e-12), 2
                );
            } else if (i == last_l_vary_idx) {
                // intensity
                errors.push_back(values[i_idx] * std::sqrt(l_cumul_sq_err));
                cov_idx++;
            } else {
                // intensity
                errors.push_back(0);
                cov_idx++;
            }
        }

        double r_cumul_sq_err = 0;

        for (auto i : std::views::iota(0, n_r)) {
            errors.push_back(std::sqrt(cov(cov_idx, cov_idx)));  // fwhm
            cov_idx++;

            size_t i_idx = 3 + 2 * n_l + 3 * i + 1;

            if (r_vary[i]) {
                // intensity
                double err = std::sqrt(cov(cov_idx, cov_idx));
                cov_idx++;

                errors.push_back(
                    values[i_idx] * std::sqrt(
                        r_cumul_sq_err + std::pow(err / (opt[i_idx] + 1e-12), 2)
                    )
                );

                r_cumul_sq_err += std::pow(
                    err / (1 - opt[i_idx] + 1e-12), 2
                );
            } else if (i == last_r_vary_idx) {
                // intensity
                errors.push_back(values[i_idx] * std::sqrt(r_cumul_sq_err));
                cov_idx++;
            } else {
                // intensity
                errors.push_back(0);
                cov_idx++;
            }

            errors.push_back(std::sqrt(cov(cov_idx, cov_idx)));  // t0
            cov_idx++;
        }
    }


    void paste_fit_result_errors(model::LifetimeModel& model) {
        model.scale_component.err = errors[0];
        model.background_component.err = errors[1];
        model.shift_component.err = errors[2];

        for (auto i : std::views::iota(0, n_l)) {
            model.lifetime_components[i].lifetime.err = errors[3 + 2 * i];
            model.lifetime_components[i].intensity.err = errors[3 + 2 * i + 1];
        }

        for (auto i : std::views::iota(0, n_r)) {
            model.resolution_components[i].fwhm.err = errors[3 + 2 * n_l + 3 * i];
            model.resolution_components[i].intensity.err = errors[3 + 2 * n_l + 3 * i + 1];
            model.resolution_components[i].t0.err = errors[3 + 2 * n_l + 3 * i + 2];
        }
    }


    FitResult post_fit() {
        auto& model = output_model.value();

        int n_l = model.lifetime_components.size();
        int n_r = model.resolution_components.size();

        auto values = get_backtransformed_variables();

        paste_fit_result_values();

        FitResult fit_result {
            std::move(model),
            summary,
            summary.termination_type,
            summary.num_successful_steps,
            std::nullopt,
            std::nullopt
        };

        if (summary.termination_type != ceres::CONVERGENCE) {
            return fit_result;
        }

        fit_result.red_chi_2 = core::fitting::compute_red_chi2(summary);

        auto cov = core::fitting::compute_covariance_multi_blocks(problem, opt);

        if (!cov.has_value()) {
            return fit_result;
        }

        fit_result.cov = cov;

        get_errors(cov.value());

        paste_fit_result_errors(fit_result.model);

        if (last_l_vary_idx > -1) {
            fit_result.model.lifetime_components[last_l_vary_idx].intensity.vary = true;
        }
        if (last_r_vary_idx > -1) {
            fit_result.model.resolution_components[last_r_vary_idx].intensity.vary = true;
        }

        return fit_result;
    }


public:
    LifetimeFit() = delete;
    LifetimeFit(LifetimeFit& other) = delete;

    LifetimeFit(
        const Eigen::VectorXd& x,
        const Eigen::VectorXd& y,
        const model::LifetimeModel& input_model,
        uint64_t counts,
        double peak_center
    ) : x(x), y(y), input_model(input_model) {
        output_model = input_model;

        n_l = input_model.lifetime_components.size();
        n_r = input_model.resolution_components.size();

        set_model_defaults(counts, peak_center);

        normalize_and_fix_params();

        get_transformed_variables();

        build_problem();
    }


    FitResult fit() {
        ceres::Solver::Options options;
        options.linear_solver_type = ceres::DENSE_QR;

        Solve(options, &problem, &summary);

        return post_fit();
    }

};


} // anonymous namespace


FitResult fit_lifetime_model(
    const Eigen::VectorXd& x,
    const Eigen::VectorXd& y,
    const model::LifetimeModel& input_model,
    uint64_t counts,
    double peak_center
) {
    auto fitter = LifetimeFit(x, y, input_model, counts, peak_center);
    return fitter.fit();
}


} // namespace compositron::pals::fitting
