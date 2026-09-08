#pragma once

#include <ceres/cost_function.h>
#include <vector>
#include <algorithm>

#include <Eigen/Core>
#include <ceres/evaluation_callback.h>
#include <ceres/jet.h>

#include "compositron/constants.hpp"
#include "compositron/pals/model.hpp"


namespace compositron::pals::fitting {


Eigen::VectorXd lifetime_spectrum_component(
    const Eigen::VectorXd& t,
    double n,
    double t0,
    double tau,
    double l_int,
    double sig,
    double r_int,
    double r_t0
);


struct FitResult {
    model::LifetimeModel model;
    ceres::Solver::Summary summary;
    ceres::TerminationType fit_status;
    size_t n_dpoints;
    int n_eval;
    std::optional<double> red_chi_2;
    std::optional<Eigen::MatrixXd> cov;
};


class LifetimeEvalCallback : public ceres::EvaluationCallback {
    const Eigen::VectorXd& t;
    const Eigen::VectorXd& y;
    const Eigen::VectorXd& w;
    const std::vector<double>& p;
    size_t n_l;
    size_t n_r;
    double available_l_intensity;
    double available_r_intensity;
    std::vector<double> abs_l_intensities;
    std::vector<double> abs_r_intensities;
    std::vector<uint8_t> vary;
    int last_l_vary_idx;
    int last_r_vary_idx;


    void transform_intensities() {
        double leftover_intensity = available_l_intensity;

        for (size_t i = 0; i < n_l; ++i) {
            int idx = 3 + 2 * i + 1;

            if (vary[idx]) {
                double next_intensity = p[idx] * leftover_intensity;
                abs_l_intensities[i] = next_intensity;
                leftover_intensity = std::max(leftover_intensity - next_intensity, 0.);
            } else if (i == last_l_vary_idx) {
                abs_l_intensities[i] = leftover_intensity;
            } else {
                abs_l_intensities[i] = p[idx];
            }
        }

        leftover_intensity = available_r_intensity;

        for (size_t i = 0; i < n_r; ++i) {
            int idx = 3 + 2 * n_l + 3 * i + 1;

            if (vary[idx]) {
                double next_intensity = p[idx] * leftover_intensity;
                abs_r_intensities[i] = next_intensity;
                leftover_intensity = std::max(leftover_intensity - next_intensity, 0.);
            } else if (i == last_r_vary_idx) {
                abs_r_intensities[i] = leftover_intensity;
            } else {
                abs_r_intensities[i] = p[idx];
            }
        }
    }


public:
    void calc_res() {
        transform_intensities();

        double n = p[0];
        double bg = p[1];
        double t0 = p[2];

        residuals_.setConstant(bg);

        for (size_t i = 0; i < n_l; ++i) {
            size_t l_base_idx = 3 + 2 * i;
            double tau = p[l_base_idx];
            double l_int = abs_l_intensities[i];

            for (size_t j = 0; j < n_r; ++j) {
                size_t r_base_idx = 3 + 2 * n_l + 3 * j;
                double sig = p[r_base_idx];
                double r_int = abs_r_intensities[j];
                double r_t0 = p[r_base_idx + 2];

                residuals_ += lifetime_spectrum_component(
                    t, n, t0, tau, l_int, sig, r_int, r_t0
                );
            }
        }

        residuals_ -= y;
        residuals_ = residuals_.cwiseProduct(w);
    }


    void calc_res_and_jac() {
        transform_intensities();

        double n = p[0];
        double bg = p[1];
        double t0 = p[2];

        residuals_.setZero();

        size_t n_dpoints = t.size();

        Eigen::MatrixXd terms(n_dpoints, n_l * n_r);
        std::vector<Eigen::MatrixXd> deriv;

        for (size_t i = 0; i < n_l * n_r; ++i) {
            deriv.emplace_back(n_dpoints, 4);
        }

        for (size_t i = 0; i < n_l; ++i) {
            size_t l_base_idx = 3 + 2 * i;
            double tau = p[l_base_idx];
            double l_int = abs_l_intensities[i];

            for (size_t j = 0; j < n_r; ++j) {
                size_t col_idx = j + i * n_r;

                size_t r_base_idx = 3 + 2 * n_l + 3 * j;
                double sig = p[r_base_idx];
                double r_int = abs_r_intensities[j];
                double r_t0 = p[r_base_idx + 2];

                double sig_over_tau_sqrt_2 = sig / (tau * SQRT_2);
                double sig_over_tau_sqrt_2_sq = sig_over_tau_sqrt_2 * sig_over_tau_sqrt_2;
                double one_over_tau = 1 / tau;
                double one_over_tau_sq = one_over_tau * one_over_tau;
                double one_over_sig_sqrt_2 = 1 / (sig * SQRT_2);
                double one_over_tau_sqrt_2_pi = 1 / (tau * SQRT_2_PI);

                Eigen::VectorXd t_shifted = t.array() - (t0 + r_t0);

                Eigen::VectorXd a = -t_shifted.array() * one_over_tau + sig_over_tau_sqrt_2_sq;
                Eigen::VectorXd b = sig_over_tau_sqrt_2 - t_shifted.array() * one_over_sig_sqrt_2;

                Eigen::VectorXd decay = n * 0.5 * l_int * r_int * one_over_tau
                    * a.array().exp()
                    * b.array().unaryExpr([](auto x){ return std::erfc(x); });

                Eigen::VectorXd exp_comb = n * l_int * r_int * one_over_tau_sqrt_2_pi
                    * (a - b.cwiseProduct(b)).array().exp();

                terms.col(col_idx) = decay;

                residuals_ += decay;

                Eigen::VectorXd df_dt0 = -exp_comb / sig + decay / tau;

                deriv[col_idx].col(0) = df_dt0;
                deriv[col_idx].col(1) = exp_comb.array() * sig * one_over_tau_sq
                    + decay.array() * (
                        t_shifted.array() * one_over_tau_sq - sig * sig / (tau * tau * tau) - one_over_tau
                    );
                deriv[col_idx].col(2) = decay.array() * sig * one_over_tau_sq
                    - exp_comb.array() * (t_shifted.array() / (sig * sig) + one_over_tau);
                deriv[col_idx].col(3) = df_dt0;
            }
        }

        Eigen::MatrixXd df_dl_int(n_dpoints, n_l);

        for (size_t i = 0; i < n_l; ++i) {
            size_t idx = 3 + 2 * i + 1;

            if (vary[idx]) {
                Eigen::VectorXd tmp = Eigen::VectorXd::Zero(n_dpoints);

                for (size_t n = i; n < n_l; ++n) {
                    size_t n_idx = 3 + 2 * n + 1;
                    if (!vary[n_idx] && n != last_l_vary_idx) continue;

                    for (size_t m = 0; m < n_r; ++m) {
                        double a = (n == i) ?
                            a = 1 / (p[idx] + 1e-5) :
                            a = -1 / (1 - p[idx] + 1e-5);

                        tmp = tmp + a * terms.col(m + n * n_r);
                    }
                }

                df_dl_int.col(i) = tmp;
            } else {
                df_dl_int.col(i).array() = 0;
            }
        }

        Eigen::MatrixXd df_dr_int(n_dpoints, n_r);

        for (size_t i = 0; i < n_r; ++i) {
            size_t idx = 3 + 2 * n_l + 3 * i + 1;

            if (vary[idx]) {
                Eigen::VectorXd tmp = Eigen::VectorXd::Zero(n_dpoints);

                for (size_t m = i; m < n_r; ++m) {
                    size_t m_idx = 3 + 2 * n_l + 3 * n_r + 1;
                    if (!vary[m_idx] && m != last_r_vary_idx) continue;

                    for (size_t n = 0; n < n_l; ++n) {
                        double a = (m == i) ?
                            a = 1 / (p[idx] + 1e-5) :
                            a = -1 / (1 - p[idx] + 1e-5);

                        tmp = tmp + a * terms.col(m + n * n_r);
                    }
                }

                df_dr_int.col(i) = tmp;
            } else {
                df_dr_int.col(i).array() = 0;
            }
        }

        if (vary[0]) {
            jacobians_.row(0) = residuals_ / n;
        }

        if (vary[1]) {
            jacobians_.row(1).array() = 1;
        }

        if (vary[2]) {
            Eigen::VectorXd sum = Eigen::VectorXd::Zero(n_dpoints);
            for (size_t i = 0; i < deriv.size(); ++i) {
                sum = sum + deriv[i].col(0);
            }
            jacobians_.row(2) = sum;
        }

        for (size_t i = 0; i < n_l; ++i) {
            size_t flat_idx = 3 + 2 * i;

            if (vary[flat_idx]) {
                Eigen::VectorXd sum = Eigen::VectorXd::Zero(n_dpoints);
                for (size_t k = i * n_r; k < (i + 1) * n_r; ++k) {
                    sum = sum + deriv[k].col(1);
                }
                jacobians_.row(flat_idx) = sum;
            }

            if (vary[flat_idx + 1]) {
                jacobians_.row(flat_idx + 1) = df_dl_int.col(i);
            }
        }

        for (size_t i = 0; i < n_r; ++i) {
            size_t flat_idx = 3 + 2 * n_l + 3 * i;

            if (vary[flat_idx]) {
                Eigen::VectorXd sum = Eigen::VectorXd::Zero(n_dpoints);
                for (size_t k = i; k < n_l * n_r; k += n_r) {
                    sum = sum + deriv[k].col(2);
                }
                jacobians_.row(flat_idx) = sum;
            }

            if (vary[flat_idx + 1]) {
                jacobians_.row(flat_idx + 1) = df_dr_int.col(i);
            }

            if (vary[flat_idx + 2]) {
                Eigen::VectorXd sum = Eigen::VectorXd::Zero(n_dpoints);
                for (size_t k = i; k < n_l * n_r; k += n_r) {
                    sum = sum + deriv[k].col(3);
                }
                jacobians_.row(flat_idx + 2) = sum;
            }
        }

        jacobians_.array().rowwise() *= w.transpose().array();

        residuals_ = residuals_.array() + bg;
        residuals_ -= y;
        residuals_ = residuals_.cwiseProduct(w);
    }


private:
    Eigen::VectorXd residuals_;
    Eigen::MatrixXd jacobians_;
    bool jacobians_are_stale_ = true;


public:
    LifetimeEvalCallback(
        const Eigen::VectorXd& t,
        const Eigen::VectorXd& y,
        const Eigen::VectorXd& w,
        const std::vector<double>& p,
        const std::vector<uint8_t>& vary,
        size_t n_l,
        size_t n_r,
        int last_l_vary_idx,
        int last_r_vary_idx,
        double available_l_intensity,
        double available_r_intensity
    ) : t(t),
        y(y),
        w(w),
        p(p),
        n_l(n_l),
        n_r(n_r),
        vary(vary),
        last_l_vary_idx(last_l_vary_idx),
        last_r_vary_idx(last_r_vary_idx),
        available_l_intensity(available_l_intensity),
        available_r_intensity(available_r_intensity)
    {
        residuals_ = Eigen::VectorXd(t.size());
        jacobians_ = Eigen::MatrixXd(p.size(), t.size());

        abs_l_intensities = std::vector<double>(n_l, 0);
        abs_r_intensities = std::vector<double>(n_r, 0);

        PrepareForEvaluation(true, true);
    }


    void PrepareForEvaluation(
        bool evaluate_jacobians, bool new_evaluation_point
    ) final {
        if (new_evaluation_point) {
            if (evaluate_jacobians) {
                calc_res_and_jac();
                jacobians_are_stale_ = false;
            } else {
                calc_res();
                jacobians_are_stale_ = true;
            }
        } else if (evaluate_jacobians && jacobians_are_stale_) {
            calc_res_and_jac();
            jacobians_are_stale_ = false;
        }
    }


    const Eigen::VectorXd& residuals() const { return residuals_; }
    const Eigen::MatrixXd& jacobians() const { return jacobians_; }
    bool jacobians_are_stale() const { return jacobians_are_stale_; }
};


class LifetimeCostFunction : public ceres::CostFunction {
public:
    int index_ = -1;
    const LifetimeEvalCallback& callback_;
    int num_residuals_;

    LifetimeCostFunction(
        int index,
        const LifetimeEvalCallback& callback
    ) : index_(index), callback_(callback) {
        set_num_residuals(1);
        for (size_t i = 0; i < callback.jacobians().rows(); ++i) {
            mutable_parameter_block_sizes()->push_back(1);
        }
    }

    bool Evaluate(
        double const* const* params,
        double* residuals,
        double** jacobians
    ) const {
        residuals[0] = callback_.residuals()(index_);

        if (!jacobians) return true;

        for (size_t i = 0; i < parameter_block_sizes().size(); ++i) {
            if (jacobians[i] == nullptr) continue;
            jacobians[i][0] = callback_.jacobians()(i, index_);
        }

        return true;
    }
};


FitResult fit_lifetime_model(
    const Eigen::VectorXd& x,
    const Eigen::VectorXd& y,
    const model::LifetimeModel& model,
    uint64_t counts,
    double peak_center
);


} // namespace compositron::pals::fitting
