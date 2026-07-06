use std::f64::consts::SQRT_2;

use nalgebra::Const;
use nalgebra::DVector;
use nalgebra::Dyn;
use nalgebra::Matrix;
use nalgebra::ViewStorage;
use statrs::consts::SQRT_2PI;
use statrs::function::erf::erfc;
use varpro::model::builder::error::ModelBuildError;
use varpro::model::SeparableModel;
use varpro::prelude::*;
use varpro::problem::*;
use varpro::solvers::levmar::LevMarSolver;

use crate::core::fitting::{FitError, FitParam, FWHM_OVER_SIGMA};
use crate::pals::model::LifetimeModel;

pub fn free_lifetime_basis(
    t: &DVector<f64>,
    t0: f64,
    tau: f64,
    sig: f64,
    r_t0: f64,
    t0_param: FitParam,
    tau_param: FitParam,
    sig_param: FitParam,
    r_t0_param: FitParam,
) -> DVector<f64> {
    let t0 = t0_param.transform(t0);
    let tau = tau_param.transform(tau);
    let sig = sig_param.transform(sig);
    let r_t0 = r_t0_param.transform(r_t0);

    let sig_over_tau_sqrt_2 = sig / (tau * SQRT_2);
    let sig_over_tau_sqrt_2_sq = sig_over_tau_sqrt_2.powi(2);
    let one_over_tau = 1. / tau;
    let one_over_sig_sqrt_2 =  1. / (sig * SQRT_2);
    t.map(|t_i| {
        let tt_i = t_i - (t0 + r_t0);
        let e = (-tt_i * one_over_tau + sig_over_tau_sqrt_2_sq).exp();
        let c = erfc(sig_over_tau_sqrt_2 - tt_i * one_over_sig_sqrt_2);
        0.5 * one_over_tau * e * c
    })
}

fn free_lifetime_dt0(
    t: &DVector<f64>,
    t0: f64,
    tau: f64,
    sig: f64,
    r_t0: f64,
    t0_param: FitParam,
    tau_param: FitParam,
    sig_param: FitParam,
    r_t0_param: FitParam,
) -> DVector<f64> {
    let diff = t0_param.diff(t0);

    let t0 = t0_param.transform(t0);
    let tau = tau_param.transform(tau);
    let sig = sig_param.transform(sig);
    let r_t0 = r_t0_param.transform(r_t0);

    let sig_over_tau_sqrt_2 = sig / (tau * SQRT_2);
    let sig_over_tau_sqrt_2_sq = sig_over_tau_sqrt_2.powi(2);
    let one_over_tau = 1. / tau;
    let one_over_sig_sqrt_2 =  1. / (sig * SQRT_2);
    let one_over_tau_sqrt_two_pi = 1. / (tau * SQRT_2PI);
    t.map(|t_i| {
        let tt_i = t_i - (t0 + r_t0);
        let a = -tt_i * one_over_tau + sig_over_tau_sqrt_2_sq;
        let b = sig_over_tau_sqrt_2 - tt_i * one_over_sig_sqrt_2;
        let e = a.exp();
        let c = erfc(b);
        let decay = 0.5 * one_over_tau * e * c;
        let comb = one_over_tau_sqrt_two_pi * (a - b.powi(2)).exp();
        (-comb / sig + decay / tau) * diff
    })
}

fn free_lifetime_dtau(
    t: &DVector<f64>,
    t0: f64,
    tau: f64,
    sig: f64,
    r_t0: f64,
    t0_param: FitParam,
    tau_param: FitParam,
    sig_param: FitParam,
    r_t0_param: FitParam,
) -> DVector<f64> {
    let diff = tau_param.diff(tau);

    let t0 = t0_param.transform(t0);
    let tau = tau_param.transform(tau);
    let sig = sig_param.transform(sig);
    let r_t0 = r_t0_param.transform(r_t0);

    let sig_over_tau_sqrt_2 = sig / (tau * SQRT_2);
    let sig_over_tau_sqrt_2_sq = sig_over_tau_sqrt_2.powi(2);
    let one_over_tau = 1. / tau;
    let one_over_sig_sqrt_2 =  1. / (sig * SQRT_2);
    let one_over_tau_sqrt_two_pi = 1. / (SQRT_2PI * tau);
    t.map(|t_i| {
        let tt_i = t_i - (t0 + r_t0);
        let a = -tt_i * one_over_tau + sig_over_tau_sqrt_2_sq;
        let b = sig_over_tau_sqrt_2 - tt_i * one_over_sig_sqrt_2;
        let e = a.exp();
        let c = erfc(b);
        let decay = 0.5 * one_over_tau * e * c;
        let comb = one_over_tau_sqrt_two_pi * (a - b.powi(2)).exp();
        (
            comb * sig * one_over_tau.powi(2)
            + decay * (
                -sig.powi(2) * one_over_tau.powi(3)
                + tt_i * one_over_tau.powi(2)
                - one_over_tau
            )
        ) * diff
    })
}

fn free_lifetime_dsig(
    t: &DVector<f64>,
    t0: f64,
    tau: f64,
    sig: f64,
    r_t0: f64,
    t0_param: FitParam,
    tau_param: FitParam,
    sig_param: FitParam,
    r_t0_param: FitParam,
) -> DVector<f64> {
    let diff = sig_param.diff(sig);

    let t0 = t0_param.transform(t0);
    let tau = tau_param.transform(tau);
    let sig = sig_param.transform(sig);
    let r_t0 = r_t0_param.transform(r_t0);

    let sig_over_tau_sqrt_2 = sig / (tau * SQRT_2);
    let sig_over_tau_sqrt_2_sq = sig_over_tau_sqrt_2.powi(2);
    let one_over_tau = 1. / tau;
    let one_over_sig_sqrt_2 =  1. / (sig * SQRT_2);
    let one_over_tau_sqrt_two_pi = 1. / (SQRT_2PI * tau);
    t.map(|t_i| {
        let tt_i = t_i - (t0 + r_t0);
        let a = -tt_i * one_over_tau + sig_over_tau_sqrt_2_sq;
        let b = sig_over_tau_sqrt_2 - tt_i * one_over_sig_sqrt_2;
        let e = a.exp();
        let c = erfc(b);
        let decay = 0.5 * one_over_tau * e * c;
        let comb = one_over_tau_sqrt_two_pi * (a - b.powi(2)).exp();
        (
            decay * sig * one_over_tau.powi(2)
            - comb * (one_over_tau + tt_i / sig.powi(2))
        ) * diff
    })
}

fn build_free_model(
    x: &[f64], lt_model: &LifetimeModel
) -> Result<SeparableModel<f64>, ModelBuildError> {
    let n_l = lt_model.lifetime_components.len();
    let n_r = lt_model.resolution_components.len();

    let mut variables = Vec::new();
    let mut initial_values = Vec::new();

    variables.push("t0".into());
    initial_values.push(2572.);

    for (i, component) in lt_model.lifetime_components.iter().enumerate() {
        if component.lifetime.vary {
            variables.push(format!("lifetime_{}", i + 1));
            initial_values.push(component.lifetime.backtransform(component.lifetime.val));
        }
    }

    for (i, component) in lt_model.resolution_components.iter().enumerate() {
        if component.fwhm.vary {
            variables.push(format!("res_sigma_{}", i + 1));
            // TODO: convert fwhm into sigma
            initial_values.push(component.fwhm.backtransform(component.fwhm.val));
        }
        if component.t0.vary {
            variables.push(format!("res_t0_{}", i + 1));
            initial_values.push(component.t0.backtransform(component.t0.val));
        }
    }

    let mut builder = SeparableModelBuilder::<f64>::new(&variables);

    for i in 0..n_l {
        for j in 0..n_r {
            let tau = format!("lifetime_{}", i + 1);
            let sig = format!("res_sigma_{}", j + 1);
            let r_t0 = format!("res_t0_{}", j + 1);
            let t0_param = lt_model.shift_component;
            let tau_param = lt_model.lifetime_components[i].lifetime;
            let sig_param = lt_model.resolution_components[j].fwhm;
            let r_t0_param = lt_model.resolution_components[j].t0;
            builder = builder.function(
                    ["t0", &tau, &sig, &r_t0],
                    move |t: &DVector<f64>, a: f64, b: f64, c: f64, d: f64| {
                        free_lifetime_basis(
                            t, a, b, c, d,
                            t0_param, tau_param, sig_param, r_t0_param
                        )
                    }
                )
                .partial_deriv("t0",
                    move |t: &DVector<f64>, a: f64, b: f64, c: f64, d: f64| {
                        free_lifetime_dt0(
                            t, a, b, c, d,
                            t0_param, tau_param, sig_param, r_t0_param
                        )
                    }
                )
                .partial_deriv(tau,
                    move |t: &DVector<f64>, a: f64, b: f64, c: f64, d: f64| {
                        free_lifetime_dtau(
                            t, a, b, c, d,
                            t0_param, tau_param, sig_param, r_t0_param
                        )
                    }
                )
                .partial_deriv(sig,
                    move |t: &DVector<f64>, a: f64, b: f64, c: f64, d: f64| {
                        free_lifetime_dsig(
                            t, a, b, c, d,
                            t0_param, tau_param, sig_param, r_t0_param
                        )
                    }
                )
                .partial_deriv(r_t0,
                    move |t: &DVector<f64>, a: f64, b: f64, c: f64, d: f64| {
                        free_lifetime_dt0(
                            t, d, b, c, a,
                            r_t0_param, tau_param, sig_param, t0_param
                        )
                    }
                );
        }
    }

    builder
        .independent_variable(DVector::from_column_slice(x))
        .initial_parameters(initial_values)
        .build()
}

fn calculate_intensities(
    lin: &Matrix<f64, Dyn, Const<1>, ViewStorage<'_, f64, Dyn, Const<1>, Const<1>, Dyn>>,
    n_l: usize,
    n_r: usize,
) -> Vec<f64> {
    let mut intensities = Vec::new();

    let n = lin.sum();

    for i in 0..n_l {
        let i_i = lin.rows(i * n_r, n_r).sum() / n;
        intensities.push(i_i);
    }

    for i in 0..n_r {
        let j_i = lin[i] / intensities[0] / n;
        intensities.push(j_i);
    }

    intensities.push(n);

    intensities
}

pub fn fit_free_intensity_lifetime_spectrum(
    x: &[f64], y: &[f64], lt_model: &mut LifetimeModel
) -> Result<(), FitError> {
    let model = build_free_model(x, lt_model)?;

    let weights = y
        .iter()
        .map(|&i| 1. / i.max(1.).sqrt())
        .collect::<Vec<f64>>();

    let problem = SeparableProblemBuilder::new(model)
        .observations(DVector::from_column_slice(y))
        .weights(DVector::from_vec(weights))
        .build()?;

    let Ok(fit_result) = LevMarSolver::default().solve(problem) else {
        return Err(FitError::RuntimeError);
    };

    let nonlin = fit_result.nonlinear_parameters();

    let Some(lin) = fit_result.linear_coefficients() else {
        return Err(FitError::RuntimeError);
    };

    let stats = varpro::statistics::FitStatistics::try_from(&fit_result)?;

    let cov = stats.covariance_matrix();

    let errs = cov.diagonal().map(|x| x.sqrt());

    let n_l = lt_model.lifetime_components.len();
    let n_r = lt_model.resolution_components.len();

    let intensities = calculate_intensities(&lin, n_l, n_r);

    for i in 0..n_l {
        println!(
            "Lifetime {}: {}",
            i + 1,
            lt_model.lifetime_components[i].lifetime.transform(nonlin[i])
        );
        println!(
            "Intensity {}: {}",
            i + 1,
            intensities[i],
        );
    }

    for i in 0..lt_model.resolution_components.len() {
        println!(
            "Sigma {}: {}",
            i + 1,
            lt_model.resolution_components[i].fwhm.transform(nonlin[2*i+n_l])
        );
        println!(
            "Intensity {}: {}",
            i + 1,
            intensities[n_l + i],
        );
        println!(
            "t0 {}: {}",
            i + 1,
            lt_model.resolution_components[i].t0.transform(nonlin[2*i+1+n_l])
        );
    }

    println!("N: {}", intensities[n_l + n_r]);

    Ok(())
}
