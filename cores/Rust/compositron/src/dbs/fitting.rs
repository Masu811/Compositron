use std::collections::HashMap;

use nalgebra::DVector;
use statrs::function::erf::erf;
use thiserror::Error;
use varpro::prelude::*;
use varpro::problem::*;
use varpro::solvers::levmar::LevMarSolver;

const TWO_OVER_SQRT_PI: f64 = 1.1283791670955126;

#[derive(Debug, Error)]
pub enum FitError {
    #[error("Data contains only NaN values")]
    AllNANValues,

    #[error("ModelBuilderError")]
    ModelBuilderError {
        #[from]
        inner: varpro::model::builder::error::ModelBuildError,
    },

    #[error("SeparableProblemBuilderError")]
    SeparableProblemBuilderError {
        #[from]
        inner: varpro::problem::SeparableProblemBuilderError,
    },

    #[error("StatisticsError")]
    StatisticsError {
        #[from]
        inner: varpro::statistics::Error<varpro::model::errors::ModelError>,
    },

    #[error("Something went wrong during least squares fit")]
    RuntimeError,
}

#[derive(Debug)]
pub struct Param {
    pub val: f64,
    pub err: f64,
}

fn gauss_basis(
    x: &DVector<f64>, x0: f64, sig: f64
) -> DVector<f64> {
    x.map(|x| {
        let u = (x - x0) / sig;
        (-0.5 * u*u).exp()
    })
}

fn gaussian_dx0(
    x: &DVector<f64>, x0: f64, sig: f64
) -> DVector<f64> {
    x.map(|x| {
        let u = (x - x0) / sig;
        (-0.5 * u*u).exp() * u / sig
    })
}

fn gaussian_dsig(
    x: &DVector<f64>, x0: f64, sig: f64
) -> DVector<f64> {
    x.map(|x| {
        let u = (x - x0) / sig;
        (-0.5 * u*u).exp() * u * u / sig
    })
}

pub fn fit_gaussian(
    x: &[f64], y: &[f64],
) -> Result<HashMap<&'static str, Param>, FitError> {
    let model = SeparableModelBuilder::<f64>::new(&["x0", "sig"])
        .function(&["x0", "sig"], gauss_basis)
        .partial_deriv("x0", gaussian_dx0)
        .partial_deriv("sig", gaussian_dsig)
        .independent_variable(DVector::from_column_slice(x))
        .initial_parameters(vec![511., 1.])
        .build()?;

    let problem = SeparableProblemBuilder::new(model)
        .observations(DVector::from_column_slice(y))
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

    let mut params = HashMap::<&'static str, Param>::new();

    params.insert("amp_1", Param { val: lin[0], err: errs[0] });
    params.insert("x0_1", Param { val: nonlin[0], err: errs[1] });
    params.insert("sig_1", Param { val: nonlin[1], err: errs[1] });

    Ok(params)
}

fn erf_basis(
    x: &DVector<f64>, x0: f64, sig: f64,
) -> DVector<f64> {
    x.map(|x| {
        let u = (x - x0) / sig;
        erf(u)
    })
}

fn erf_dx0(
    x: &DVector<f64>, x0: f64, sig: f64,
) -> DVector<f64> {
    x.map(|x| {
        let u = (x - x0) / sig;
        let e = (-u * u).exp();
        -TWO_OVER_SQRT_PI * e / sig
    })
}

fn erf_dsig(
    x: &DVector<f64>, x0: f64, sig: f64,
) -> DVector<f64> {
    x.map(|x| {
        let u = (x - x0) / sig;
        let e = (-u * u).exp();
        TWO_OVER_SQRT_PI * e * u / sig
    })
}

fn linear_basis(x: &DVector<f64>) -> DVector<f64> {
    x.clone()
}

fn const_basis(x: &DVector<f64>) -> DVector<f64> {
    DVector::from_element(x.len(), 1.)
}

/// Returns: Gauss amp, Gauss x0, Gauss sigma,
/// erf amp, linear coeff, const coeff
pub fn fit_erf_linear_1_gauss(
    x: &[f64], y: &[f64],
) -> Result<HashMap<&'static str, Param>, FitError> {
    let model = SeparableModelBuilder::<f64>::new(&["x0", "sig"])
        .function(&["x0", "sig"], gauss_basis)
        .partial_deriv("x0", gaussian_dx0)
        .partial_deriv("sig", gaussian_dsig)
        .function(&["x0", "sig"], erf_basis)
        .partial_deriv("x0", erf_dx0)
        .partial_deriv("sig", erf_dsig)
        .invariant_function(linear_basis)
        .invariant_function(const_basis)
        .independent_variable(DVector::from_column_slice(x))
        .initial_parameters(vec![511., 1.])
        .build()?;

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

    let mut params = HashMap::<&'static str, Param>::new();

    params.insert("amp_1", Param { val: lin[0], err: errs[0] });
    params.insert("x0_1", Param { val: nonlin[0], err: errs[4] });
    params.insert("sig_1", Param { val: nonlin[1], err: errs[5] });
    params.insert("erf_amp", Param { val: lin[1], err: errs[1] });
    params.insert("lin", Param { val: lin[2], err: errs[2] });
    params.insert("const", Param { val: lin[3], err: errs[3] });

    Ok(params)
}

pub fn fit_erf_linear_2_gauss(
    x: &[f64], y: &[f64],
) -> Result<HashMap<&'static str, Param>, FitError> {
    let model = SeparableModelBuilder::<f64>::new(
            &["x0_1", "sig_1", "x0_2", "sig_2"]
        )
        .function(&["x0_1", "sig_1"], gauss_basis)
        .partial_deriv("x0_1", gaussian_dx0)
        .partial_deriv("sig_1", gaussian_dsig)
        .function(&["x0_2", "sig_2"], gauss_basis)
        .partial_deriv("x0_2", gaussian_dx0)
        .partial_deriv("sig_2", gaussian_dsig)
        .function(&["x0_1", "sig_1"], erf_basis)
        .partial_deriv("x0_1", erf_dx0)
        .partial_deriv("sig_1", erf_dsig)
        .invariant_function(linear_basis)
        .invariant_function(const_basis)
        .independent_variable(DVector::from_column_slice(x))
        .initial_parameters(vec![511., 1., 510., 1.5])
        .build()?;

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

    let mut params = HashMap::<&'static str, Param>::new();

    params.insert("amp_1", Param { val: lin[0], err: errs[0] });
    params.insert("x0_1", Param { val: nonlin[0], err: errs[5] });
    params.insert("sig_1", Param { val: nonlin[1], err: errs[6] });
    params.insert("amp_2", Param { val: lin[1], err: errs[1] });
    params.insert("x0_2", Param { val: nonlin[2], err: errs[7] });
    params.insert("sig_2", Param { val: nonlin[3], err: errs[8] });
    params.insert("erf_amp", Param { val: lin[2], err: errs[2] });
    params.insert("lin", Param { val: lin[3], err: errs[3] });
    params.insert("const", Param { val: lin[4], err: errs[4] });

    Ok(params)
}

pub fn fit_erf_linear_3_gauss(
    x: &[f64], y: &[f64],
) -> Result<HashMap<&'static str, Param>, FitError> {
    let model = SeparableModelBuilder::<f64>::new(
            &["x0_1", "sig_1", "x0_2", "sig_2", "x0_3", "sig_3"]
        )
        .function(&["x0_1", "sig_1"], gauss_basis)
        .partial_deriv("x0_1", gaussian_dx0)
        .partial_deriv("sig_1", gaussian_dsig)
        .function(&["x0_2", "sig_2"], gauss_basis)
        .partial_deriv("x0_2", gaussian_dx0)
        .partial_deriv("sig_2", gaussian_dsig)
        .function(&["x0_3", "sig_3"], gauss_basis)
        .partial_deriv("x0_3", gaussian_dx0)
        .partial_deriv("sig_3", gaussian_dsig)
        .function(&["x0_1", "sig_1"], erf_basis)
        .partial_deriv("x0_1", erf_dx0)
        .partial_deriv("sig_1", erf_dsig)
        .invariant_function(linear_basis)
        .invariant_function(const_basis)
        .independent_variable(DVector::from_column_slice(x))
        .initial_parameters(vec![511., 1., 510., 1.5, 509., 2.0])
        .build()?;

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

    let mut params = HashMap::<&'static str, Param>::new();

    params.insert("amp_1", Param { val: lin[0], err: errs[0] });
    params.insert("x0_1", Param { val: nonlin[0], err: errs[6] });
    params.insert("sig_1", Param { val: nonlin[1], err: errs[7] });
    params.insert("amp_2", Param { val: lin[1], err: errs[1] });
    params.insert("x0_2", Param { val: nonlin[2], err: errs[8] });
    params.insert("sig_2", Param { val: nonlin[3], err: errs[9] });
    params.insert("amp_3", Param { val: lin[2], err: errs[2] });
    params.insert("x0_3", Param { val: nonlin[2], err: errs[10] });
    params.insert("sig_3", Param { val: nonlin[3], err: errs[11] });
    params.insert("erf_amp", Param { val: lin[3], err: errs[3] });
    params.insert("lin", Param { val: lin[4], err: errs[4] });
    params.insert("const", Param { val: lin[5], err: errs[5] });

    Ok(params)
}

pub fn background(
    x: &[f64], amp: f64, x0: f64, sig: f64, lin: f64, off: f64
) -> Vec<f64> {
    x.iter().map(|x| {
        let u = (x - x0) / sig;
        amp * erf(u) + lin * x + off
    }).collect()
}
