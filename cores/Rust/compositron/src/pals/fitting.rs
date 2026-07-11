use std::collections::HashMap;
use std::f64::{self, consts::SQRT_2};

use nalgebra::{DMatrix, DVector, Dyn, Owned};
use statrs::consts::SQRT_2PI;
use statrs::function::erf::erfc;
use levenberg_marquardt::{LeastSquaresProblem, LevenbergMarquardt, TerminationReason};
use thiserror::Error;

use crate::pals::model::{LifetimeModel};
use crate::core::fitting::{FitParam, FWHM_OVER_SIGMA};

enum ParamType {
    Varied(usize),
    Fixed(usize),
    LastVariedIntensity(usize),
}

struct Problem<'a> {
    x: DVector<f64>,
    y: DVector<f64>,
    w: DVector<f64>,
    p: DVector<f64>,
    n_l: usize,
    n_r: usize,
    transformed_params: DVector<f64>,
    fixed_params: DVector<f64>,
    param_types: Vec<ParamType>,
    varied_param_names: Vec<String>,
    params: Vec<&'a FitParam>,
    available_l_intensity: f64,
    available_r_intensity: f64,
    absolute_l_intensities: Vec<f64>,
    absolute_r_intensities: Vec<f64>,
    transforms: Vec<fn (f64, &FitParam) -> f64>,
    backtransforms: Vec<fn (f64, &FitParam) -> f64>,
    diffs: Vec<fn (f64, &FitParam) -> f64>,
}

impl<'a> Problem<'a> {
    fn set_init(&mut self) {
        for i in 0..self.p.len() {
            self.p[i] = self.backtransforms[i](
                self.transformed_params[i], self.params[i]
            );
        }
    }

    fn get_param(&self, i: usize) -> f64 {
        match self.param_types[i] {
            ParamType::Varied(idx) => self.transformed_params[idx],
            ParamType::Fixed(idx) => self.fixed_params[idx],
            ParamType::LastVariedIntensity(idx) => self.fixed_params[idx],
        }
    }

    fn make_intensities(&mut self) {
        let mut leftover_intensity = self.available_l_intensity;

        for i in 0..self.n_l {
            match self.param_types[3 + 2 * i + 1] {
                ParamType::Varied(idx) => {
                    let next_intensity = self.transformed_params[idx] * leftover_intensity;
                    self.absolute_l_intensities[i] = next_intensity;
                    leftover_intensity = (leftover_intensity - next_intensity).max(0.);
                },
                ParamType::Fixed(idx) => {
                    self.absolute_l_intensities[i] = self.fixed_params[idx];
                },
                ParamType::LastVariedIntensity(_) => {
                    self.absolute_l_intensities[i] = leftover_intensity;
                },
            }
        }

        let mut leftover_intensity = self.available_r_intensity;

        for i in 0..self.n_r {
            match self.param_types[3 + 2 * self.n_l + 3 * i + 1] {
                ParamType::Varied(idx) => {
                    let next_intensity = self.transformed_params[idx] * leftover_intensity;
                    self.absolute_r_intensities[i] = next_intensity;
                    leftover_intensity = (leftover_intensity - next_intensity).max(0.);
                },
                ParamType::Fixed(idx) => {
                    self.absolute_r_intensities[i] = self.fixed_params[idx];
                },
                ParamType::LastVariedIntensity(_) => {
                    self.absolute_r_intensities[i] = leftover_intensity;
                },
            }
        }
    }
}

impl<'a> LeastSquaresProblem<f64, Dyn, Dyn> for Problem<'a> {
    type ParameterStorage = Owned<f64, Dyn>;
    type ResidualStorage = Owned<f64, Dyn>;
    type JacobianStorage = Owned<f64, Dyn, Dyn>;

    fn set_params(&mut self, params: &DVector<f64>) {
        self.p.copy_from(params);
        for i in 0..params.len() {
            self.transformed_params[i] = self.transforms[i](
                params[i], self.params[i]
            );
        }
        self.make_intensities();
    }

    fn params(&self) -> DVector<f64> {
        self.p.clone()
    }

    fn residuals(&self) -> Option<DVector<f64>> {
        let mut f = DVector::from_element(self.x.len(), 0.);

        let n_l = self.n_l;
        let n_r = self.n_r;

        let n = self.get_param(0);
        let bg = self.get_param(1);
        let t0 = self.get_param(2);

        let mut p = vec![0.; 2 * n_l + 3 * n_r];

        for i in 0..n_l {
            // tau
            p[2 * i] = self.get_param(3 + 2 * i);
            // l_int
            p[2 * i + 1] = self.absolute_l_intensities[i];
        }

        for i in 0..n_r {
            // sig
            p[2 * n_l + 3 * i] = self.get_param(3 + 2 * n_l + 3 * i);
            // r_int
            p[2 * n_l + 3 * i + 1] = self.absolute_r_intensities[i];
            // r_t0
            p[2 * n_l + 3 * i + 2] = self.get_param(3 + 2 * n_l + 3 * i + 2);
        }

        for i in 0..n_l {
            let tau = p[2 * i];
            let l_int = p[2 * i + 1];

            for j in 0..n_r {
                let sig = p[2 * n_l + 3 * j];
                let r_int = p[2 * n_l + 3 * j + 1];
                let r_t0 = p[2 * n_l + 3 * j + 2];

                let sig_over_tau_sqrt_2 = sig / (tau * SQRT_2);
                let sig_over_tau_sqrt_2_sq = sig_over_tau_sqrt_2.powi(2);
                let one_over_tau = 1. / tau;
                let one_over_sig_sqrt_2 =  1. / (sig * SQRT_2);
                f += self.x.map(|t_i| {
                    let tt_i = t_i - (t0 + r_t0);
                    let e = (-tt_i * one_over_tau + sig_over_tau_sqrt_2_sq).exp();
                    let c = erfc(sig_over_tau_sqrt_2 - tt_i * one_over_sig_sqrt_2);
                    0.5 * l_int * r_int * one_over_tau * e * c
                });
            }
        }

        f *= n;
        if bg > 0. {
            f.add_scalar_mut(bg);
        }

        let r = (&self.y - f).component_mul(&self.w);

        if r.iter().any(|x| x.is_nan()) {
            // TODO: log
        }

        Some(r)
    }

    fn jacobian(&self) -> Option<DMatrix<f64>> {
        let mut f = DVector::from_element(self.x.len(), 0.);

        let n_l = self.n_l;
        let n_r = self.n_r;

        let n_dpoints = self.x.len();
        let mut terms = DMatrix::zeros(n_dpoints, n_l * n_r);
        let mut deriv = Vec::with_capacity(n_l * n_r);

        let n = self.get_param(0);
        let t0 = self.get_param(2);

        let mut p = vec![0.; 2 * n_l + 3 * n_r];

        for i in 0..n_l {
            // tau
            p[2 * i] = self.get_param(3 + 2 * i);
            // l_int
            p[2 * i + 1] = self.absolute_l_intensities[i];
        }

        for i in 0..n_r {
            // sig
            p[2 * n_l + 3 * i] = self.get_param(3 + 2 * n_l + 3 * i);
            // r_int
            p[2 * n_l + 3 * i + 1] = self.absolute_r_intensities[i];
            // r_t0
            p[2 * n_l + 3 * i + 2] = self.get_param(3 + 2 * n_l + 3 * i + 2);
        }

        for i in 0..n_l {
            let tau = p[2 * i];
            let l_int = p[2 * i + 1];

            for j in 0..n_r {
                let sig = p[2 * n_l + 3 * j];
                let r_int = p[2 * n_l + 3 * j + 1];
                let r_t0 = p[2 * n_l + 3 * j + 2];

                let sig_over_tau_sqrt_2 = sig / (tau * SQRT_2);
                let sig_over_tau_sqrt_2_sq = sig_over_tau_sqrt_2.powi(2);
                let one_over_tau = 1. / tau;
                let one_over_tau_sqrt_2_pi = 1. / (tau * SQRT_2PI);
                let one_over_sig_sqrt_2 =  1. / (sig * SQRT_2);

                let tt = self.x.map(|t_i| t_i - (t0 + r_t0));

                let a = tt.map(|tt_i| {
                    -tt_i * one_over_tau + sig_over_tau_sqrt_2_sq
                });

                let b = tt.map(|tt_i| {
                    sig_over_tau_sqrt_2 - tt_i * one_over_sig_sqrt_2
                });

                let decay = DVector::from_fn(a.len(), |i, _| {
                    let e = a[i].exp();
                    let c = erfc(b[i]);
                    n * 0.5 * l_int * r_int * one_over_tau * e * c
                });

                let exp_comb = DVector::from_fn(a.len(), |i, _| {
                    n * l_int * r_int * one_over_tau_sqrt_2_pi
                        * (a[i] - b[i].powi(2)).exp()
                });

                terms.set_column(j + i * n_r, &decay);

                f += &decay;

                let df_dt0 = -&exp_comb / sig + &decay / tau;
                let df_dtau = &exp_comb * sig / tau.powi(2)
                    + &decay.component_mul(&(&tt / tau.powi(2)).add_scalar(
                        -sig.powi(2) / tau.powi(3) - one_over_tau
                    ));
                let df_dsig = &decay * sig / tau.powi(2)
                    - &exp_comb.component_mul(
                        &(&tt / sig.powi(2)).add_scalar(one_over_tau)
                    );

                let mut d = DMatrix::zeros(n_dpoints, 4);

                d.set_column(0, &df_dt0);
                d.set_column(1, &df_dtau);
                d.set_column(2, &df_dsig);
                d.set_column(3, &df_dt0);

                deriv.push(d);
            }
        }

        let mut df_dl_int = Vec::new();

        for i in 0..n_l {
            match self.param_types[3 + 2 * i + 1] {
                ParamType::Varied(idx) => {
                    let mut tmp = DVector::zeros(n_dpoints);

                    for n in i..n_l {
                        for m in 0..n_r {
                            let a = if n == i {
                                1. / (self.transformed_params[idx] + 1e-5)
                            } else {
                                -1. / (1. - self.transformed_params[idx] + 1e-5)
                            };

                            tmp.axpy(a, &terms.column(n * n_r + m), 1.);
                        }
                    }

                    df_dl_int.push(tmp);
                },
                _ => df_dl_int.push(DVector::zeros(n_dpoints)),
            }
        }

        let mut df_dr_int = Vec::new();

        for i in 0..n_r {
            match self.param_types[3 + 2 * n_l + 3 * i + 1] {
                ParamType::Varied(idx) => {
                    let mut tmp = DVector::zeros(n_dpoints);

                    for m in i..n_r {
                        for n in 0..n_l {
                            let a = if m == i {
                                1. / (self.transformed_params[idx] + 1e-5)
                            } else {
                                -1. / (1. - self.transformed_params[idx] + 1e-5)
                            };

                            tmp.axpy(a, &terms.column(n * n_r + m), 1.);
                        }
                    }

                    df_dr_int.push(tmp);
                },
                _ => df_dr_int.push(DVector::zeros(n_dpoints)),
            }
        }

        let mut derivatives = HashMap::<String, DVector<f64>>::new();

        derivatives.insert("N".into(), &f / n);
        derivatives.insert(
            "background".into(), DVector::from_element(n_dpoints, 1.)
        );
        derivatives.insert("t0".into(), deriv.iter().fold(
            DVector::zeros(n_dpoints),
            |acc, mat| acc + mat.column(0),
        ));

        for i in 0..n_l {
            let mut tmp = DVector::zeros(n_dpoints);
            for j in i * n_r..(i + 1) * n_r {
                tmp += &deriv[j].column(1);
            }
            derivatives.insert(format!("lifetime_{}", i + 1), tmp);
            derivatives.insert(
                format!("intensity_{}", i + 1), df_dl_int[i].clone()
            );
        }

        for i in 0..n_r {
            let mut tmp = DVector::zeros(n_dpoints);
            for j in (i..n_l * n_r).step_by(n_r) {
                tmp += &deriv[j].column(2);
            }
            derivatives.insert(format!("res_fwhm_{}", i + 1), tmp);

            let mut tmp = DVector::zeros(n_dpoints);
            for j in (i..n_l * n_r).step_by(n_r) {
                tmp += &deriv[j].column(3);
            }
            derivatives.insert(format!("res_t0_{}", i + 1), tmp);
            derivatives.insert(
                format!("res_intensity_{}", i + 1), df_dr_int[i].clone()
            );
        }

        let mut j = DMatrix::zeros(n_dpoints, self.p.len());

        for (i, param) in self.varied_param_names.iter().enumerate() {
            j.set_column(i, &(
                -derivatives.get(param.into()).unwrap().component_mul(
                    &(&self.w * self.diffs[i](self.p[i], self.params[i]))
                )
            ));

            if j.iter().any(|x| x.is_nan()) {
                // TODO: log
            }
        }

        Some(j)
    }
}

fn compute_red_chi2(
    problem: &impl levenberg_marquardt::LeastSquaresProblem<f64, Dyn, Dyn>,
) -> Option<f64> {
    let r = problem.residuals()?;
    let j = problem.jacobian()?;

    let m = r.len();
    let n = j.ncols();

    if m <= n {
        return None;
    }

    let dof = (m - n) as f64;
    let red_chi2 = r.norm_squared() / dof;

    Some(red_chi2)
}

fn compute_covariance(
    problem: &Problem, red_chi2: f64
) -> Option<DMatrix<f64>> {
    let mut j = problem.jacobian()?;

    let d = DVector::from_iterator(
        problem.p.len(),
        problem.p
            .iter()
            .enumerate()
            .map(|(i, &p)| problem.diffs[i](p, problem.params[i]))
    );

    for i in 0..problem.p.len() {
        j.column_mut(i).scale_mut(1. / (d[i] + 1e-6));
    }

    let jtj = j.transpose() * &j;

    let jtj_inv = jtj.try_inverse()?;

    Some(jtj_inv * red_chi2)
}

struct ProblemTemplate<'a> {
    transformed_params: Vec<f64>,
    fixed_params: Vec<f64>,
    param_types: Vec<ParamType>,
    varied_param_names: Vec<String>,
    params: Vec<&'a FitParam>,
    available_l_intensity: f64,
    available_r_intensity: f64,
    transforms: Vec<fn (f64, &FitParam) -> f64>,
    backtransforms: Vec<fn (f64, &FitParam) -> f64>,
    diffs: Vec<fn (f64, &FitParam) -> f64>,
}

enum ParamKind {
    LIntensity,
    RIntensity,
    Sigma,
    Other,
}

impl<'a> ProblemTemplate<'a> {
    fn gobble(
        &mut self,
        name: String,
        param: &'a FitParam,
        param_kind: ParamKind,
    ) {
        if param.vary {
            self.transformed_params.push(match param_kind {
                ParamKind::Sigma => param.val / FWHM_OVER_SIGMA,
                _ => param.val,
            });
            self.param_types.push(
                ParamType::Varied(self.transformed_params.len() - 1)
            );
            self.varied_param_names.push(name);
            self.params.push(param);
            self.transforms.push(param.transform);
            self.backtransforms.push(param.backtransform);
            self.diffs.push(param.diff);
        } else {
            self.fixed_params.push(match param_kind {
                ParamKind::Sigma => param.val / FWHM_OVER_SIGMA,
                _ => param.val,
            });
            self.param_types.push(
                ParamType::Fixed(self.fixed_params.len() - 1)
            );
            match param_kind {
                ParamKind::LIntensity => self.available_l_intensity -= param.val,
                ParamKind::RIntensity => self.available_r_intensity -= param.val,
                _ => {},
            }
        }
    }
}

fn prepare_fit(
    model: &mut LifetimeModel,
    x: &[f64],
    counts: f64,
    peak_center: f64,
) -> (Option<usize>, Option<usize>) {
    if model.scale_component.default_initialized {
        model.scale_component.val = 2. * counts;
        model.scale_component = model.scale_component.min(0.);
    }

    if model.shift_component.default_initialized {
        let l = x[0];
        let r = x[x.len() - 1];
        let d = r - l;
        model.shift_component = model.shift_component.min(l - d);
        model.shift_component = model.shift_component.max(r + d);
        model.shift_component.val = peak_center;
    }

    if model.background_component.default_initialized {
        model.background_component.vary = false;
    }

    let n_l = model.lifetime_components.len();
    let n_r = model.resolution_components.len();

    let l_int_sum = model.lifetime_components
        .iter()
        .fold(0., |cum, comp| cum + comp.intensity.val);
    let r_int_sum = model.resolution_components
        .iter()
        .fold(0., |cum, comp| cum + comp.intensity.val);

    model.lifetime_components.iter_mut().for_each(|comp| comp.intensity.val /= l_int_sum);
    model.resolution_components.iter_mut().for_each(|comp| comp.intensity.val /= r_int_sum);

    let mut last_l_vary_idx = None;

    for i in 1..=n_l {
        let l_int_idx = n_l - i;
        let param = &mut model.lifetime_components[l_int_idx].intensity;
        if param.vary {
            last_l_vary_idx = Some(l_int_idx);
            param.val = 0.;
            param.vary = false;
            break;
        }
    }

    let mut last_r_vary_idx = None;

    for i in 1..=n_r {
        let r_int_idx = n_r - i;
        let param = &mut model.resolution_components[r_int_idx].intensity;
        if param.vary {
            last_r_vary_idx = Some(r_int_idx);
            param.val = 0.;
            param.vary = false;
            break;
        }
    }

    if model.resolution_components.iter().all(|comp| comp.t0.vary) {
        model.resolution_components[0].t0.val = 0.;
        model.resolution_components[0].t0.vary = false;
    }

    (last_l_vary_idx, last_r_vary_idx)
}

fn build_problem<'a>(
    x: &[f64], y: &[f64], model: &'a LifetimeModel
) -> Problem<'a> {
    let w = DVector::from_iterator(
        y.len(), y.iter().map(|y_i| 1. / y_i.max(1.).sqrt())
    );

    let n_l = model.lifetime_components.len();
    let n_r = model.resolution_components.len();

    let mut tmp = ProblemTemplate {
        transformed_params: Vec::new(),
        fixed_params: Vec::new(),
        param_types: Vec::new(),
        varied_param_names: Vec::new(),
        params: Vec::new(),
        available_l_intensity: 1.,
        available_r_intensity: 1.,
        transforms: Vec::new(),
        backtransforms: Vec::new(),
        diffs: Vec::new(),
    };

    tmp.gobble("N".into(), &model.scale_component, ParamKind::Other);
    tmp.gobble("background".into(), &model.background_component, ParamKind::Other);
    tmp.gobble("t0".into(), &model.shift_component, ParamKind::Other);

    for i in 0..n_l {
        let c = &model.lifetime_components[i];
        tmp.gobble(format!("lifetime_{}", i + 1), &c.lifetime, ParamKind::Other);
        tmp.gobble(format!("intensity_{}", i + 1), &c.intensity, ParamKind::LIntensity);
    }

    for i in 0..n_r {
        let c = &model.resolution_components[i];
        tmp.gobble(format!("res_fwhm_{}", i + 1), &c.fwhm, ParamKind::Sigma);
        tmp.gobble(format!("res_intensity_{}", i + 1), &c.intensity, ParamKind::RIntensity);
        tmp.gobble(format!("res_t0_{}", i + 1), &c.t0, ParamKind::Other);
    }

    let problem = Problem {
        x: DVector::from_column_slice(x),
        y: DVector::from_column_slice(y),
        p: DVector::from_element(tmp.transformed_params.len(), f64::NAN),
        w,
        n_l: model.lifetime_components.len(),
        n_r: model.resolution_components.len(),
        transformed_params: DVector::from_vec(tmp.transformed_params),
        fixed_params: DVector::from_vec(tmp.fixed_params),
        param_types: tmp.param_types,
        varied_param_names: tmp.varied_param_names,
        params: tmp.params,
        available_l_intensity: tmp.available_l_intensity,
        available_r_intensity: tmp.available_r_intensity,
        absolute_l_intensities: vec![0.; n_l],
        absolute_r_intensities: vec![0.; n_r],
        transforms: tmp.transforms,
        backtransforms: tmp.backtransforms,
        diffs: tmp.diffs,
    };

    problem
}

fn paste_val_and_err(
    param: &mut FitParam,
    param_idx: usize,
    param_types: &Vec<ParamType>,
    transformed_params: &DVector<f64>,
    errors: Option<&DVector<f64>>
) {
    match param_types[param_idx] {
        ParamType::Varied(idx) => {
            param.val = transformed_params[idx];
            if let Some(errs) = &errors {
                param.err = errs[idx];
            }
        },
        ParamType::Fixed(_) => {
            param.err = 0.;
        },
        _ => {},
    }
}

fn transform_intensity_errors(
    model: &mut LifetimeModel,
    transformed_params: &DVector<f64>,
    param_types: &Vec<ParamType>,
    absolute_l_intensities: &Vec<f64>,
    absolute_r_intensities: &Vec<f64>,
    last_l_vary_idx: Option<usize>,
    last_r_vary_idx: Option<usize>,
    errors: Option<&DVector<f64>>,
) {
    let Some(errs) = errors else { return; };

    let mut l_cumul_sq_err = 0.;

    if let Some(last_l_vary_idx) = last_l_vary_idx {
        for (i, comp) in model.lifetime_components.iter_mut().enumerate() {
            let i_idx = 3 + 2 * i + 1;
            match param_types[i_idx] {
                ParamType::Varied(idx) => {
                    comp.intensity.err = absolute_l_intensities[i] * (
                        l_cumul_sq_err +
                        (errs[idx] / (transformed_params[idx] + 1e-12)).powi(2)
                    ).sqrt();
                    l_cumul_sq_err += (
                        errs[idx] / (1. - transformed_params[idx] + 1e-12)
                    ).powi(2);
                },
                ParamType::LastVariedIntensity(_) => {
                    comp.intensity.err = absolute_l_intensities[last_l_vary_idx] * l_cumul_sq_err.sqrt();
                    comp.intensity.vary = true;
                    break;
                },
                ParamType::Fixed(_) => {
                    comp.intensity.err = 0.;
                },
            }
        }
    }

    let mut r_cumul_sq_err = 0.;
    let n_l = model.lifetime_components.len();

    if let Some(last_r_vary_idx) = last_r_vary_idx {
        for (i, comp) in model.resolution_components.iter_mut().enumerate() {
            let i_idx = 3 + 2 * n_l + 3 * i + 1;
            match param_types[i_idx] {
                ParamType::Varied(idx) => {
                    comp.intensity.err = absolute_r_intensities[i] * (
                        r_cumul_sq_err +
                        (errs[idx] / (transformed_params[idx] + 1e-12)).powi(2)
                    ).sqrt();
                    r_cumul_sq_err += (
                        errs[idx] / (1. - transformed_params[idx] + 1e-12)
                    ).powi(2);
                },
                ParamType::LastVariedIntensity(_) => {
                    comp.intensity.err = absolute_r_intensities[last_r_vary_idx] * r_cumul_sq_err.sqrt();
                    comp.intensity.vary = true;
                    break;
                },
                ParamType::Fixed(_) => {},
            }
        }
    }
}

fn post_fit(
    model: &mut LifetimeModel,
    param_types: &Vec<ParamType>,
    transformed_params: &DVector<f64>,
    absolute_l_intensities: &Vec<f64>,
    absolute_r_intensities: &Vec<f64>,
    last_l_vary_idx: Option<usize>,
    last_r_vary_idx: Option<usize>,
    errors: Option<&DVector<f64>>,
) {
    let n_l = model.lifetime_components.len();
    let n_r = model.resolution_components.len();

    for i in 0..n_l {
        model.lifetime_components[i].intensity.val = absolute_l_intensities[i];
    }

    for i in 0..n_r {
        model.resolution_components[i].intensity.val = absolute_r_intensities[i];
    }

    paste_val_and_err(
        &mut model.scale_component, 0, param_types, transformed_params, errors
    );
    paste_val_and_err(
        &mut model.background_component, 1, param_types, transformed_params, errors
    );
    paste_val_and_err(
        &mut model.shift_component, 2, param_types, transformed_params, errors
    );

    for i in 0..n_l {
        let comp = &mut model.lifetime_components[i];
        let idx = 3 + 2 * i;
        paste_val_and_err(
            &mut comp.lifetime, idx, param_types, transformed_params, errors
        );
    }

    for i in 0..n_r {
        let comp = &mut model.resolution_components[i];
        let idx = 3 + 2 * n_l + 3 * i;
        paste_val_and_err(
            &mut comp.fwhm, idx, param_types, transformed_params, errors
        );
        comp.fwhm.val *= FWHM_OVER_SIGMA;
        comp.fwhm.err *= FWHM_OVER_SIGMA;
        paste_val_and_err(
            &mut comp.t0, idx + 2, param_types, transformed_params, errors
        );
    }

    transform_intensity_errors(
        model,
        transformed_params,
        param_types,
        absolute_l_intensities,
        absolute_r_intensities,
        last_l_vary_idx,
        last_r_vary_idx,
        errors,
    );
}

pub struct FitResult {
    pub model: LifetimeModel,
    pub n_dpoints: usize,
    pub fit_status: FitStatus,
    pub n_eval: usize,
    pub red_chi_2: Option<f64>,
    pub cov: Option<DMatrix<f64>>,
}

#[derive(Debug, Error)]
pub enum FitError {
    #[error("Fit failed: {info}")]
    Failure {
        info: String,
    },
}

pub struct FitStatus {
    termination: TerminationReason,
}

impl FitStatus {
    pub fn repr(&self) -> &str {
        match self.termination {
            TerminationReason::User(_) => {
                "The residual or Jacobian computation failed"
            },
            TerminationReason::Numerical(_) => {
                "Encountered NaN or inf"
            },
            TerminationReason::ResidualsZero => {
                "The residuals are literally zero"
            },
            TerminationReason::Orthogonal => {
                "gtol termination criterion fulfilled"
            },
            TerminationReason::Converged { ftol, xtol } => {
                if ftol && !xtol {
                    "ftol termination criterion fulfilled"
                } else if !ftol && xtol {
                    "xtol termination criterion fulfilled"
                } else {
                    "ftol and xtol criterion fulfilled"
                }
            },
            TerminationReason::NoImprovementPossible(_) => {
                "The bound for `ftol`, `xtol` or `gtol` was set so low that \
                the test passed with the machine epsilon but not with the \
                actual bound"
            },
            TerminationReason::LostPatience => {
                "Maximum number of function evaluations was hit"
            },
            TerminationReason::NoParameters => {
                "The number of parameters is zero"
            },
            TerminationReason::NoResiduals => {
                "The number of residuals is zero"
            },
            TerminationReason::WrongDimensions(_) => {
                "The shape of the computed residuals or Jacobian is not correct"
            }
        }
    }
}

impl std::fmt::Display for FitStatus {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.repr())
    }
}

pub fn fit_lifetime_spectrum(
    x: &[f64],
    y: &[f64],
    mut model: LifetimeModel,
    counts: f64,
    peak_center: f64,
) -> Result<FitResult, FitError> {
    let (last_l_vary_idx, last_r_vary_idx) = prepare_fit(
        &mut model, x, counts, peak_center
    );

    let mut problem = build_problem(x, y, &model);

    if let Some(idx) = last_l_vary_idx {
        problem.param_types[3 + 2 * idx + 1] = ParamType::LastVariedIntensity(0);
    }

    if let Some(idx) = last_r_vary_idx {
        problem.param_types[3 + 2 * problem.n_l + 3 * idx + 1] = ParamType::LastVariedIntensity(0);
    }

    problem.set_init();
    problem.make_intensities();

    let lm = LevenbergMarquardt::new()
        .with_ftol(1.49012e-08)
        .with_xtol(1.49012e-08)
        .with_gtol(0.)
        .with_scale_diag(true);

    let (problem, report) = lm.minimize(problem);

    let fit_status = FitStatus { termination: report.termination };

    if !fit_status.termination.was_successful() {
        return Err(FitError::Failure { info: fit_status.repr().into() });
    }

    let red_chi2 = compute_red_chi2(&problem);

    let cov = match red_chi2 {
        None => None,
        Some(chi2) => compute_covariance(&problem, chi2)
    };

    let std_err = match &cov {
        None => None,
        Some(c) => Some(c.diagonal().map(|v| v.sqrt())),
    };

    let param_types = problem.param_types;
    let transformed_params = problem.transformed_params;
    let absolute_l_intensities = problem.absolute_l_intensities;
    let absolute_r_intensities = problem.absolute_r_intensities;

    post_fit(
        &mut model,
        &param_types,
        &transformed_params,
        &absolute_l_intensities,
        &absolute_r_intensities,
        last_l_vary_idx,
        last_r_vary_idx,
        std_err.as_ref(),
    );

    Ok(FitResult {
        model,
        n_dpoints: x.len(),
        fit_status: fit_status,
        n_eval: report.number_of_evaluations,
        red_chi_2: red_chi2,
        cov: cov,
    })
}
