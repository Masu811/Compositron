use std::collections::HashMap;
use std::f64::{self, consts::SQRT_2};
use std::{fs::File, io:: {BufWriter, Write}};

use approx::assert_relative_eq;
use levenberg_marquardt::{differentiate_numerically, LeastSquaresProblem, LevenbergMarquardt};
use nalgebra::{DMatrix, DVector, Dyn, Owned};
use statrs::consts::SQRT_2PI;
use statrs::function::erf::erfc;
use thiserror::Error;

use crate::{core::fitting::FitParam, pals::model::LifetimeModel};

#[derive(Debug, Error)]
pub enum FitError {

}

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
            println!("NaN values in model function detected");
            println!("internal params: {:?}", self.p);
            println!("transformed params: {:?}", self.transformed_params);
            println!("l intensities: {:?}", self.absolute_l_intensities);
            println!("r intensities: {:?}", self.absolute_r_intensities);
            panic!();
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
                    n * 0.5 * l_int * r_int * one_over_tau_sqrt_2_pi
                        * (a[i] - b[i].powi(2)).exp()
                });

                terms.set_column(j + i * n_r, &decay);

                f += &decay;

                let df_dt0 = -&exp_comb / sig + &decay / tau;
                let df_dtau = &exp_comb * sig / tau.powi(2)
                    + &decay.component_mul(&(&tt / tau.powi(2)).add_scalar(
                        sig.powi(2) / tau.powi(3) - one_over_tau
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
            j.set_column(
                i,
                &(
                    // -derivatives.get(param.into()).unwrap().component_mul(
                    //     &(&self.w * self.diffs[i](self.p[i], self.params[i]))
                    // )
                    -derivatives.get(param.into()).unwrap().component_mul(
                        &(&self.w * self.diffs[i](self.p[i], self.params[i]))
                    )
                )
            );

            if j.iter().any(|x| x.is_nan()) {
                println!("NaN values in Jacobian detected");
                println!("Concerns parameter {}", self.varied_param_names[i]);
                println!("internal params: {:?}", self.p);
                println!("transformed params: {:?}", self.transformed_params);
                println!("l intensities: {:?}", self.absolute_l_intensities);
                println!("r intensities: {:?}", self.absolute_r_intensities);
                panic!();
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
    problem: &impl levenberg_marquardt::LeastSquaresProblem<f64, Dyn, Dyn>,
    red_chi2: f64,
) -> Option<DMatrix<f64>> {
    let j = problem.jacobian()?;

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

fn gobble<'a>(
    name: String,
    param: &'a FitParam,
    tmp: &mut ProblemTemplate<'a>,
    is_l_intensity: bool,
    is_r_intensity: bool,
) {
    if param.vary {
        tmp.transformed_params.push(param.val);
        tmp.param_types.push(ParamType::Varied(tmp.transformed_params.len() - 1));
        tmp.varied_param_names.push(name);
        tmp.params.push(param);
        tmp.transforms.push(param.transform);
        tmp.backtransforms.push(param.backtransform);
        tmp.diffs.push(param.diff);
    } else {
        tmp.fixed_params.push(param.val);
        tmp.param_types.push(ParamType::Fixed(tmp.fixed_params.len() - 1));
        if is_l_intensity {
            tmp.available_l_intensity -= param.val;
        } else if is_r_intensity {
            tmp.available_r_intensity -= param.val;
        }
    }
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

    gobble("N".into(), &model.scale_component, &mut tmp, false, false);
    gobble("background".into(), &model.background_component, &mut tmp, false, false);
    gobble("t0".into(), &model.shift_component, &mut tmp, false, false);

    for i in 0..n_l {
        let c = &model.lifetime_components[i];
        gobble(format!("lifetime_{}", i + 1), &c.lifetime, &mut tmp, false, false);
        gobble(format!("intensity_{}", i + 1), &c.intensity, &mut tmp, true, false);
    }

    for i in 0..n_r {
        let c = &model.resolution_components[i];
        gobble(format!("res_fwhm_{}", i + 1), &c.fwhm, &mut tmp, false, false);
        gobble(format!("res_intensity_{}", i + 1), &c.intensity, &mut tmp, false, true);
        gobble(format!("res_t0_{}", i + 1), &c.t0, &mut tmp, false, false);
    }

    let mut problem = Problem {
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

    problem.set_init();
    problem.make_intensities();

    problem
}

fn prepare_fit(model: LifetimeModel) -> LifetimeModel {
    // Make initial guesses for missing parameters
    //
    // Slap on necessary boundaries
    //
    // Set Last Varied Intensity
    model
}

fn post_fit(model: LifetimeModel) -> LifetimeModel {
    // Transform intensities
    //
    // Put values and errors into model
    model
}

pub struct FitResult {
    model: LifetimeModel,
    termination: levenberg_marquardt::TerminationReason,
    n_eval: usize,
    chi_2: Option<f64>,
    cov: Option<DMatrix<f64>>,
}

pub fn fit_lifetime_spectrum(
    x: &[f64], y: &[f64], model: LifetimeModel
) -> Result<FitResult, FitError> {
    let mut problem = build_problem(x, y, &model);

    let jacobian_trait = problem.jacobian().unwrap();
    let jacobian_numerical = differentiate_numerically(&mut problem).unwrap();

    let nrows = jacobian_numerical.nrows();
    let ncols = jacobian_numerical.ncols();

    // let diff = jacobian_numerical - jacobian_trait;

    let mut file = BufWriter::new(File::create("num.csv").unwrap());

    for i in 0..nrows {
        for j in 0..ncols {
            if j > 0 {
                write!(file, ",").unwrap();
            }
            write!(file, "{}", jacobian_numerical[(i, j)]).unwrap();
        }
        writeln!(file).unwrap();
    }

    /////

    let mut file = BufWriter::new(File::create("ana.csv").unwrap());

    for i in 0..nrows {
        for j in 0..ncols {
            if j > 0 {
                write!(file, ",").unwrap();
            }
            write!(file, "{}", jacobian_trait[(i, j)]).unwrap();
        }
        writeln!(file).unwrap();
    }

    panic!();

    let (problem, report) = LevenbergMarquardt::new().minimize(problem);

    // TODO: if not converged, return err

    let red_chi2 = compute_red_chi2(&problem);

    let cov = match red_chi2 {
        None => None,
        Some(chi2) => compute_covariance(&problem, chi2)
    };

    let std_err = match cov {
        None => None,
        Some(c) => Some(c.diagonal().map(|v| v.sqrt())),
    };

    // if let Some(c) = std_err {
    //     println!("std err: {}", c);
    // }

    println!("Converged: {:?}", report.termination);
    println!("N Iter:    {}", report.number_of_evaluations);
    println!("n:       {}", problem.get_param(0));
    println!("bg:      {}", problem.get_param(1));
    println!("t0:      {}", problem.get_param(2));
    for i in 0..problem.n_l {
        println!("tau_{}:   {}", i + 1, problem.get_param(3 + 2 * i));
        println!("int_{}:   {}", i + 1, problem.absolute_l_intensities[i]);
    }
    for i in 0..problem.n_r {
        println!("fwhm_{}:  {}", i + 1, problem.get_param(3 + 2 * problem.n_l + 3 * i));
        println!("r_int_{}: {}", i + 1, problem.get_param(3 + 2 * problem.n_l + 3 * i + 1));
        println!("r_t0_{}:  {}", i + 1, problem.absolute_r_intensities[i]);
    }

    Ok(FitResult {
        model,
        termination: report.termination,
        n_eval: report.number_of_evaluations,
        chi_2: red_chi2,
        cov: cov,
    })
}
