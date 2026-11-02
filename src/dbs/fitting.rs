use levenberg_marquardt::{LeastSquaresProblem, LevenbergMarquardt};
use nalgebra::{DMatrix, DVector, Dyn, Owned};
use statrs::function::erf::erfc;

use crate::core::fitting::{BoundedLeastSquaresProblem, FitParam, FitStatistics, FitStatus, LMFitError, SimpleFitParam};

use crate::constants::TWO_OVER_SQRT_PI;


struct GaussProblem {
    x: DVector<f64>,
    y: DVector<f64>,
    p: DVector<f64>,
    params: Vec<FitParam>,
    transformed_params: DVector<f64>,
    transforms: Vec<fn (f64, &FitParam) -> f64>,
    diffs: Vec<fn (f64, &FitParam) -> f64>,
}


impl LeastSquaresProblem<f64, Dyn, Dyn> for GaussProblem {
    type ParameterStorage = Owned<f64, Dyn>;
    type ResidualStorage = Owned<f64, Dyn>;
    type JacobianStorage = Owned<f64, Dyn, Dyn>;

    fn set_params(&mut self, params: &DVector<f64>) {
        self.p.copy_from(params);
        for i in 0..params.len() {
            self.transformed_params[i] = self.transforms[i](
                params[i], &self.params[i]
            );
        }
    }

    fn params(&self) -> DVector<f64> {
        self.p.clone()
    }

    fn residuals(&self) -> Option<DVector<f64>> {
        let amp = self.transformed_params[0];
        let x0 = self.transformed_params[1];
        let sig = self.transformed_params[2];

        let f = self.x.map(|x| {
            let u = (x - x0) / sig;
            amp * (-0.5 * u*u).exp()
        });

        let r = &self.y - f;

        Some(r)
    }

    fn jacobian(&self) -> Option<DMatrix<f64>> {
        let mut j = DMatrix::zeros(self.y.len(), 3);

        let amp = self.transformed_params[0];
        let x0 = self.transformed_params[1];
        let sig = self.transformed_params[2];

        let diff_amp = self.diffs[0](self.p[0], &self.params[0]);
        let diff_x0 = self.diffs[1](self.p[1], &self.params[1]);
        let diff_sig = self.diffs[2](self.p[2], &self.params[2]);

        for (i, x) in self.x.iter().enumerate() {
            let u = (x - x0) / sig;
            let z = u / sig;

            let mut e = (-0.5 * u*u).exp();

            j[(i, 0)] = e * diff_amp;

            e *= amp;

            j[(i, 1)] = e * z * diff_x0;
            j[(i, 2)] = e * u * z * diff_sig;
        }

        Some(j)
    }
}


impl BoundedLeastSquaresProblem for GaussProblem {
    fn diffs(&self) -> &Vec<fn (f64, &FitParam) -> f64> {
        &self.diffs
    }

    fn param_data(&self, i: usize) -> &FitParam {
        &self.params[i]
    }
}


pub struct GaussParams {
    pub amplitude: SimpleFitParam,
    pub center: SimpleFitParam,
    pub sigma: SimpleFitParam,
}


pub struct GaussFitResult {
    pub params: GaussParams,
    pub stats: FitStatistics,
}


pub fn fit_gauss<'a>(
    x: &[f64], y: &[f64], init: &[f64]
) -> Result<GaussFitResult, LMFitError> {
    let mut params = Vec::new();

    params.push(FitParam::new(init[0]).min(0.));  // amp_1
    params.push(FitParam::new(init[1]));  // x0_1
    params.push(FitParam::new(init[2]).min(0.1));  // sig_1

    let transforms = params.iter().map(|param| param.transform).collect();
    let backtransforms = params.iter().map(|param| param.backtransform).collect::<Vec<_>>();
    let diffs = params.iter().map(|param| param.diff).collect();

    let mut p = Vec::new();

    for i in 0..init.len() {
        p.push(backtransforms[i](init[i], &params[i]));
    }

    let problem = GaussProblem {
        x: DVector::from_column_slice(x),
        y: DVector::from_column_slice(y),
        p: DVector::from_column_slice(&p),
        params,
        transformed_params: DVector::from_column_slice(init),
        transforms,
        diffs,
    };

    let (result, report) = LevenbergMarquardt::new().minimize(problem);

    if !report.termination.was_successful() {
        let fit_status = FitStatus { termination: report.termination };
        return Err(LMFitError::Failure { info: fit_status.repr().into() });
    }

    let stats = result.stats(report).unwrap();

    let errs = stats.err.clone().unwrap_or(DVector::from_element(p.len(), f64::NAN));
    let vals = result.transformed_params;

    let params = GaussParams {
        amplitude: SimpleFitParam { val: vals[0], err: errs[0] },
        center: SimpleFitParam { val: vals[1], err: errs[1] },
        sigma: SimpleFitParam { val: vals[2], err: errs[2] },
    };

    Ok(GaussFitResult { params, stats })
}


struct ErfLinear1GaussProblem {
    x: DVector<f64>,
    y: DVector<f64>,
    w: DVector<f64>,
    p: DVector<f64>,
    params: Vec<FitParam>,
    transformed_params: DVector<f64>,
    transforms: Vec<fn (f64, &FitParam) -> f64>,
    diffs: Vec<fn (f64, &FitParam) -> f64>,
}


impl LeastSquaresProblem<f64, Dyn, Dyn> for ErfLinear1GaussProblem {
    type ParameterStorage = Owned<f64, Dyn>;
    type ResidualStorage = Owned<f64, Dyn>;
    type JacobianStorage = Owned<f64, Dyn, Dyn>;

    fn set_params(&mut self, params: &DVector<f64>) {
        self.p.copy_from(params);
        for i in 0..params.len() {
            self.transformed_params[i] = self.transforms[i](
                params[i], &self.params[i]
            );
        }
    }

    fn params(&self) -> DVector<f64> {
        self.p.clone()
    }

    fn residuals(&self) -> Option<DVector<f64>> {
        let amp_1 = self.transformed_params[0];
        let x0_1 = self.transformed_params[1];
        let sig_1 = self.transformed_params[2];
        let erf_amp = self.transformed_params[3];
        let lin = self.transformed_params[4];
        let off = self.transformed_params[5];

        let f = self.x.map(|x| {
            let u = (x - x0_1) / sig_1;
            amp_1 * (-0.5 * u*u).exp() + erf_amp * erfc(u) + lin * x + off
        });

        let r = (f - &self.y).component_mul(&self.w);

        Some(r)
    }

    fn jacobian(&self) -> Option<DMatrix<f64>> {
        let mut j = DMatrix::zeros(self.y.len(), 6);

        let amp_1 = self.transformed_params[0];
        let x0_1 = self.transformed_params[1];
        let sig_1 = self.transformed_params[2];
        let erf_amp = self.transformed_params[3];

        let diff_amp_1 = self.diffs[0](self.p[0], &self.params[0]);
        let diff_x0_1 = self.diffs[1](self.p[1], &self.params[1]);
        let diff_sig_1 = self.diffs[2](self.p[2], &self.params[2]);
        let diff_erf_amp = self.diffs[3](self.p[3], &self.params[3]);
        let diff_lin = self.diffs[4](self.p[4], &self.params[4]);
        let diff_off = self.diffs[5](self.p[5], &self.params[5]);

        for (i, x) in self.x.iter().enumerate() {
            let u_1 = (x - x0_1) / sig_1;
            let z_1 = u_1 / sig_1;
            let h_1 = erf_amp * TWO_OVER_SQRT_PI * (-u_1*u_1).exp() / sig_1;

            let mut e_1 = (-0.5 * u_1*u_1).exp();

            j[(i, 0)] = e_1 * diff_amp_1;

            e_1 *= amp_1;

            j[(i, 1)] = (e_1 * z_1 + h_1) * diff_x0_1;
            j[(i, 2)] = (e_1 * u_1 * z_1 + h_1 * u_1) * diff_sig_1;
            j[(i, 3)] = erfc(u_1) * diff_erf_amp;
            j[(i, 4)] = *x * diff_lin;
            j[(i, 5)] = diff_off;
        }

        for mut col in j.column_iter_mut() {
            col.component_mul_assign(&self.w);
        }

        Some(j)
    }
}


impl BoundedLeastSquaresProblem for ErfLinear1GaussProblem {
    fn diffs(&self) -> &Vec<fn (f64, &FitParam) -> f64> {
        &self.diffs
    }

    fn param_data(&self, i: usize) -> &FitParam {
        &self.params[i]
    }
}


pub struct ErfLinear1GaussParams {
    pub gauss_amplitude: SimpleFitParam,
    pub center: SimpleFitParam,
    pub sigma: SimpleFitParam,
    pub erf_amplitude: SimpleFitParam,
    pub bg_linear_coeff: SimpleFitParam,
    pub bg_const_coeff: SimpleFitParam,
}


pub struct ErfLinear1GaussFitResult {
    pub params: ErfLinear1GaussParams,
    pub stats: FitStatistics,
}


pub fn fit_erf_linear_1_gauss(
    x: &[f64], y: &[f64], init: &[f64]
) -> Result<ErfLinear1GaussFitResult, LMFitError> {
    let mut params = Vec::new();

    params.push(FitParam::new(init[0]).min(0.));  // amp_1
    params.push(FitParam::new(init[1]));  // x0_1
    params.push(FitParam::new(init[2]).min(0.1));  // sig_1
    params.push(FitParam::new(init[3]).min(0.));  // erf_amp
    params.push(FitParam::new(init[4]));  // lin
    params.push(FitParam::new(init[5]));  // off

    let transforms = params.iter().map(|param| param.transform).collect();
    let backtransforms = params.iter().map(|param| param.backtransform).collect::<Vec<_>>();
    let diffs = params.iter().map(|param| param.diff).collect();

    let mut p = Vec::new();

    for i in 0..init.len() {
        p.push(backtransforms[i](init[i], &params[i]));
    }

    let problem = ErfLinear1GaussProblem {
        x: DVector::from_column_slice(x),
        y: DVector::from_column_slice(y),
        w: DVector::from_iterator(y.len(), y.iter().map(|y_i| 1. / y_i.max(1.).sqrt())),
        p: DVector::from_column_slice(init),
        params,
        transformed_params: DVector::from_column_slice(init),
        transforms,
        diffs,
    };

    let (result, report) = LevenbergMarquardt::new().minimize(problem);

    if !report.termination.was_successful() {
        let fit_status = FitStatus { termination: report.termination };
        return Err(LMFitError::Failure { info: fit_status.repr().into() });
    }

    let stats = result.stats(report).unwrap();

    let errs = stats.err.clone().unwrap_or(DVector::from_element(p.len(), f64::NAN));
    let vals = result.transformed_params;

    let params = ErfLinear1GaussParams {
        gauss_amplitude: SimpleFitParam { val: vals[0], err: errs[0] },
        center: SimpleFitParam { val: vals[1], err: errs[1] },
        sigma: SimpleFitParam { val: vals[2], err: errs[2] },
        erf_amplitude: SimpleFitParam { val: vals[3], err: errs[3] },
        bg_linear_coeff: SimpleFitParam { val: vals[4], err: errs[4] },
        bg_const_coeff: SimpleFitParam { val: vals[5], err: errs[5] },
    };

    Ok(ErfLinear1GaussFitResult { params, stats })
}


struct ErfLinear2GaussProblem {
    x: DVector<f64>,
    y: DVector<f64>,
    w: DVector<f64>,
    p: DVector<f64>,
    params: Vec<FitParam>,
    transformed_params: DVector<f64>,
    transforms: Vec<fn (f64, &FitParam) -> f64>,
    diffs: Vec<fn (f64, &FitParam) -> f64>,
}


impl LeastSquaresProblem<f64, Dyn, Dyn> for ErfLinear2GaussProblem {
    type ParameterStorage = Owned<f64, Dyn>;
    type ResidualStorage = Owned<f64, Dyn>;
    type JacobianStorage = Owned<f64, Dyn, Dyn>;

    fn set_params(&mut self, params: &DVector<f64>) {
        self.p.copy_from(params);
        for i in 0..params.len() {
            self.transformed_params[i] = self.transforms[i](
                params[i], &self.params[i]
            );
        }
    }

    fn params(&self) -> DVector<f64> {
        self.p.clone()
    }

    fn residuals(&self) -> Option<DVector<f64>> {
        let amp_1 = self.transformed_params[0];
        let x0_1 = self.transformed_params[1];
        let sig_1 = self.transformed_params[2];
        let amp_2 = self.transformed_params[3];
        let x0_2 = self.transformed_params[4];
        let sig_2 = self.transformed_params[5];
        let erf_amp = self.transformed_params[6];
        let lin = self.transformed_params[7];
        let off = self.transformed_params[8];

        let f = self.x.map(|x| {
            let u_1 = (x - x0_1) / sig_1;
            let u_2 = (x - x0_2) / sig_2;
            amp_1 * (-0.5 * u_1*u_1).exp()
                + amp_2 * (-0.5 * u_2*u_2).exp()
                + erf_amp * erfc(u_1) + lin * x + off
        });

        let r = (f - &self.y).component_mul(&self.w);

        Some(r)
    }

    fn jacobian(&self) -> Option<DMatrix<f64>> {
        let mut j = DMatrix::zeros(self.y.len(), 9);

        let amp_1 = self.transformed_params[0];
        let x0_1 = self.transformed_params[1];
        let sig_1 = self.transformed_params[2];
        let amp_2 = self.transformed_params[3];
        let x0_2 = self.transformed_params[4];
        let sig_2 = self.transformed_params[5];
        let erf_amp = self.transformed_params[6];

        let diff_amp_1 = self.diffs[0](self.p[0], &self.params[0]);
        let diff_x0_1 = self.diffs[1](self.p[1], &self.params[1]);
        let diff_sig_1 = self.diffs[2](self.p[2], &self.params[2]);
        let diff_amp_2 = self.diffs[3](self.p[3], &self.params[3]);
        let diff_x0_2 = self.diffs[4](self.p[4], &self.params[4]);
        let diff_sig_2 = self.diffs[5](self.p[5], &self.params[5]);
        let diff_erf_amp = self.diffs[6](self.p[6], &self.params[6]);
        let diff_lin = self.diffs[7](self.p[7], &self.params[7]);
        let diff_off = self.diffs[8](self.p[8], &self.params[8]);

        for (i, x) in self.x.iter().enumerate() {
            let u_1 = (x - x0_1) / sig_1;
            let z_1 = u_1 / sig_1;
            let h_1 = erf_amp * TWO_OVER_SQRT_PI * (-u_1*u_1).exp() / sig_1;
            let u_2 = (x - x0_2) / sig_2;
            let z_2 = u_2 / sig_2;

            let mut e_1 = (-0.5 * u_1*u_1).exp();
            let mut e_2 = (-0.5 * u_2*u_2).exp();

            j[(i, 0)] = e_1 * diff_amp_1;
            j[(i, 3)] = e_2 * diff_amp_2;

            e_1 *= amp_1;
            e_2 *= amp_2;

            j[(i, 1)] = (e_1 * z_1 + h_1) * diff_x0_1;
            j[(i, 2)] = (e_1 * u_1 * z_1 + h_1 * u_1) * diff_sig_1;
            j[(i, 4)] = e_2 * z_2 * diff_x0_2;
            j[(i, 5)] = e_2 * u_2 * z_2 * diff_sig_2;
            j[(i, 6)] = erfc(u_1) * diff_erf_amp;
            j[(i, 7)] = *x * diff_lin;
            j[(i, 8)] = diff_off;
        }

        for mut col in j.column_iter_mut() {
            col.component_mul_assign(&self.w);
        }

        Some(j)
    }
}


impl BoundedLeastSquaresProblem for ErfLinear2GaussProblem {
    fn diffs(&self) -> &Vec<fn (f64, &FitParam) -> f64> {
        &self.diffs
    }

    fn param_data(&self, i: usize) -> &FitParam {
        &self.params[i]
    }
}


pub struct ErfLinear2GaussParams {
    pub gauss_1_amplitude: SimpleFitParam,
    pub gauss_1_center: SimpleFitParam,
    pub gauss_1_sigma: SimpleFitParam,
    pub gauss_2_amplitude: SimpleFitParam,
    pub gauss_2_center: SimpleFitParam,
    pub gauss_2_sigma: SimpleFitParam,
    pub erf_amplitude: SimpleFitParam,
    pub bg_linear_coeff: SimpleFitParam,
    pub bg_const_coeff: SimpleFitParam,
}


pub struct ErfLinear2GaussFitResult {
    pub params: ErfLinear2GaussParams,
    pub stats: FitStatistics,
}


pub fn fit_erf_linear_2_gauss(
    x: &[f64], y: &[f64], init: &[f64]
) -> Result<ErfLinear2GaussFitResult, LMFitError> {
    let mut params = Vec::new();

    params.push(FitParam::new(init[0]).min(0.));  // amp_1
    params.push(FitParam::new(init[1]));  // x0_1
    params.push(FitParam::new(init[2]).min(0.1));  // sig_1
    params.push(FitParam::new(init[3]).min(0.));  // amp_2
    params.push(FitParam::new(init[4]));  // x0_2
    params.push(FitParam::new(init[5]).min(0.1));  // sig_2
    params.push(FitParam::new(init[6]).min(0.));  // erf_amp
    params.push(FitParam::new(init[7]));  // lin
    params.push(FitParam::new(init[8]));  // off

    let transforms = params.iter().map(|param| param.transform).collect();
    let backtransforms = params.iter().map(|param| param.backtransform).collect::<Vec<_>>();
    let diffs = params.iter().map(|param| param.diff).collect();

    let mut p = Vec::new();

    for i in 0..init.len() {
        p.push(backtransforms[i](init[i], &params[i]));
    }

    let problem = ErfLinear2GaussProblem {
        x: DVector::from_column_slice(x),
        y: DVector::from_column_slice(y),
        w: DVector::from_iterator(y.len(), y.iter().map(|y_i| 1. / y_i.max(1.).sqrt())),
        p: DVector::from_column_slice(init),
        params,
        transformed_params: DVector::from_column_slice(init),
        transforms,
        diffs,
    };

    let (result, report) = LevenbergMarquardt::new().with_patience(1000).minimize(problem);

    if !report.termination.was_successful() {
        let fit_status = FitStatus { termination: report.termination };
        return Err(LMFitError::Failure { info: fit_status.repr().into() });
    }

    let stats = result.stats(report).unwrap();

    let errs = stats.err.clone().unwrap_or(DVector::from_element(p.len(), f64::NAN));
    let vals = result.transformed_params;

    let params = ErfLinear2GaussParams {
        gauss_1_amplitude: SimpleFitParam { val: vals[0], err: errs[0] },
        gauss_1_center: SimpleFitParam { val: vals[1], err: errs[1] },
        gauss_1_sigma: SimpleFitParam { val: vals[2], err: errs[2] },
        gauss_2_amplitude: SimpleFitParam { val: vals[3], err: errs[3] },
        gauss_2_center: SimpleFitParam { val: vals[4], err: errs[4] },
        gauss_2_sigma: SimpleFitParam { val: vals[5], err: errs[5] },
        erf_amplitude: SimpleFitParam { val: vals[6], err: errs[6] },
        bg_linear_coeff: SimpleFitParam { val: vals[7], err: errs[7] },
        bg_const_coeff: SimpleFitParam { val: vals[8], err: errs[8] },
    };

    Ok(ErfLinear2GaussFitResult { params, stats })
}


struct ErfLinear3GaussProblem {
    x: DVector<f64>,
    y: DVector<f64>,
    w: DVector<f64>,
    p: DVector<f64>,
    params: Vec<FitParam>,
    transformed_params: DVector<f64>,
    transforms: Vec<fn (f64, &FitParam) -> f64>,
    diffs: Vec<fn (f64, &FitParam) -> f64>,
}


impl LeastSquaresProblem<f64, Dyn, Dyn> for ErfLinear3GaussProblem {
    type ParameterStorage = Owned<f64, Dyn>;
    type ResidualStorage = Owned<f64, Dyn>;
    type JacobianStorage = Owned<f64, Dyn, Dyn>;

    fn set_params(&mut self, params: &DVector<f64>) {
        self.p.copy_from(params);
        for i in 0..params.len() {
            self.transformed_params[i] = self.transforms[i](
                params[i], &self.params[i]
            );
        }
    }

    fn params(&self) -> DVector<f64> {
        self.p.clone()
    }

    fn residuals(&self) -> Option<DVector<f64>> {
        let amp_1 = self.transformed_params[0];
        let x0_1 = self.transformed_params[1];
        let sig_1 = self.transformed_params[2];
        let amp_2 = self.transformed_params[3];
        let x0_2 = self.transformed_params[4];
        let sig_2 = self.transformed_params[5];
        let amp_3 = self.transformed_params[6];
        let x0_3 = self.transformed_params[7];
        let sig_3 = self.transformed_params[8];
        let erf_amp = self.transformed_params[9];
        let lin = self.transformed_params[10];
        let off = self.transformed_params[11];

        let f = self.x.map(|x| {
            let u_1 = (x - x0_1) / sig_1;
            let u_2 = (x - x0_2) / sig_2;
            let u_3 = (x - x0_3) / sig_3;
            amp_1 * (-0.5 * u_1*u_1).exp()
                + amp_2 * (-0.5 * u_2*u_2).exp()
                + amp_3 * (-0.5 * u_3*u_3).exp()
                + erf_amp * erfc(u_1) + lin * x + off
        });

        let r = (f - &self.y).component_mul(&self.w);

        Some(r)
    }

    fn jacobian(&self) -> Option<DMatrix<f64>> {
        let mut j = DMatrix::zeros(self.y.len(), 12);

        let amp_1 = self.transformed_params[0];
        let x0_1 = self.transformed_params[1];
        let sig_1 = self.transformed_params[2];
        let amp_2 = self.transformed_params[3];
        let x0_2 = self.transformed_params[4];
        let sig_2 = self.transformed_params[5];
        let amp_3 = self.transformed_params[6];
        let x0_3 = self.transformed_params[7];
        let sig_3 = self.transformed_params[8];
        let erf_amp = self.transformed_params[9];

        let diff_amp_1 = self.diffs[0](self.p[0], &self.params[0]);
        let diff_x0_1 = self.diffs[1](self.p[1], &self.params[1]);
        let diff_sig_1 = self.diffs[2](self.p[2], &self.params[2]);
        let diff_amp_2 = self.diffs[3](self.p[3], &self.params[3]);
        let diff_x0_2 = self.diffs[4](self.p[4], &self.params[4]);
        let diff_sig_2 = self.diffs[5](self.p[5], &self.params[5]);
        let diff_amp_3 = self.diffs[6](self.p[6], &self.params[6]);
        let diff_x0_3 = self.diffs[7](self.p[7], &self.params[7]);
        let diff_sig_3 = self.diffs[8](self.p[8], &self.params[8]);
        let diff_erf_amp = self.diffs[9](self.p[9], &self.params[9]);
        let diff_lin = self.diffs[10](self.p[10], &self.params[10]);
        let diff_off = self.diffs[11](self.p[11], &self.params[11]);

        for (i, x) in self.x.iter().enumerate() {
            let u_1 = (x - x0_1) / sig_1;
            let z_1 = u_1 / sig_1;
            let h_1 = erf_amp * TWO_OVER_SQRT_PI * (-u_1*u_1).exp() / sig_1;
            let u_2 = (x - x0_2) / sig_2;
            let z_2 = u_2 / sig_2;
            let u_3 = (x - x0_3) / sig_3;
            let z_3 = u_3 / sig_3;

            let mut e_1 = (-0.5 * u_1*u_1).exp();
            let mut e_2 = (-0.5 * u_2*u_2).exp();
            let mut e_3 = (-0.5 * u_3*u_3).exp();

            j[(i, 0)] = e_1 * diff_amp_1;
            j[(i, 3)] = e_2 * diff_amp_2;
            j[(i, 6)] = e_3 * diff_amp_3;

            e_1 *= amp_1;
            e_2 *= amp_2;
            e_3 *= amp_3;

            j[(i, 1)] = (e_1 * z_1 + h_1) * diff_x0_1;
            j[(i, 2)] = (e_1 * u_1 * z_1 + h_1 * u_1) * diff_sig_1;
            j[(i, 4)] = e_2 * z_2 * diff_x0_2;
            j[(i, 5)] = e_2 * u_2 * z_2 * diff_sig_2;
            j[(i, 7)] = e_3 * z_3 * diff_x0_3;
            j[(i, 8)] = e_3 * u_3 * z_3 * diff_sig_3;
            j[(i, 9)] = erfc(u_1) * diff_erf_amp;
            j[(i, 10)] = *x * diff_lin;
            j[(i, 11)] = diff_off;
        }

        for mut col in j.column_iter_mut() {
            col.component_mul_assign(&self.w);
        }

        Some(j)
    }
}


impl BoundedLeastSquaresProblem for ErfLinear3GaussProblem {
    fn diffs(&self) -> &Vec<fn (f64, &FitParam) -> f64> {
        &self.diffs
    }

    fn param_data(&self, i: usize) -> &FitParam {
        &self.params[i]
    }
}


pub struct ErfLinear3GaussParams {
    pub gauss_1_amplitude: SimpleFitParam,
    pub gauss_1_center: SimpleFitParam,
    pub gauss_1_sigma: SimpleFitParam,
    pub gauss_2_amplitude: SimpleFitParam,
    pub gauss_2_center: SimpleFitParam,
    pub gauss_2_sigma: SimpleFitParam,
    pub gauss_3_amplitude: SimpleFitParam,
    pub gauss_3_center: SimpleFitParam,
    pub gauss_3_sigma: SimpleFitParam,
    pub erf_amplitude: SimpleFitParam,
    pub bg_linear_coeff: SimpleFitParam,
    pub bg_const_coeff: SimpleFitParam,
}


pub struct ErfLinear3GaussFitResult {
    pub params: ErfLinear3GaussParams,
    pub stats: FitStatistics,
}


pub fn fit_erf_linear_3_gauss(
    x: &[f64], y: &[f64], init: &[f64]
) -> Result<ErfLinear3GaussFitResult, LMFitError> {
    let mut params = Vec::new();

    params.push(FitParam::new(init[0]).min(0.));  // amp_1
    params.push(FitParam::new(init[1]));  // x0_1
    params.push(FitParam::new(init[2]).min(0.1));  // sig_1
    params.push(FitParam::new(init[3]).min(0.));  // amp_2
    params.push(FitParam::new(init[4]));  // x0_2
    params.push(FitParam::new(init[5]).min(0.1));  // sig_2
    params.push(FitParam::new(init[6]).min(0.));  // amp_3
    params.push(FitParam::new(init[7]));  // x0_3
    params.push(FitParam::new(init[8]).min(0.1));  // sig_3
    params.push(FitParam::new(init[9]).min(0.));  // erf_amp
    params.push(FitParam::new(init[10]));  // lin
    params.push(FitParam::new(init[11]));  // off

    let transforms = params.iter().map(|param| param.transform).collect();
    let backtransforms = params.iter().map(|param| param.backtransform).collect::<Vec<_>>();
    let diffs = params.iter().map(|param| param.diff).collect();

    let mut p = Vec::new();

    for i in 0..init.len() {
        p.push(backtransforms[i](init[i], &params[i]));
    }

    let problem = ErfLinear3GaussProblem {
        x: DVector::from_column_slice(x),
        y: DVector::from_column_slice(y),
        w: DVector::from_iterator(y.len(), y.iter().map(|y_i| 1. / y_i.max(1.).sqrt())),
        p: DVector::from_column_slice(init),
        params,
        transformed_params: DVector::from_column_slice(init),
        transforms,
        diffs,
    };

    let (result, report) = LevenbergMarquardt::new().minimize(problem);

    if !report.termination.was_successful() {
        let fit_status = FitStatus { termination: report.termination };
        return Err(LMFitError::Failure { info: fit_status.repr().into() });
    }

    let stats = result.stats(report).unwrap();

    let errs = stats.err.clone().unwrap_or(DVector::from_element(p.len(), f64::NAN));
    let vals = result.transformed_params;

    let params = ErfLinear3GaussParams {
        gauss_1_amplitude: SimpleFitParam { val: vals[0], err: errs[0] },
        gauss_1_center: SimpleFitParam { val: vals[1], err: errs[1] },
        gauss_1_sigma: SimpleFitParam { val: vals[2], err: errs[2] },
        gauss_2_amplitude: SimpleFitParam { val: vals[3], err: errs[3] },
        gauss_2_center: SimpleFitParam { val: vals[4], err: errs[4] },
        gauss_2_sigma: SimpleFitParam { val: vals[5], err: errs[5] },
        gauss_3_amplitude: SimpleFitParam { val: vals[6], err: errs[6] },
        gauss_3_center: SimpleFitParam { val: vals[7], err: errs[7] },
        gauss_3_sigma: SimpleFitParam { val: vals[8], err: errs[8] },
        erf_amplitude: SimpleFitParam { val: vals[9], err: errs[9] },
        bg_linear_coeff: SimpleFitParam { val: vals[10], err: errs[10] },
        bg_const_coeff: SimpleFitParam { val: vals[11], err: errs[11] },
    };

    Ok(ErfLinear3GaussFitResult { params, stats })
}


pub fn erf_linear_background(
    x: &[f64], amp: f64, x0: f64, sig: f64, lin: f64, off: f64
) -> Vec<f64> {
    x.iter().map(|x| {
        let u = (x - x0) / sig;
        amp * erfc(u) + lin * x + off
    }).collect()
}
