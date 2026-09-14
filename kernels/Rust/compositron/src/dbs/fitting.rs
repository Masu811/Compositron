use std::collections::HashMap;

use levenberg_marquardt::{LeastSquaresProblem, LevenbergMarquardt};
use nalgebra::{DMatrix, DVector, Dyn, Owned};
use statrs::function::erf::erfc;

use crate::core::fitting::{FitStatus, LMFitError, SimpleFitParam};

use crate::constants::TWO_OVER_SQRT_PI;

struct GaussProblem {
    x: DVector<f64>,
    y: DVector<f64>,
    p: DVector<f64>,
}

impl LeastSquaresProblem<f64, Dyn, Dyn> for GaussProblem {
    type ParameterStorage = Owned<f64, Dyn>;
    type ResidualStorage = Owned<f64, Dyn>;
    type JacobianStorage = Owned<f64, Dyn, Dyn>;

    fn set_params(&mut self, params: &DVector<f64>) {
        self.p.copy_from(params);
    }

    fn params(&self) -> DVector<f64> {
        self.p.clone()
    }

    fn residuals(&self) -> Option<DVector<f64>> {
        let amp = self.p[0];
        let x0 = self.p[1];
        let sig = self.p[2];

        let f = self.x.map(|x| {
            let u = (x - x0) / sig;
            amp * (-0.5 * u*u).exp()
        });

        let r = &self.y - f;

        Some(r)
    }

    fn jacobian(&self) -> Option<DMatrix<f64>> {
        let mut j = DMatrix::zeros(self.y.len(), 3);

        let amp = self.p[0];
        let x0 = self.p[1];
        let sig = self.p[2];

        for (i, x) in self.x.iter().enumerate() {
            let u = (x - x0) / sig;
            let z = u / sig;

            let mut e = (-0.5 * u*u).exp();

            j[(i, 0)] = e;

            e *= amp;

            j[(i, 1)] = e * z;
            j[(i, 2)] = e * u * z;
        }

        Some(j)
    }
}


pub fn fit_gauss(
    x: &[f64], y: &[f64], init: &[f64]
) -> Result<HashMap<&'static str, SimpleFitParam>, LMFitError> {
    let problem = GaussProblem {
        x: DVector::from_column_slice(x),
        y: DVector::from_column_slice(y),
        p: DVector::from_column_slice(init),
    };

    let (result, report) = LevenbergMarquardt::new().minimize(problem);

    let fit_status = FitStatus { termination: report.termination };

    if !fit_status.termination.was_successful() {
        return Err(LMFitError::Failure { info: fit_status.repr().into() });
    }

    let mut params = HashMap::new();

    params.insert("amp_1", SimpleFitParam { val: result.p[0], err: f64::NAN });
    params.insert("x0_1", SimpleFitParam { val: result.p[1], err: f64::NAN });
    params.insert("sig_1", SimpleFitParam { val: result.p[2], err: f64::NAN });

    Ok(params)
}


struct ErfLinear1GaussProblem {
    x: DVector<f64>,
    y: DVector<f64>,
    p: DVector<f64>,
}

impl LeastSquaresProblem<f64, Dyn, Dyn> for ErfLinear1GaussProblem {
    type ParameterStorage = Owned<f64, Dyn>;
    type ResidualStorage = Owned<f64, Dyn>;
    type JacobianStorage = Owned<f64, Dyn, Dyn>;

    fn set_params(&mut self, params: &DVector<f64>) {
        self.p.copy_from(params);
    }

    fn params(&self) -> DVector<f64> {
        self.p.clone()
    }

    fn residuals(&self) -> Option<DVector<f64>> {
        let amp_1 = self.p[0];
        let x0_1 = self.p[1];
        let sig_1 = self.p[2];
        let erf_amp = self.p[3];
        let lin = self.p[4];
        let off = self.p[5];

        let f = self.x.map(|x| {
            let u = (x - x0_1) / sig_1;
            amp_1 * (-0.5 * u*u).exp() + erf_amp * erfc(u) + lin * x + off
        });

        let r = &self.y - f;

        Some(r)
    }

    fn jacobian(&self) -> Option<DMatrix<f64>> {
        let mut j = DMatrix::zeros(self.y.len(), 6);

        let amp_1 = self.p[0];
        let x0_1 = self.p[1];
        let sig_1 = self.p[2];
        let erf_amp = self.p[3];

        for (i, x) in self.x.iter().enumerate() {
            let u_1 = (x - x0_1) / sig_1;
            let z_1 = u_1 / sig_1;
            let h_1 = erf_amp * TWO_OVER_SQRT_PI * (-u_1*u_1).exp() / sig_1;

            let mut e_1 = (-0.5 * u_1*u_1).exp();

            j[(i, 0)] = e_1;

            e_1 *= amp_1;

            j[(i, 1)] = e_1 * z_1 + h_1;
            j[(i, 2)] = e_1 * u_1 * z_1 + h_1 * u_1;
            j[(i, 3)] = erfc(u_1);
            j[(i, 4)] = *x;
            j[(i, 5)] = 1.;
        }

        Some(j)
    }
}


pub fn fit_erf_linear_1_gauss(
    x: &[f64], y: &[f64], init: &[f64]
) -> Result<HashMap<&'static str, SimpleFitParam>, LMFitError> {
    let problem = ErfLinear1GaussProblem {
        x: DVector::from_column_slice(x),
        y: DVector::from_column_slice(y),
        p: DVector::from_column_slice(init),
    };

    let (result, report) = LevenbergMarquardt::new().minimize(problem);

    let fit_status = FitStatus { termination: report.termination };

    if !fit_status.termination.was_successful() {
        return Err(LMFitError::Failure { info: fit_status.repr().into() });
    }

    let mut params = HashMap::new();

    params.insert("amp_1", SimpleFitParam { val: result.p[0], err: f64::NAN });
    params.insert("x0_1", SimpleFitParam { val: result.p[1], err: f64::NAN });
    params.insert("sig_1", SimpleFitParam { val: result.p[2], err: f64::NAN });
    params.insert("erf_amp", SimpleFitParam { val: result.p[3], err: f64::NAN });
    params.insert("lin", SimpleFitParam { val: result.p[4], err: f64::NAN });
    params.insert("const", SimpleFitParam { val: result.p[5], err: f64::NAN });

    Ok(params)
}


struct ErfLinear2GaussProblem {
    x: DVector<f64>,
    y: DVector<f64>,
    p: DVector<f64>,
}

impl LeastSquaresProblem<f64, Dyn, Dyn> for ErfLinear2GaussProblem {
    type ParameterStorage = Owned<f64, Dyn>;
    type ResidualStorage = Owned<f64, Dyn>;
    type JacobianStorage = Owned<f64, Dyn, Dyn>;

    fn set_params(&mut self, params: &DVector<f64>) {
        self.p.copy_from(params);
    }

    fn params(&self) -> DVector<f64> {
        self.p.clone()
    }

    fn residuals(&self) -> Option<DVector<f64>> {
        let amp_1 = self.p[0];
        let x0_1 = self.p[1];
        let sig_1 = self.p[2];
        let amp_2 = self.p[3];
        let x0_2 = self.p[4];
        let sig_2 = self.p[5];
        let erf_amp = self.p[6];
        let lin = self.p[7];
        let off = self.p[8];

        let f = self.x.map(|x| {
            let u_1 = (x - x0_1) / sig_1;
            let u_2 = (x - x0_2) / sig_2;
            amp_1 * (-0.5 * u_1*u_1).exp()
                + amp_2 * (-0.5 * u_2*u_2).exp()
                + erf_amp * erfc(u_1) + lin * x + off
        });

        let r = &self.y - f;

        Some(r)
    }

    fn jacobian(&self) -> Option<DMatrix<f64>> {
        let mut j = DMatrix::zeros(self.y.len(), 9);

        let amp_1 = self.p[0];
        let x0_1 = self.p[1];
        let sig_1 = self.p[2];
        let amp_2 = self.p[3];
        let x0_2 = self.p[4];
        let sig_2 = self.p[5];
        let erf_amp = self.p[6];

        for (i, x) in self.x.iter().enumerate() {
            let u_1 = (x - x0_1) / sig_1;
            let z_1 = u_1 / sig_1;
            let h_1 = erf_amp * TWO_OVER_SQRT_PI * (-u_1*u_1).exp() / sig_1;
            let u_2 = (x - x0_2) / sig_2;
            let z_2 = u_2 / sig_2;

            let mut e_1 = (-0.5 * u_1*u_1).exp();
            let mut e_2 = (-0.5 * u_2*u_2).exp();

            j[(i, 0)] = e_1;
            j[(i, 3)] = e_2;

            e_1 *= amp_1;
            e_2 *= amp_2;

            j[(i, 1)] = e_1 * z_1 + h_1;
            j[(i, 2)] = e_1 * u_1 * z_1 + h_1 * u_1;
            j[(i, 4)] = e_2 * z_2;
            j[(i, 5)] = e_2 * u_2 * z_2;
            j[(i, 6)] = erfc(u_1);
            j[(i, 7)] = *x;
            j[(i, 8)] = 1.;
        }

        Some(j)
    }
}


pub fn fit_erf_linear_2_gauss(
    x: &[f64], y: &[f64], init: &[f64]
) -> Result<HashMap<&'static str, SimpleFitParam>, LMFitError> {
    let problem = ErfLinear2GaussProblem {
        x: DVector::from_column_slice(x),
        y: DVector::from_column_slice(y),
        p: DVector::from_column_slice(init),
    };

    let (result, report) = LevenbergMarquardt::new().minimize(problem);

    let fit_status = FitStatus { termination: report.termination };

    if !fit_status.termination.was_successful() {
        return Err(LMFitError::Failure { info: fit_status.repr().into() });
    }

    let mut params = HashMap::new();

    params.insert("amp_1", SimpleFitParam { val: result.p[0], err: f64::NAN });
    params.insert("x0_1", SimpleFitParam { val: result.p[1], err: f64::NAN });
    params.insert("sig_1", SimpleFitParam { val: result.p[2], err: f64::NAN });
    params.insert("amp_2", SimpleFitParam { val: result.p[3], err: f64::NAN });
    params.insert("x0_2", SimpleFitParam { val: result.p[4], err: f64::NAN });
    params.insert("sig_2", SimpleFitParam { val: result.p[5], err: f64::NAN });
    params.insert("erf_amp", SimpleFitParam { val: result.p[6], err: f64::NAN });
    params.insert("lin", SimpleFitParam { val: result.p[7], err: f64::NAN });
    params.insert("const", SimpleFitParam { val: result.p[8], err: f64::NAN });

    Ok(params)
}


struct ErfLinear3GaussProblem {
    x: DVector<f64>,
    y: DVector<f64>,
    p: DVector<f64>,
}

impl LeastSquaresProblem<f64, Dyn, Dyn> for ErfLinear3GaussProblem {
    type ParameterStorage = Owned<f64, Dyn>;
    type ResidualStorage = Owned<f64, Dyn>;
    type JacobianStorage = Owned<f64, Dyn, Dyn>;

    fn set_params(&mut self, params: &DVector<f64>) {
        self.p.copy_from(params);
    }

    fn params(&self) -> DVector<f64> {
        self.p.clone()
    }

    fn residuals(&self) -> Option<DVector<f64>> {
        let amp_1 = self.p[0];
        let x0_1 = self.p[1];
        let sig_1 = self.p[2];
        let amp_2 = self.p[3];
        let x0_2 = self.p[4];
        let sig_2 = self.p[5];
        let amp_3 = self.p[6];
        let x0_3 = self.p[7];
        let sig_3 = self.p[8];
        let erf_amp = self.p[9];
        let lin = self.p[10];
        let off = self.p[11];

        let f = self.x.map(|x| {
            let u_1 = (x - x0_1) / sig_1;
            let u_2 = (x - x0_2) / sig_2;
            let u_3 = (x - x0_3) / sig_3;
            amp_1 * (-0.5 * u_1*u_1).exp()
                + amp_2 * (-0.5 * u_2*u_2).exp()
                + amp_3 * (-0.5 * u_3*u_3).exp()
                + erf_amp * erfc(u_1) + lin * x + off
        });

        let r = &self.y - f;

        Some(r)
    }

    fn jacobian(&self) -> Option<DMatrix<f64>> {
        let mut j = DMatrix::zeros(self.y.len(), 12);

        let amp_1 = self.p[0];
        let x0_1 = self.p[1];
        let sig_1 = self.p[2];
        let amp_2 = self.p[3];
        let x0_2 = self.p[4];
        let sig_2 = self.p[5];
        let amp_3 = self.p[6];
        let x0_3 = self.p[7];
        let sig_3 = self.p[8];
        let erf_amp = self.p[9];

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

            j[(i, 0)] = e_1;
            j[(i, 3)] = e_2;
            j[(i, 6)] = e_3;

            e_1 *= amp_1;
            e_2 *= amp_2;
            e_3 *= amp_3;

            j[(i, 1)] = e_1 * z_1 + h_1;
            j[(i, 2)] = e_1 * u_1 * z_1 + h_1 * u_1;
            j[(i, 4)] = e_2 * z_2;
            j[(i, 5)] = e_2 * u_2 * z_2;
            j[(i, 7)] = e_3 * z_3;
            j[(i, 8)] = e_3 * u_3 * z_3;
            j[(i, 9)] = erfc(u_1);
            j[(i, 10)] = *x;
            j[(i, 11)] = 1.;
        }

        Some(j)
    }
}


pub fn fit_erf_linear_3_gauss(
    x: &[f64], y: &[f64], init: &[f64]
) -> Result<HashMap<&'static str, SimpleFitParam>, LMFitError> {
    let problem = ErfLinear3GaussProblem {
        x: DVector::from_column_slice(x),
        y: DVector::from_column_slice(y),
        p: DVector::from_column_slice(init),
    };

    let (result, report) = LevenbergMarquardt::new().minimize(problem);

    let fit_status = FitStatus { termination: report.termination };

    if !fit_status.termination.was_successful() {
        return Err(LMFitError::Failure { info: fit_status.repr().into() });
    }

    let mut params = HashMap::new();

    params.insert("amp_1", SimpleFitParam { val: result.p[0], err: f64::NAN });
    params.insert("x0_1", SimpleFitParam { val: result.p[1], err: f64::NAN });
    params.insert("sig_1", SimpleFitParam { val: result.p[2], err: f64::NAN });
    params.insert("amp_2", SimpleFitParam { val: result.p[3], err: f64::NAN });
    params.insert("x0_2", SimpleFitParam { val: result.p[4], err: f64::NAN });
    params.insert("sig_2", SimpleFitParam { val: result.p[5], err: f64::NAN });
    params.insert("amp_3", SimpleFitParam { val: result.p[6], err: f64::NAN });
    params.insert("x0_3", SimpleFitParam { val: result.p[7], err: f64::NAN });
    params.insert("sig_3", SimpleFitParam { val: result.p[8], err: f64::NAN });
    params.insert("erf_amp", SimpleFitParam { val: result.p[9], err: f64::NAN });
    params.insert("lin", SimpleFitParam { val: result.p[10], err: f64::NAN });
    params.insert("const", SimpleFitParam { val: result.p[11], err: f64::NAN });

    Ok(params)
}


pub fn erf_linear_background(
    x: &[f64], amp: f64, x0: f64, sig: f64, lin: f64, off: f64
) -> Vec<f64> {
    x.iter().map(|x| {
        let u = (x - x0) / sig;
        amp * erfc(u) + lin * x + off
    }).collect()
}
