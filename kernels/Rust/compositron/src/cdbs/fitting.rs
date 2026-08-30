use std::collections::HashMap;

use nalgebra::{DMatrix, DVector, Dyn, Const, Owned};
use levenberg_marquardt::{LeastSquaresProblem, LevenbergMarquardt};

use crate::core::fitting::{FitStatus, LMFitError, SimpleFitParam};

fn gauss2d(
    x: &DVector<f64>,
    y: &DVector<f64>,
    amp: f64,
    x0: f64,
    y0: f64,
    sig_x: f64,
    sig_y: f64,
    phi: f64,
) -> DVector<f64> {
    let cos_phi = phi.cos();
    let sin_phi = phi.sin();

    let mat = DMatrix::from_fn(y.len(), x.len(), |i, j| {
        let x_j = (x[j] - x0) * cos_phi - (y[i] - y0) * sin_phi;
        let y_i = (x[j] - x0) * sin_phi + (y[i] - y0) * cos_phi;
        amp * (-0.5 * (x_j / sig_x).powi(2) - 0.5 * (y_i / sig_y).powi(2)).exp()
    });

    mat.reshape_generic(Dyn(y.len() * x.len()), Const::<1>)
}

fn gauss2d_damp(
    _: &DVector<f64>,
    _: &DVector<f64>,
    f: &DVector<f64>,
    amp: f64,
    _: f64,
    _: f64,
    _: f64,
    _: f64,
    _: f64,
) -> DVector<f64> {
    f / amp
}

fn gauss2d_dx0(
    x: &DVector<f64>,
    y: &DVector<f64>,
    f: &DVector<f64>,
    _: f64,
    x0: f64,
    y0: f64,
    sig_x: f64,
    sig_y: f64,
    phi: f64,
) -> DVector<f64> {
    let cos_phi = phi.cos();
    let sin_phi = phi.sin();

    let mat = DMatrix::from_fn(y.len(), x.len(), |i, j| {
        (x[j] - x0) * ((cos_phi / sig_x).powi(2) + (sin_phi / sig_y).powi(2))
         - (y[i] - y0) * cos_phi * sin_phi
            * ((1. / sig_x).powi(2) - (1. / sig_y).powi(2))
    });

    f.component_mul(&mat.reshape_generic(Dyn(y.len() * x.len()), Const::<1>))
}

fn gauss2d_dy0(
    x: &DVector<f64>,
    y: &DVector<f64>,
    f: &DVector<f64>,
    _: f64,
    x0: f64,
    y0: f64,
    sig_x: f64,
    sig_y: f64,
    phi: f64,
) -> DVector<f64> {
    let cos_phi = phi.cos();
    let sin_phi = phi.sin();

    let mat = DMatrix::from_fn(y.len(), x.len(), |i, j| {
        (y[i] - y0) * ((sin_phi / sig_x).powi(2) + (cos_phi / sig_y).powi(2))
         - (x[j] - x0) * cos_phi * sin_phi
            * ((1. / sig_x).powi(2) - (1. / sig_y).powi(2))
    });

    f.component_mul(&mat.reshape_generic(Dyn(y.len() * x.len()), Const::<1>))
}

fn gauss2d_dsig_x(
    x: &DVector<f64>,
    y: &DVector<f64>,
    f: &DVector<f64>,
    _: f64,
    x0: f64,
    y0: f64,
    sig_x: f64,
    _: f64,
    phi: f64,
) -> DVector<f64> {
    let cos_phi = phi.cos();
    let sin_phi = phi.sin();

    let mat = DMatrix::from_fn(y.len(), x.len(), |i, j| {
        let x_j = (x[j] - x0) * cos_phi - (y[i] - y0) * sin_phi;
        x_j.powi(2) / sig_x.powi(3)
    });

    f.component_mul(&mat.reshape_generic(Dyn(y.len() * x.len()), Const::<1>))
}

fn gauss2d_dsig_y(
    x: &DVector<f64>,
    y: &DVector<f64>,
    f: &DVector<f64>,
    _: f64,
    x0: f64,
    y0: f64,
    _: f64,
    sig_y: f64,
    phi: f64,
) -> DVector<f64> {
    let cos_phi = phi.cos();
    let sin_phi = phi.sin();

    let mat = DMatrix::from_fn(y.len(), x.len(), |i, j| {
        let y_i = (x[j] - x0) * sin_phi + (y[i] - y0) * cos_phi;
        y_i.powi(2) / sig_y.powi(3)
    });

    f.component_mul(&mat.reshape_generic(Dyn(y.len() * x.len()), Const::<1>))
}

fn gauss2d_dphi(
    x: &DVector<f64>,
    y: &DVector<f64>,
    f: &DVector<f64>,
    _: f64,
    x0: f64,
    y0: f64,
    sig_x: f64,
    sig_y: f64,
    phi: f64,
) -> DVector<f64> {
    let cos_phi = phi.cos();
    let sin_phi = phi.sin();

    let mat = DMatrix::from_fn(y.len(), x.len(), |i, j| {
        let x_j = (x[j] - x0) * cos_phi - (y[i] - y0) * sin_phi;
        let y_i = (x[j] - x0) * sin_phi + (y[i] - y0) * cos_phi;
        x_j * y_i * (1. / sig_x.powi(2) - 1. / sig_y.powi(2))
    });

    f.component_mul(&mat.reshape_generic(Dyn(y.len() * x.len()), Const::<1>))
}

struct Problem {
    x: DVector<f64>,
    y: DVector<f64>,
    z: DVector<f64>,
    p: DVector<f64>,
}

impl LeastSquaresProblem<f64, Dyn, Dyn> for Problem {
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
        let f = gauss2d(
            &self.x,
            &self.y,
            self.p[0],
            self.p[1],
            self.p[2],
            self.p[3],
            self.p[4],
            self.p[5],
        );

        let r = &self.z - f;

        Some(r)
    }

    fn jacobian(&self) -> Option<DMatrix<f64>> {
        let mut j = DMatrix::zeros(self.z.len(), 6);

        let f = gauss2d(
            &self.x,
            &self.y,
            self.p[0],
            self.p[1],
            self.p[2],
            self.p[3],
            self.p[4],
            self.p[5],
        );

        for (i, func) in [
            gauss2d_damp,
            gauss2d_dx0,
            gauss2d_dy0,
            gauss2d_dsig_x,
            gauss2d_dsig_y,
            gauss2d_dphi,
        ].iter().enumerate() {
            j.set_column(i, &(-func(
                &self.x,
                &self.y,
                &f,
                self.p[0],
                self.p[1],
                self.p[2],
                self.p[3],
                self.p[4],
                self.p[5],
            )));
        }

        Some(j)
    }
}

pub fn fit_gauss2d(
    x: &DVector<f64>, y: &DVector<f64>, z: &DMatrix<f64>, init: &[f64]
) -> Result<HashMap<&'static str, SimpleFitParam>, LMFitError> {
    let z = z.to_owned().reshape_generic(Dyn(x.len() * y.len()), Const::<1>);

    let problem = Problem {
        x: x.to_owned(),
        y: y.to_owned(),
        z: z,
        p: DVector::from_column_slice(init),
    };

    let (result, report) = LevenbergMarquardt::new().minimize(problem);

    let fit_status = FitStatus { termination: report.termination };

    if !fit_status.termination.was_successful() {
        return Err(LMFitError::Failure { info: fit_status.repr().into() });
    }

    let mut params = HashMap::new();

    params.insert("amp", SimpleFitParam { val: result.p[0], err: f64::NAN });
    params.insert("x0", SimpleFitParam { val: result.p[1], err: f64::NAN });
    params.insert("y0", SimpleFitParam { val: result.p[2], err: f64::NAN });
    params.insert("sig_x", SimpleFitParam { val: result.p[3], err: f64::NAN });
    params.insert("sig_y", SimpleFitParam { val: result.p[4], err: f64::NAN });
    params.insert("phi", SimpleFitParam { val: result.p[5], err: f64::NAN });

    Ok(params)
}
