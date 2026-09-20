use nalgebra::{DMatrix, DVector, Dyn, Const, Owned};
use levenberg_marquardt::{LeastSquaresProblem, LevenbergMarquardt};

use crate::core::fitting::{
    BoundedLeastSquaresProblem, FitParam, FitStatistics, FitStatus, LMFitError,
    SimpleFitParam
};


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


struct Gauss2DProblem {
    x: DVector<f64>,
    y: DVector<f64>,
    z: DVector<f64>,
    p: DVector<f64>,
    params: Vec<FitParam>,
    transformed_params: DVector<f64>,
    transforms: Vec<fn (f64, &FitParam) -> f64>,
    diffs: Vec<fn (f64, &FitParam) -> f64>,
}


impl LeastSquaresProblem<f64, Dyn, Dyn> for Gauss2DProblem {
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
        let y0 = self.transformed_params[2];
        let sig_x = self.transformed_params[3];
        let sig_y = self.transformed_params[4];
        let phi = self.transformed_params[5];

        let f = gauss2d(&self.x, &self.y, amp, x0, y0, sig_x, sig_y, phi);

        let r = &self.z - f;

        Some(r)
    }

    fn jacobian(&self) -> Option<DMatrix<f64>> {
        let mut jac = DMatrix::zeros(self.z.len(), 6);

        let amp = self.transformed_params[0];
        let x0 = self.transformed_params[1];
        let y0 = self.transformed_params[2];
        let sig_x = self.transformed_params[3];
        let sig_y = self.transformed_params[4];
        let phi = self.transformed_params[5];

        let diff_amp = self.diffs[0](self.p[0], &self.params[0]);
        let diff_x0 = self.diffs[1](self.p[1], &self.params[1]);
        let diff_y0 = self.diffs[2](self.p[2], &self.params[2]);
        let diff_sig_x = self.diffs[3](self.p[3], &self.params[3]);
        let diff_sig_y = self.diffs[4](self.p[4], &self.params[4]);
        let diff_phi = self.diffs[5](self.p[5], &self.params[5]);

        let cos_phi = phi.cos();
        let sin_phi = phi.sin();

        let one_over_sig_x = 1. / sig_x;
        let one_over_sig_y = 1. / sig_y;
        let one_over_sig_x_cubed = 1. / sig_x.powi(3);
        let one_over_sig_y_cubed = 1. / sig_y.powi(3);

        let h_x = (cos_phi / sig_x).powi(2) + (sin_phi / sig_y).powi(2);
        let h_y = (sin_phi / sig_x).powi(2) + (cos_phi / sig_y).powi(2);
        let h_phi = (1. / sig_x).powi(2) - (1. / sig_y).powi(2);
        let h = cos_phi * sin_phi * h_phi;

        for j in 0..self.x.len() {
            for i in 0..self.y.len() {
                let flat_idx = j * self.y.len() + i;

                let x_j = self.x[j] - x0;
                let y_i = self.y[i] - y0;

                let x = x_j * cos_phi - y_i * sin_phi;
                let y = x_j * sin_phi + y_i * cos_phi;

                let mut f = (
                    -0.5 * (x * one_over_sig_x).powi(2)
                    - 0.5 * (y * one_over_sig_y).powi(2)
                ).exp();

                jac[(flat_idx, 0)] = f * diff_amp;

                f *= amp;

                jac[(flat_idx, 1)] = f * (x_j * h_x - y_i * h) * diff_x0;
                jac[(flat_idx, 2)] = f * (y_i * h_y - x_j * h) * diff_y0;
                jac[(flat_idx, 3)] = f * x.powi(2) * one_over_sig_x_cubed * diff_sig_x;
                jac[(flat_idx, 4)] = f * y.powi(2) * one_over_sig_y_cubed * diff_sig_y;
                jac[(flat_idx, 5)] = f * x * y * h_phi * diff_phi;
            }
        }

        Some(jac)
    }
}


impl BoundedLeastSquaresProblem for Gauss2DProblem {
    fn diffs(&self) -> &Vec<fn (f64, &FitParam) -> f64> {
        &self.diffs
    }

    fn param_data(&self, i: usize) -> &FitParam {
        &self.params[i]
    }
}


pub struct Gauss2DParams {
    pub amplitude: SimpleFitParam,
    pub x0: SimpleFitParam,
    pub y0: SimpleFitParam,
    pub sig_x: SimpleFitParam,
    pub sig_y: SimpleFitParam,
    pub phi: SimpleFitParam,
}


pub struct Gauss2DFitResult {
    pub params: Gauss2DParams,
    pub stats: FitStatistics,
}


pub fn fit_gauss2d(
    x: &DVector<f64>, y: &DVector<f64>, z: &DMatrix<f64>, init: &[f64]
) -> Result<Gauss2DFitResult, LMFitError> {
    let z = z.to_owned().reshape_generic(Dyn(x.len() * y.len()), Const::<1>);

    let mut params = Vec::new();

    params.push(FitParam::new(init[0]).min(0.));  // amp
    params.push(FitParam::new(init[1]));  // x0
    params.push(FitParam::new(init[2]));  // y0
    params.push(FitParam::new(init[3]).min(0.1));  // sig_x
    params.push(FitParam::new(init[4]).min(0.1));  // sig_y
    params.push(FitParam::new(init[5]));  // phi

    let transforms = params.iter().map(|param| param.transform).collect();
    let backtransforms = params.iter().map(|param| param.backtransform).collect::<Vec<_>>();
    let diffs = params.iter().map(|param| param.diff).collect();

    let mut p = Vec::new();

    for i in 0..init.len() {
        p.push(backtransforms[i](init[i], &params[i]));
    }

    let problem = Gauss2DProblem {
        x: x.to_owned(),
        y: y.to_owned(),
        z: z,
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

    let params = Gauss2DParams {
        amplitude: SimpleFitParam { val: vals[0], err: errs[0] },
        x0: SimpleFitParam { val: vals[1], err: errs[1] },
        y0: SimpleFitParam { val: vals[2], err: errs[2] },
        sig_x: SimpleFitParam { val: vals[3], err: errs[3] },
        sig_y: SimpleFitParam { val: vals[4], err: errs[4] },
        phi: SimpleFitParam { val: vals[5], err: errs[5] },
    };

    Ok(Gauss2DFitResult { params, stats })
}
