use levenberg_marquardt::TerminationReason;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum VarproFitError {
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
pub struct SimpleFitParam {
    pub val: f64,
    pub err: f64,
}

fn min_transform(value: f64, param: &FitParam) -> f64 {
    param.min - 1. + (value.powi(2) + 1.).sqrt()
}

fn min_backtransform(value: f64, param: &FitParam) -> f64 {
    ((value - param.min + 1.).powi(2) - 1.).sqrt()
}

fn min_diff(value: f64, _: &FitParam) -> f64 {
    value / (value.powi(2) + 1.).sqrt()
}

fn max_transform(value: f64, param: &FitParam) -> f64 {
    param.max + 1. - (value.powi(2) + 1.).sqrt()
}

fn max_backtransform(value: f64, param: &FitParam) -> f64 {
    ((param.max - value + 1.).powi(2) - 1.).sqrt()
}

fn max_diff(value: f64, _: &FitParam) -> f64 {
    -value / (value.powi(2) + 1.).sqrt()
}

fn minmax_transform(value: f64, param: &FitParam) -> f64 {
    param.min + (param.max - param.min) / 2. * (1. + value.sin())
}

fn minmax_backtransform(value: f64, param: &FitParam) -> f64 {
    (2. * (value - param.min) / (param.max - param.min) - 1.).asin()
}

fn minmax_diff(value: f64, param: &FitParam) -> f64 {
    (param.max - param.min) / 2. * value.cos()
}

fn identity_transform(value: f64, _: &FitParam) -> f64 {
    value
}

fn identity_diff(_: f64, _: &FitParam) -> f64 {
    1.
}

#[derive(Clone, Copy)]
pub struct FitParam {
    pub val: f64,
    pub err: f64,
    pub min: f64,
    pub max: f64,
    pub vary: bool,
    pub transform: fn (f64, &FitParam) -> f64,
    pub backtransform: fn (f64, &FitParam) -> f64,
    pub diff: fn (f64, &FitParam) -> f64,
    pub default_initialized: bool,
}

impl std::fmt::Debug for FitParam {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "FitParam(val={}, err={}, min={}, max={}, vary={})",
            self.val,
            self.err,
            self.min,
            self.max,
            self.vary,
        )
    }
}

impl FitParam {
    pub fn new<T: Into<f64>>(value: T) -> Self {
        FitParam {
            val: value.into(),
            err: f64::NAN,
            min: f64::NEG_INFINITY,
            max: f64::INFINITY,
            vary: true,
            transform: identity_transform,
            backtransform: identity_transform,
            diff: identity_diff,
            default_initialized: false,
        }
    }

    pub fn default<T: Into<f64>>(value: T) -> Self {
        let mut p = FitParam::new(value);
        p.default_initialized = true;
        p
    }

    pub fn min<T: Into<f64>>(mut self, min: T) -> Self {
        self.min = min.into();

        if self.max.is_infinite() {
            self.transform = min_transform;
            self.backtransform = min_backtransform;
            self.diff = min_diff;
        } else {
            self.transform = minmax_transform;
            self.backtransform = minmax_backtransform;
            self.diff = minmax_diff;
        }

        self
    }

    pub fn max<T: Into<f64>>(mut self, max: T) -> Self {
        self.max = max.into();

        if self.min.is_infinite() {
            self.transform = max_transform;
            self.backtransform = max_backtransform;
            self.diff = max_diff;
        } else {
            self.transform = minmax_transform;
            self.backtransform = minmax_backtransform;
            self.diff = minmax_diff;
        }

        self
    }

    pub fn fixed(mut self) -> Self {
        self.vary = false;
        self
    }

    pub fn transform(&self, value: f64) -> f64 {
        (self.transform)(value, self)
    }

    pub fn backtransform(&self, value: f64) -> f64 {
        (self.backtransform)(value, self)
    }

    pub fn diff(&self, value: f64) -> f64 {
        (self.diff)(value, self)
    }
}

#[derive(Debug, Error)]
pub enum LMFitError {
    #[error("Fit failed: {info}")]
    Failure {
        info: String,
    },
}

pub struct FitStatus {
    pub termination: TerminationReason,
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
