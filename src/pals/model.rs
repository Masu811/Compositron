use std::collections::HashMap;

use thiserror::Error;

use crate::core::fitting::FitParam;

#[derive(Debug)]
pub struct LifetimeComponent {
    pub lifetime: FitParam,
    pub intensity: FitParam,
}

#[derive(Debug)]
pub struct ResolutionComponent {
    pub fwhm: FitParam,
    pub intensity: FitParam,
    pub t0: FitParam,
}

#[derive(Debug)]
pub struct LifetimeModel {
    pub lifetime_components: Vec<LifetimeComponent>,
    pub resolution_components: Vec<ResolutionComponent>,
    pub background_component: FitParam,
    pub scale_component: FitParam,
    pub shift_component: FitParam,
}

#[derive(Debug, Error)]
pub enum ParseModelError {
    #[error("Lifetime Model is missing {parameter}")]
    MissingParameter {
        parameter: String,
    }
}

pub fn parse_model(
    mut model: HashMap<String, FitParam>
) -> Result<LifetimeModel, ParseModelError> {
    let mut lifetime_components = Vec::new();
    let mut n_l = 1;

    while let Some(lifetime) = model.remove(&format!("lifetime_{n_l}")) {
        let i_key = format!("intensity_{n_l}");
        lifetime_components.push(LifetimeComponent {
            lifetime,
            intensity: model.remove(&i_key)
                .unwrap_or(FitParam::default(0.5))
        });
        n_l += 1;
    }

    let mut resolution_components = Vec::new();
    let mut n_r = 1;

    while let Some(fwhm) = model.remove(&format!("res_fwhm_{n_r}")) {
        let i_key = format!("res_intensity_{n_r}");
        let t_key = format!("res_t0_{n_r}");
        resolution_components.push(ResolutionComponent {
            fwhm,
            intensity: model.remove(&i_key)
                .unwrap_or(FitParam::default(0.5)),
            t0: model.remove(&t_key)
                .unwrap_or(FitParam::default(0.)),
        });
        n_r += 1;
    }

    if n_l == 1 || n_r == 1 {
        if n_l == 1 && n_r == 1 {
            return Err(ParseModelError::MissingParameter {
                parameter: "lifetime and resolution".into(),
            });
        } else if n_l == 1 {
            return Err(ParseModelError::MissingParameter {
                parameter: "lifetime".into(),
            });
        } else {
            return Err(ParseModelError::MissingParameter {
                parameter: "resolution".into(),
            });
        }
    }

    Ok(LifetimeModel {
        lifetime_components: lifetime_components,
        resolution_components: resolution_components,
        background_component: model.remove("background")
            .unwrap_or(FitParam::default(0.).fixed()),
        scale_component: model.remove("N")
            .unwrap_or(FitParam::default(1e7).min(0.)),
        shift_component: model.remove("t0")
            .unwrap_or(FitParam::default(0.)),
    })
}

#[macro_export]
macro_rules! parse_parameter_tuple {
    ($value:literal) => {
        $crate::core::fitting::FitParam::new($value)
    };
    (($value:expr)) => {
        $crate::core::fitting::FitParam::new($value)
    };
    (($value:expr, true)) => {
        $crate::core::fitting::FitParam::new($value)
    };
    (($value:expr, false)) => {
        $crate::core::fitting::FitParam::new($value).fixed()
    };
    (($value:expr, $min:expr)) => {
        $crate::core::fitting::FitParam::new($value).min($min)
    };
    (($value:expr, true, $min:expr)) => {
        $crate::core::fitting::FitParam::new($value).min($min)
    };
    (($value:expr, false, $min:expr)) => {
        $crate::core::fitting::FitParam::new($value).min($min).fixed()
    };
    (($value:expr, $min:expr, $max:expr)) => {
        $crate::core::fitting::FitParam::new($value).min($min).max($max)
    };
    (($value:expr, true, $min:expr, $max:expr)) => {
        $crate::core::fitting::FitParam::new($value).min($min).max($max)
    };
    (($value:expr, false, $min:expr, $max:expr)) => {
        $crate::core::fitting::FitParam::new($value).min($min).max($max).fixed()
    };
}

#[macro_export]
macro_rules! lifetime_model {
    ($($name:literal: $tuple:tt),*) => {
        $crate::pals::model::parse_model(
            std::collections::HashMap::<String, $crate::core::fitting::FitParam>::from([
                $(($name.into(), $crate::parse_parameter_tuple!($tuple))),*
            ])
        )
    };
}
