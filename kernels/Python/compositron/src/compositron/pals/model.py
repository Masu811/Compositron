from dataclasses import dataclass

import numpy as np

from ..core.fitting import FitParam


FitParamInitializerList = (
    float
    | tuple[float]
    | tuple[float, float]
    | tuple[float, bool]
    | tuple[float, bool, float]
    | tuple[float, float, float]
    | tuple[float, bool, float, float]
)


def parse_parameter_tuple(
    params: FitParamInitializerList,
    default_initialized: bool = False,
) -> FitParam:
    param = FitParam(np.nan)
    param.default_initialized = default_initialized

    if isinstance(params, (float, int)):
        param.val = params
        return param

    param.val = params[0]

    if len(params) == 1:
        return param

    if isinstance(params[1], (float, int)):
        param.min = params[1]

        if len(params) == 2:
            return param

        param.max = params[2]

        return param
    else:
        param.vary = params[1]

        if len(params) == 2:
            return param

        param.min = params[2]

        if len(params) == 3:
            return param

        param.max = params[3]

        return param


@dataclass
class LifetimeComponent:
    lifetime: FitParam
    intensity: FitParam


@dataclass
class ResolutionComponent:
    fwhm: FitParam
    intensity: FitParam
    t0: FitParam


@dataclass
class LifetimeModel:
    lifetime_components: list[LifetimeComponent]
    resolution_components: list[ResolutionComponent]
    scale_component: FitParam
    background_component: FitParam
    shift_component: FitParam


    def __init__(self, model: dict[str, FitParamInitializerList]) -> None:
        model = {key.lower(): val for key, val in model.items()}

        self.lifetime_components = []
        self.resolution_components = []

        n_l = 1

        while True:
            if not f"lifetime_{n_l}" in model:
                break

            lifetime = parse_parameter_tuple(model.pop(f"lifetime_{n_l}"))

            if f"intensity_{n_l}" in model:
                intensity = parse_parameter_tuple(model.pop(f"intensity_{n_l}"))
            else:
                intensity = parse_parameter_tuple(0.5, True)

            self.lifetime_components.append(LifetimeComponent(
                lifetime, intensity
            ))

            n_l += 1

        n_r = 1

        while True:
            if not f"res_fwhm_{n_r}" in model:
                break

            fwhm = parse_parameter_tuple(model.pop(f"res_fwhm_{n_r}"))

            if f"res_intensity_{n_r}" in model:
                intensity = parse_parameter_tuple(model.pop(f"res_intensity_{n_r}"))
            else:
                intensity = parse_parameter_tuple(0.5, True)

            if f"res_t0_{n_r}" in model:
                t0 = parse_parameter_tuple(model.pop(f"res_t0_{n_r}"))
            else:
                t0 = parse_parameter_tuple(0., True)

            self.resolution_components.append(ResolutionComponent(
                fwhm, intensity, t0
            ))

            n_r += 1

        if "background" in model:
            self.background_component = parse_parameter_tuple(model["background"])
        else:
            self.background_component = parse_parameter_tuple((0., False), True)

        if "N" in model:
            self.scale_component = parse_parameter_tuple(model["N"])
        else:
            self.scale_component = parse_parameter_tuple(1e7, True)

        if "t0" in model:
            self.shift_component = parse_parameter_tuple(model["t0"])
        else:
            self.shift_component = parse_parameter_tuple(0., True)
