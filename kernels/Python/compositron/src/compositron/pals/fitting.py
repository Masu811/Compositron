from copy import deepcopy
from dataclasses import dataclass
from itertools import product
from typing import Callable

import lmfit
from lmfit import Parameter
import numpy as np
import numpy.typing as npt
from scipy.special import erfc

from compositron.core.fitting import FitParam

from ..constants import FWHM_OVER_SIGMA, SQRT_2, SQRT_2_PI
from .model import LifetimeComponent, LifetimeModel, ResolutionComponent


def lifetime_spectrum_component(
    x: npt.NDArray[np.floating],
    n: float,
    t0: float,
    tau: float,
    l_int: float,
    sig: float,
    r_int: float,
    r_t0: float,
) -> npt.NDArray[np.floating]:
    sig_over_tau_sqrt_2 = sig / (tau * SQRT_2)
    sig_over_tau_sqrt_2_sq = sig_over_tau_sqrt_2 ** 2
    one_over_tau = 1 / tau
    one_over_sig_sqrt_2 = 1 / (sig * SQRT_2)

    t_shifted = x - (t0 + r_t0)
    e = np.exp(-t_shifted * one_over_tau + sig_over_tau_sqrt_2_sq)
    c = erfc(sig_over_tau_sqrt_2 - t_shifted * one_over_sig_sqrt_2)

    return 0.5 * n * l_int * r_int * one_over_tau * e * c


@dataclass
class LifetimeFit:
    n_l: int
    n_r: int
    model: LifetimeModel
    l_vary: list[bool]
    r_vary: list[bool]
    l_last_vary_idx: int
    r_last_vary_idx: int
    available_l_intensity: float
    available_r_intensity: float

    def __init__(self, model: LifetimeModel) -> None:
        self.model = model

        self.n_l = len(self.model.lifetime_components)
        self.n_r = len(self.model.resolution_components)

        self.l_vary = [
            lt.intensity.vary for lt in self.model.lifetime_components
        ]
        l_vary_idcs = [i for i, vary in enumerate(self.l_vary) if vary]
        self.l_last_vary_idx = (
            max(l_vary_idcs) if len(l_vary_idcs) > 0 else 0
        )

        self.r_vary = [
            res.intensity.vary for res in self.model.resolution_components
        ]
        r_vary_idcs = [j for j, vary in enumerate(self.r_vary) if vary]
        self.r_last_vary_idx = (
            max(r_vary_idcs) if len(r_vary_idcs) > 0 else 0
        )

        self.available_l_intensity = 1 - sum(
            lt_comp.intensity.val
            for lt_comp in model.lifetime_components if not lt_comp.intensity.vary
        )
        self.available_r_intensity = 1 - sum(
            res_comp.intensity.val
            for res_comp in model.resolution_components if not res_comp.intensity.vary
        )


    def make_model_func(self) -> tuple[Callable, Callable]:
        def func(x, **params):
            n = params["N"]
            t0 = params["t0"]
            bg = params["background"]
            y = np.full_like(x, bg, dtype=float)

            l_int = [
                params[f"h_{i}"] if v else params[f"intensity_{i}"]
                for i, v in enumerate(self.l_vary, 1)
            ]
            l_cumul = sum(i for i, v in zip(l_int, self.l_vary) if not v)
            for i in range(self.n_l):
                if self.l_vary[i]:
                    l_int[i] *= max(1 - l_cumul, 0)
                    l_cumul += l_int[i]

            r_int = [
                params[f"res_h_{j}"] if v else params[f"res_intensity_{j}"]
                for j, v in enumerate(self.r_vary, 1)
            ]
            r_cumul = sum(i for i, v in zip(r_int, self.r_vary) if not v)
            for j in range(self.n_r):
                if self.r_vary[j]:
                    r_int[j] *= max(1 - r_cumul, 0)
                    r_cumul += r_int[j]

            for i, j in product(range(1, self.n_l+1), range(1, self.n_r+1)):
                tau = params[f"lifetime_{i}"]
                l_i = l_int[i-1]
                sig = params[f"res_sigma_{j}"]
                r_i = r_int[j-1]
                r_t0 = params[f"res_t0_{j}"]

                y += lifetime_spectrum_component(
                    x, n, t0, tau, l_i, sig, r_i, r_t0,
                )

            return y

        def dfunc(params, _, weights, x):
            n = params["N"]
            t0 = params["t0"]

            l_int = [
                params[f"h_{i}"] if v else params[f"intensity_{i}"]
                for i, v in enumerate(self.l_vary, 1)
            ]
            l_cumul = sum(i for i, v in zip(l_int, self.l_vary) if not v)
            for i in range(self.n_l):
                if self.l_vary[i]:
                    l_int[i] *= max(1 - l_cumul, 0)
                    l_cumul += l_int[i]

            r_int = [
                params[f"res_h_{j}"] if v else params[f"res_intensity_{j}"]
                for j, v in enumerate(self.r_vary, 1)
            ]
            r_cumul = sum(i for i, v in zip(r_int, self.r_vary) if not v)
            for j in range(self.n_r):
                if self.r_vary[j]:
                    r_int[j] *= max(1 - r_cumul, 0)
                    r_cumul += r_int[j]

            terms = {}
            deriv = {}

            y = np.zeros_like(x)

            for i, j in product(range(1, self.n_l+1), range(1, self.n_r+1)):
                l_tau = params[f"lifetime_{i}"]
                l_i = l_int[i-1]
                r_sigma = params[f"res_sigma_{j}"]
                r_i = r_int[j-1]
                r_t0 = params[f"res_t0_{j}"]
                # e becomes numerically unstable for lifetimes ~ tcal[1]
                _t = x - (t0 + r_t0)
                a = -_t / l_tau + (r_sigma/(SQRT_2*l_tau))**2
                b = r_sigma / (l_tau * SQRT_2) - _t / (r_sigma * SQRT_2)
                e = np.exp(a)
                c = erfc(b)
                decay = n * l_i * r_i / (2*l_tau) * e * c
                y += decay
                exp_comb = n * (l_i * r_i) / (SQRT_2_PI * l_tau) * np.exp(a - b**2)

                df_dt0 = (
                    -exp_comb / r_sigma
                    + decay / l_tau
                )
                df_dtau = (
                    exp_comb * r_sigma / l_tau**2
                    + decay * (-r_sigma**2/l_tau**3 + _t/l_tau**2 - 1/l_tau)
                )
                df_dsigma = (
                    decay * r_sigma / l_tau**2
                    - exp_comb * (1/l_tau + _t/r_sigma**2)
                )
                df_dr_t0 = df_dt0

                terms[f"{i}_{j}"] = decay

                deriv[f"{i}_{j}"] = [
                    df_dt0,
                    df_dtau,
                    df_dsigma,
                    df_dr_t0,
                ]

            df_dN = y / n
            df_dbackground = np.ones_like(x, dtype=float)
            df_dt0 = np.sum([
                deriv[f"{i}_{j}"][0]
                for j in range(1, self.n_r+1)
                for i in range(1, self.n_l+1)
            ], axis=0)

            l = self.l_vary.copy()
            l[self.l_last_vary_idx] = False
            r = self.r_vary.copy()
            r[self.r_last_vary_idx] = False

            df_dh_i = []

            for k in range(1, self.n_l+1):
                df_dh_k = []
                if not l[k-1]:
                    df_dh_i.append(np.zeros_like(x, dtype=float))
                    continue
                for n, m in product(range(1, self.n_l+1), range(1, self.n_r+1)):
                    if n < k:
                        continue
                    elif n == k:
                        df_dh_k.append(
                            terms[f"{n}_{m}"] / (params[f"h_{k}"] + 1e-5)
                        )
                    else:
                        df_dh_k.append(
                            -terms[f"{n}_{m}"] / (1 - params[f"h_{k}"] + 1e-5)
                        )

                df_dh_i.append(np.sum(df_dh_k, axis=0))

            df_dres_h_j = []

            for k in range(1, self.n_r+1):
                df_dres_h_k = []
                if not r[k-1]:
                    df_dres_h_j.append(np.zeros_like(x, dtype=float))
                    continue
                for n, m in product(range(1, self.n_l+1), range(1, self.n_r+1)):
                    if m < k:
                        continue
                    elif m == k:
                        df_dres_h_k.append(
                            terms[f"{n}_{m}"] / (params[f"res_h_{k}"] + 1e-5)
                        )
                    else:
                        df_dres_h_k.append(
                            -terms[f"{n}_{m}"] / (1 - params[f"res_h_{k}"] + 1e-4)
                        )

                df_dres_h_j.append(np.sum(df_dres_h_k, axis=0))

            out = {
                "N": df_dN,
                "t0": df_dt0,
                "background": df_dbackground,
                **{
                    f"lifetime_{i}": np.sum(
                        [deriv[f"{i}_{j}"][1] for j in range(1, self.n_r+1)], axis=0
                    ) for i in range(1, self.n_l+1)
                },
                **{
                    f"h_{i}": df_dh_i[i-1] for i in range(1, self.n_l+1)
                },
                **{
                    f"res_sigma_{j}": np.sum(
                        [deriv[f"{i}_{j}"][2] for i in range(1, self.n_l+1)], axis=0
                    ) for j in range(1, self.n_r+1)
                },
                **{
                    f"res_h_{j}": df_dres_h_j[j-1] for j in range(1, self.n_r+1)
                },
                **{
                    f"res_t0_{j}": np.sum(
                        [deriv[f"{i}_{j}"][3] for i in range(1, self.n_l+1)], axis=0
                    ) for j in range(1, self.n_r+1)
                }
            }

            out = np.array([out[name] for name, param in params.items() if param.vary])

            return -out * weights

        return func, dfunc


    def make_model(
        self,
        x: npt.NDArray[np.floating],
        counts: int,
        peak_center: float,
    ) -> tuple[lmfit.Model, lmfit.Parameters, Callable, Callable]:
        func, dfunc = self.make_model_func()

        model = lmfit.Model(func)
        params = lmfit.Parameters()

        #### General Parameters

        if self.model.scale_component.default_initialized:
            params.add(Parameter("N", 2 * counts, min=0))
        else:
            p = self.model.scale_component
            params.add(Parameter("N", p.val, p.vary, p.min, p.max))

        if self.model.shift_component.default_initialized:
            l = float(np.min(x))
            r = float(np.max(x))
            d = r - l
            params.add(Parameter("t0", peak_center, min=l-d, max=r+d))
        else:
            p = self.model.shift_component
            params.add(Parameter("t0", p.val, p.vary, p.min, p.max))

        if self.model.background_component.default_initialized:
            params.add(Parameter("background", 0, vary=False))
        else:
            p = self.model.background_component
            params.add(Parameter("background", p.val, p.vary, p.min, p.max))

        #### Lifetime Components

        remaining_l_int = self.available_l_intensity

        for i, l_comp in enumerate(self.model.lifetime_components, 1):
            lt = l_comp.lifetime
            params.add(Parameter(
                f"lifetime_{i}", lt.val, lt.vary, lt.min, lt.max
            ))

            l_int = l_comp.intensity
            params.add(Parameter(f"intensity_{i}", l_int.val, False))

            if i == self.l_last_vary_idx + 1:
                params.add(Parameter(f"h_{i}", 1., False))
            else:
                if l_int.vary:
                    q = l_int.val / remaining_l_int
                    remaining_l_int -= l_int.val
                else:
                    q = 1

                params.add(Parameter(
                    f"h_{i}",
                    value = q,
                    vary = l_int.vary,
                    min = min(max(l_int.min, 0), 1),
                    max = max(min(l_int.max, 1), 0),
                ))

        #### Resolution Components

        remaining_r_int = self.available_r_intensity

        for i, r_comp in enumerate(self.model.resolution_components, 1):
            fwhm = r_comp.fwhm
            params.add(Parameter(
                f"res_sigma_{i}",
                fwhm.val / FWHM_OVER_SIGMA,
                fwhm.vary,
                fwhm.min / FWHM_OVER_SIGMA,
                fwhm.max / FWHM_OVER_SIGMA,
            ))

            r_int = r_comp.intensity
            params.add(Parameter(f"res_intensity_{i}", r_int.val, False))

            if i == self.r_last_vary_idx + 1:
                params.add(Parameter(f"res_h_{i}", 1., False))
            else:
                if r_int.vary:
                    q = r_int.val / remaining_r_int
                    remaining_r_int -= r_int.val
                else:
                    q = 1

                params.add(Parameter(
                    f"res_h_{i}",
                    value = q,
                    vary = r_int.vary,
                    min = min(max(r_int.min, 0), 1),
                    max = max(min(r_int.max, 1), 0),
                ))

            t0 = r_comp.t0
            params.add(Parameter(
                f"res_t0_{i}", t0.val, t0.vary, t0.min, t0.max
            ))

        if all(r_comp.t0.vary for r_comp in self.model.resolution_components):
            params["res_t0_1"].vary = False

        return model, params, func, dfunc


    def post_fit(self, fit_result: lmfit.model.ModelResult) -> FitResult:
        def may_be_nan(x):
            return x if x is not None else np.nan

        params = fit_result.params

        l_int = [
            params[f"h_{i}"].value if v else params[f"intensity_{i}"].value
            for i, v in enumerate(self.l_vary, 1)
        ]
        h = l_int.copy()
        l_dint = [
            may_be_nan(params[f"h_{i}"].stderr)
            for i in range(1, self.n_l+1)
        ]
        dh = l_dint.copy()
        l_cumul = sum(i for i, v in zip(l_int, self.l_vary) if not v)
        l_cumul_sq_err = 0
        for i in range(self.n_l):
            if self.l_vary[i]:
                l_int[i] *= max(1 - l_cumul, 0)
                l_cumul += l_int[i]
                l_dint[i] = l_int[i] * np.sqrt(
                    l_cumul_sq_err + (dh[i] / (h[i] + 1e-12))**2
                )
                if i == self.l_last_vary_idx:
                    break
                l_cumul_sq_err += (dh[i] / (1 - h[i] + 1e-12))**2

        for i in range(1, self.n_l+1):
            params[f"intensity_{i}"] = param = params.pop(f"h_{i}")
            param.name = f"intensity_{i}"
            param.value = l_int[i-1]
            param.vary = self.l_vary[i-1]
            param.stderr = l_dint[i-1]

        r_int = [
            params[f"res_h_{j}"].value if v else params[f"res_intensity_{j}"].value
            for j, v in enumerate(self.r_vary, 1)
        ]
        h = r_int.copy()
        r_dint = [
            may_be_nan(params[f"res_h_{j}"].stderr)
            for j in range(1, self.n_r+1)
        ]
        dh = r_dint.copy()
        r_cumul = sum(i for i, v in zip(r_int, self.r_vary) if not v)
        r_cumul_sq_err = 0
        for j in range(self.n_r):
            if self.r_vary[j]:
                r_int[j] *= max(1 - r_cumul, 0)
                r_cumul += r_int[j]
                r_dint[j] = r_int[j] * np.sqrt(
                    r_cumul_sq_err + (dh[j] / (h[j] + 1e-12))**2
                )
                if j == self.r_last_vary_idx:
                    break
                r_cumul_sq_err += (dh[j] / (1 - h[j] + 1e-12))**2

        for j in range(1, self.n_r+1):
            params[f"res_intensity_{j}"] = param = params.pop(f"res_h_{j}")
            param.name = f"res_intensity_{j}"
            param.value = r_int[j-1]
            param.vary = self.r_vary[j-1]
            param.stderr = r_dint[j-1]

        lt_model = LifetimeModel({})

        n = params["N"]
        lt_model.scale_component = FitParam(
            n.value, may_be_nan(n.stderr), n.min, n.max, n.vary
        )

        bg = params["background"]
        lt_model.background_component = FitParam(
            bg.value, may_be_nan(bg.stderr), bg.min, bg.max, bg.vary
        )

        t0 = params["t0"]
        lt_model.shift_component = FitParam(
            t0.value, may_be_nan(t0.stderr), t0.min, t0.max, t0.vary
        )

        for i in range(self.n_l):
            lt = params[f"lifetime_{i+1}"]
            l_int = params[f"intensity_{i+1}"]

            lt_model.lifetime_components.append(LifetimeComponent(
                lifetime = FitParam(
                    lt.value, may_be_nan(lt.stderr), lt.min, lt.max, lt.vary
                ),
                intensity = FitParam(
                    l_int.value,
                    may_be_nan(l_int.stderr),
                    l_int.min,
                    l_int.max,
                    self.l_vary[i],
                )
            ))

        for i in range(self.n_r):
            fwhm = params[f"res_sigma_{i+1}"]
            r_int = params[f"res_intensity_{i+1}"]
            t0 = params[f"res_t0_{i+1}"]

            lt_model.resolution_components.append(ResolutionComponent(
                fwhm = FitParam(
                    fwhm.value * FWHM_OVER_SIGMA,
                    may_be_nan(fwhm.stderr) * FWHM_OVER_SIGMA,
                    fwhm.min * FWHM_OVER_SIGMA,
                    fwhm.max * FWHM_OVER_SIGMA,
                    fwhm.vary,
                ),
                intensity = FitParam(
                    r_int.value,
                    may_be_nan(r_int.stderr),
                    r_int.min,
                    r_int.max,
                    self.r_vary[i],
                ),
                t0 = FitParam(
                    t0.value, may_be_nan(t0.stderr), t0.min, t0.max, t0.vary
                )
            ))

        return FitResult(
            lt_model,
            fit_result.ndata,
            fit_result.success,
            fit_result.nfev,
            fit_result.redchi,
            fit_result.covar
        )


    def fit(
        self,
        x: npt.NDArray[np.floating],
        y: npt.NDArray[np.floating],
        counts: int,
        peak_center: float,
    ) -> FitResult:
        fit_model, params, _, jac = self.make_model(x, counts, peak_center)

        w = 1 / np.sqrt(np.maximum(y, 1))

        fit_result = fit_model.fit(
            y, x=x, params=params, weights=w,
            fit_kws={"Dfun": jac, "col_deriv": True}
        )

        return self.post_fit(fit_result)


@dataclass
class FitResult:
    model: LifetimeModel
    n_dpoints: int
    fit_status: int
    n_eval: int
    red_chi_2: None | float
    cov: None | npt.NDArray[np.floating]


def fit_lifetime_spectrum(
    x: npt.NDArray[np.floating],
    y: npt.NDArray[np.floating],
    model: LifetimeModel,
    counts: int,
    peak_center: float,
) -> FitResult:
    return LifetimeFit(model).fit(x, y, counts, peak_center)
