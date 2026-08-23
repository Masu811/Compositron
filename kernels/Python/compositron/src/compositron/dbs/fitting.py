import numpy as np
import numpy.typing as npt
from scipy.special import erf
from lmfit import Model

from ..constants import TWO_OVER_SQRT_PI
from ..core.fitting import SimpleFitParam


def gauss(
    x: npt.NDArray[np.floating],
    amp: float,
    x0: float,
    sig: float,
) -> npt.NDArray[np.floating]:
    return amp * np.exp(-0.5 * ((x - x0) / sig)**2)


def fit_gaussian(
    x: npt.NDArray[np.floating],
    y: npt.NDArray[np.floating],
    init: list[float | np.floating],
) -> dict[str, SimpleFitParam]:
    model = Model(gauss)

    def jac(params, _, weights, x):
        amp = params["amp"]
        x0 = params["x0"]
        sig = params["sig"]

        j = np.empty((3, x.size), dtype=np.float64)

        e = gauss(x, amp, x0, sig)
        u = (x - x0) / sig
        z = u / sig

        j[0, :] = e / amp
        j[1, :] = e * z
        j[2, :] = e * u * z

        j *= weights

        return j

    w = 1. / np.sqrt(np.maximum(y, 1.))

    params = model.make_params(
        amp = init[0],
        x0 = init[1],
        sig = init[2],
    )

    result = model.fit(
        y, x=x, weights=w, params=params, fit_kws={"Dfun": jac, "col_deriv": True}
    )

    params = result.params

    return {
        "amp_1": SimpleFitParam(params["amp"].value, params["amp"].stderr),
        "x0_1": SimpleFitParam(params["x0"].value, params["x0"].stderr),
        "sig_1": SimpleFitParam(params["sig"].value, params["sig"].stderr),
    }


def erf_linear_background(
    x: npt.NDArray[np.floating],
    amp: float,
    x0: float,
    sig: float,
    lin: float,
    off: float,
) -> npt.NDArray[np.floating]:
    return amp * erf((x - x0) / sig) + lin * x + off


def erf_linear_1_gauss(
    x: npt.NDArray[np.floating],
    amp_1: float,
    x0_1: float,
    sig_1: float,
    erf_amp: float,
    lin: float,
    off: float,
) -> npt.NDArray[np.floating]:
    return (
        gauss(x, amp_1, x0_1, sig_1)
        + erf_linear_background(x, erf_amp, x0_1, sig_1, lin, off)
    )


def fit_erf_linear_1_gauss(
    x: npt.NDArray[np.floating],
    y: npt.NDArray[np.floating],
    init: list[float | np.floating],
) -> dict[str, SimpleFitParam]:
    model = Model(erf_linear_1_gauss)

    def jac(params, _, weights, x):
        amp_1 = params["amp_1"]
        x0_1 = params["x0_1"]
        sig_1 = params["sig_1"]

        j = np.empty((6, x.size), dtype=np.float64)

        e = gauss(x, amp_1, x0_1, sig_1)
        u = (x - x0_1) / sig_1
        z = u / sig_1

        j[0, :] = e / amp_1
        j[1, :] = e * z - TWO_OVER_SQRT_PI * np.exp(-u**2) / sig_1
        j[2, :] = e * u * z + TWO_OVER_SQRT_PI * np.exp(-u**2) * u / sig_1
        j[3, :] = erf(u)
        j[4, :] = 1.
        j[5, :] = 0.

        j *= weights

        return j

    w = 1. / np.sqrt(np.maximum(y, 1.))

    params = model.make_params(
        amp_1 = init[0],
        x0_1 = init[1],
        sig_1 = init[2],
        erf_amp = init[3],
        lin = init[4],
        off = init[5],
    )

    result = model.fit(
        y, x=x, weights=w, params=params, fit_kws={"Dfun": jac, "col_deriv": True}
    )

    params = result.params

    return {
        "amp_1": SimpleFitParam(params["amp_1"].value, params["amp_1"].stderr),
        "x0_1": SimpleFitParam(params["x0_1"].value, params["x0_1"].stderr),
        "sig_1": SimpleFitParam(params["sig_1"].value, params["sig_1"].stderr),
        "erf_amp": SimpleFitParam(params["erf_amp"].value, params["erf_amp"].stderr),
        "lin": SimpleFitParam(params["lin"].value, params["lin"].stderr),
        "const": SimpleFitParam(params["off"].value, params["off"].stderr),
    }


def erf_linear_2_gauss(
    x: npt.NDArray[np.floating],
    amp_1: float,
    x0_1: float,
    sig_1: float,
    amp_2: float,
    x0_2: float,
    sig_2: float,
    erf_amp: float,
    lin: float,
    off: float,
) -> npt.NDArray[np.floating]:
    return (
        gauss(x, amp_1, x0_1, sig_1)
        + gauss(x, amp_2, x0_2, sig_2)
        + erf_linear_background(x, erf_amp, x0_1, sig_1, lin, off)
    )


def fit_erf_linear_2_gauss(
    x: npt.NDArray[np.floating],
    y: npt.NDArray[np.floating],
    init: list[float | np.floating],
) -> dict[str, SimpleFitParam]:
    model = Model(erf_linear_1_gauss)

    def jac(params, _, weights, x):
        amp_1 = params["amp_1"]
        x0_1 = params["x0_1"]
        sig_1 = params["sig_1"]
        amp_2 = params["amp_2"]
        x0_2 = params["x0_2"]
        sig_2 = params["sig_2"]

        j = np.empty((9, x.size), dtype=np.float64)

        e_1 = gauss(x, amp_1, x0_1, sig_1)
        e_2 = gauss(x, amp_2, x0_2, sig_2)
        u_1 = (x - x0_1) / sig_1
        u_2 = (x - x0_2) / sig_2
        z_1 = u_1 / sig_1
        z_2 = u_2 / sig_2

        j[0, :] = e_1 / amp_1
        j[1, :] = e_1 * z_1 - TWO_OVER_SQRT_PI * np.exp(-u_1**2) / sig_1
        j[2, :] = e_1 * u_1 * z_1 + TWO_OVER_SQRT_PI * np.exp(-u_1**2) * u_1 / sig_1
        j[3, :] = e_2 / amp_2
        j[4, :] = e_2 * z_2
        j[5, :] = e_2 * u_2 * z_2
        j[6, :] = erf(u_1)
        j[7, :] = 1.
        j[8, :] = 0.

        j *= weights

        return j

    w = 1. / np.sqrt(np.maximum(y, 1.))

    params = model.make_params(
        amp_1 = init[0],
        x0_1 = init[1],
        sig_1 = init[2],
        amp_2 = init[3],
        x0_2 = init[4],
        sig_2 = init[5],
        erf_amp = init[6],
        lin = init[7],
        off = init[8],
    )

    result = model.fit(
        y, x=x, weights=w, params=params, fit_kws={"Dfun": jac, "col_deriv": True}
    )

    params = result.params

    return {
        "amp_1": SimpleFitParam(params["amp_1"].value, params["amp_1"].stderr),
        "x0_1": SimpleFitParam(params["x0_1"].value, params["x0_1"].stderr),
        "sig_1": SimpleFitParam(params["sig_1"].value, params["sig_1"].stderr),
        "amp_2": SimpleFitParam(params["amp_2"].value, params["amp_2"].stderr),
        "x0_2": SimpleFitParam(params["x0_2"].value, params["x0_2"].stderr),
        "sig_2": SimpleFitParam(params["sig_2"].value, params["sig_2"].stderr),
        "erf_amp": SimpleFitParam(params["erf_amp"].value, params["erf_amp"].stderr),
        "lin": SimpleFitParam(params["lin"].value, params["lin"].stderr),
        "const": SimpleFitParam(params["off"].value, params["off"].stderr),
    }


def erf_linear_3_gauss(
    x: npt.NDArray[np.floating],
    amp_1: float,
    x0_1: float,
    sig_1: float,
    amp_2: float,
    x0_2: float,
    sig_2: float,
    amp_3: float,
    x0_3: float,
    sig_3: float,
    erf_amp: float,
    lin: float,
    off: float,
) -> npt.NDArray[np.floating]:
    return (
        gauss(x, amp_1, x0_1, sig_1)
        + gauss(x, amp_2, x0_2, sig_2)
        + gauss(x, amp_3, x0_3, sig_3)
        + erf_linear_background(x, erf_amp, x0_1, sig_1, lin, off)
    )


def fit_erf_linear_3_gauss(
    x: npt.NDArray[np.floating],
    y: npt.NDArray[np.floating],
    init: list[float | np.floating],
) -> dict[str, SimpleFitParam]:
    model = Model(erf_linear_1_gauss)

    def jac(params, _, weights, x):
        amp_1 = params["amp_1"]
        x0_1 = params["x0_1"]
        sig_1 = params["sig_1"]
        amp_2 = params["amp_2"]
        x0_2 = params["x0_2"]
        sig_2 = params["sig_2"]
        amp_3 = params["amp_3"]
        x0_3 = params["x0_3"]
        sig_3 = params["sig_3"]

        j = np.empty((12, x.size), dtype=np.float64)

        e_1 = gauss(x, amp_1, x0_1, sig_1)
        e_2 = gauss(x, amp_2, x0_2, sig_2)
        e_3 = gauss(x, amp_3, x0_3, sig_3)
        u_1 = (x - x0_1) / sig_1
        u_2 = (x - x0_2) / sig_2
        u_3 = (x - x0_3) / sig_3
        z_1 = u_1 / sig_1
        z_2 = u_2 / sig_2
        z_3 = u_3 / sig_3

        j[0, :] = e_1 / amp_1
        j[1, :] = e_1 * z_1 - TWO_OVER_SQRT_PI * np.exp(-u_1**2) / sig_1
        j[2, :] = e_1 * u_1 * z_1 + TWO_OVER_SQRT_PI * np.exp(-u_1**2) * u_1 / sig_1
        j[3, :] = e_2 / amp_2
        j[4, :] = e_2 * z_2
        j[5, :] = e_1 * u_2 * z_2
        j[6, :] = e_3 / amp_3
        j[7, :] = e_2 * z_3
        j[8, :] = e_3 * u_3 * z_3
        j[9, :] = erf(u_1)
        j[10, :] = 1.
        j[11, :] = 0.

        j *= weights

        return j

    w = 1. / np.sqrt(np.maximum(y, 1.))

    params = model.make_params(
        amp_1 = init[0],
        x0_1 = init[1],
        sig_1 = init[2],
        amp_2 = init[3],
        x0_2 = init[4],
        sig_2 = init[5],
        amp_3 = init[6],
        x0_3 = init[7],
        sig_3 = init[8],
        erf_amp = init[9],
        lin = init[10],
        off = init[11],
    )

    result = model.fit(
        y, x=x, weights=w, params=params, fit_kws={"Dfun": jac, "col_deriv": True}
    )

    params = result.params

    return {
        "amp_1": SimpleFitParam(params["amp_1"].value, params["amp_1"].stderr),
        "x0_1": SimpleFitParam(params["x0_1"].value, params["x0_1"].stderr),
        "sig_1": SimpleFitParam(params["sig_1"].value, params["sig_1"].stderr),
        "amp_2": SimpleFitParam(params["amp_2"].value, params["amp_2"].stderr),
        "x0_2": SimpleFitParam(params["x0_2"].value, params["x0_2"].stderr),
        "sig_2": SimpleFitParam(params["sig_2"].value, params["sig_2"].stderr),
        "erf_amp": SimpleFitParam(params["erf_amp"].value, params["erf_amp"].stderr),
        "lin": SimpleFitParam(params["lin"].value, params["lin"].stderr),
        "const": SimpleFitParam(params["off"].value, params["off"].stderr),
    }
