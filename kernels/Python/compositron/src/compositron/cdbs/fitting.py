import numpy as np
import numpy.typing as npt
from lmfit import Model

from ..core.fitting import SimpleFitParam


def gauss2d(
    x: npt.NDArray[np.floating],
    y: npt.NDArray[np.floating],
    amp: float,
    x0: float,
    y0: float,
    sig_x: float,
    sig_y: float,
    phi: float,
) -> npt.NDArray[np.floating]:
    cos_phi = np.cos(phi)
    sin_phi = np.sin(phi)

    dx = x - x0
    dy = y - y0

    return amp * np.exp(
        -0.5 * ((dx * cos_phi - dy * sin_phi) / sig_x) ** 2
        -0.5 * ((dx * sin_phi + dy * cos_phi) / sig_y) ** 2
    )


def fit_gauss2d(
    x: npt.NDArray[np.floating],
    y: npt.NDArray[np.floating],
    z: npt.NDArray[np.floating],
    init: list[float],
) -> dict[str, SimpleFitParam]:
    model = Model(gauss2d, independent_vars=["x", "y"])

    w = 1. / np.sqrt(np.maximum(z, 1.))

    params = model.make_params(
        amp = init[0],
        x0 = init[1],
        y0 = init[2],
        sig_x = init[3],
        sig_y = init[4],
        phi = init[5],
    )

    result = model.fit(
        z, x=x, y=y, weights=w, params=params,
        #fit_kws={"Dfun": jac, "col_deriv": True}
    )

    params = result.params

    return {
        "amp": SimpleFitParam(params["amp"].value, params["amp"].stderr),
        "x0": SimpleFitParam(params["x0"].value, params["x0"].stderr),
        "y0": SimpleFitParam(params["y0"].value, params["y0"].stderr),
        "sig_x": SimpleFitParam(params["sig_x"].value, params["sig_x"].stderr),
        "sig_y": SimpleFitParam(params["sig_y"].value, params["sig_y"].stderr),
        "phi": SimpleFitParam(params["phi"].value, params["phi"].stderr),
    }
