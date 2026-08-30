from dataclasses import dataclass

import numpy as np
import numpy.typing as npt

from ..constants import FWHM_OVER_SIGMA
from ..core.utils import Spectrum, TimingDetectorPair, to_min_type_spectrum
from ..core.fitting import FitParam, SimpleFitParam
from ..dbs.fitting import fit_gauss
from .fitting import fit_lifetime_spectrum, FitResult
from .model import LifetimeModel


@dataclass
class PALSpectrum:
    spectrum: Spectrum
    detpair: TimingDetectorPair
    counts: int
    peak: None | npt.NDArray[np.floating]
    peak_bnd_idcs: None | tuple[int, int]
    peak_counts: None | float
    peak_params: dict[str, SimpleFitParam]
    fit_result: None | FitResult


    def __init__(
        self, spectrum: npt.ArrayLike, detpair: TimingDetectorPair,
    ) -> None:
        self.spectrum = to_min_type_spectrum(spectrum)
        self.detpair = detpair
        self.counts = np.sum(spectrum, dtype=int)
        self.peak_bnd_idcs = None
        self.peak_counts = None
        self.peak_params = {}
        self.fit_result = None


    def _get_peak_center(self) -> float:
        argmax = np.argmax(self.spectrum)

        peak_center_window = 15

        roi = argmax - peak_center_window, argmax + peak_center_window

        y = self.spectrum[roi[0]:roi[1]].astype(float)
        x = np.arange(*roi, dtype=float)

        peak_params = fit_gauss(x, y, [self.spectrum[argmax], float(argmax), 100])

        peak_height = peak_params["amp_1"]
        peak_center = peak_params["x0_1"]
        peak_width = peak_params["sig_1"]

        self.peak_params["center_idx"] = peak_center
        self.peak_params["height"] = peak_height
        self.peak_params["fwhm"] = SimpleFitParam(
            peak_width.val * FWHM_OVER_SIGMA, peak_width.err * FWHM_OVER_SIGMA,
        )

        return peak_center.val


    def fit(
        self,
        model: LifetimeModel,
        fit_start_idx: int,
        fit_end_idx: int,
    ) -> None:
        if "center_idx" not in self.peak_params:
            self._get_peak_center()

        peak_center = self.peak_params["center_idx"].val

        y = self.spectrum[fit_start_idx:fit_end_idx].astype(float)

        tcal = self.detpair.corrected_tcal or self.detpair.tcal

        x = tcal.from_index(np.arange(fit_start_idx, fit_end_idx))

        self.fit_result = fit_lifetime_spectrum(
            x, y, model, self.counts, tcal.from_index(peak_center)
        )


    def fit_report(self, sep: str = "") -> None:
        """Print a fit report, PALSFIT style.

        Paramters
        ---------
        sep : str, optional
            Column separator character. Default is "".
        """
        err = "No fit has been performed yet"
        assert self.fit_result is not None, err

        result = self.fit_result
        model = result.model

        n_l = len(model.lifetime_components)
        n_r = len(model.resolution_components)

        headers = ["Parameter", "Value", "Stderr (abs)", "Stderr (rel)",
                   "Min", "Max", "Fixed", ""]
        formats = ["<", ">", ">", ">", ">", ">", "<", "<"]

        lt_norm = sum(l.intensity.val for l in model.lifetime_components)
        res_norm = sum(r.intensity.val for r in model.resolution_components)

        lt_norm_okay = abs(lt_norm - 1) < 1e-3
        res_norm_okay = abs(res_norm - 1) < 1e-3

        def print_params(params: list[FitParam], names: list[str]):
            table_data = np.array([
                [
                    name,
                    (
                        f"{100 * p.val:.1f}" if "int" in name else
                        f"{p.val:.2f}" if "life" in name else
                        f"{p.val:.4g}"
                    ),
                    (
                        (f"{100 * p.err:.1f}" if "int" in name else f"{p.err:.4g}")
                        if p.err is not None else "nan"
                    ),
                    f"{abs(100 * p.err / p.val):.3f} %" if p.val != 0 and p.err is not None else "nan",
                    f"{100 * p.min:.1f}" if "int" in name else f"{p.min:.1f}",
                    f"{100 * p.max:.1f}" if "int" in name else f"{p.max:.1f}",
                    "" if p.vary else "Fixed",
                    "At Boundary" if abs(p.val- p.min) < 1e-3 or abs(p.val- p.max) < 1e-3 else "",
                ]
                for name, p in zip(names, params)
            ], dtype=str)

            widths = np.max([np.max(np.vectorize(len)(table_data), axis=0),
                            np.vectorize(len)(headers)], axis=0)

            for i, head in enumerate(headers):
                print(f"{sep} {head:^{widths[i]}} ", end="")
            print(sep)
            for i, head in enumerate(headers):
                print(f"{sep} {'-' * widths[i]} ", end="")
            print(sep)

            for row in table_data:
                for col, fmt, wid in zip(row, formats, widths):
                    print(f"{sep} {col:{fmt}{wid}} ", end="")
                print(sep)

        print("#"*100)
        print(f"{' Fit report ':^100}")
        print("#"*100)
        print()
        print(f"Statistics: {self.counts:,}")
        print(f"Number of Data Points: {result.n_dpoints}")
        print()
        print(f"Lifetime Components:   {n_l}")
        print(f"Resolution Components: {n_r}")
        print(f"Background Components: {int(model.background_component is not None)}")
        print()
        print(f"Fit Method: Levenberg-Marquardt")
        print(f"Fit Status: {'Success' if result.fit_status else 'Failure'}")
        print(f"Number of Function Evaluations: {result.n_eval}")
        print(f"Reduced Chi squared: {result.red_chi_2}")
        print()
        print("-"*100)
        print(f"{'Parameters':^100}")
        print("-"*100)
        print()
        print(f"{' General Parameters ':=^100}")
        print()

        print_params(
            [model.scale_component, model.shift_component, model.background_component],
            ["N", "t0", "background", "mean_lifetime"],
        )

        print()
        print(f"{' Lifetime Components ':=^100}")
        print()

        if not lt_norm_okay:
            print("!!! intensities don't add up to 100% !!!")
            print()

        params = {
            f"lifetime_{i}": l_comp.lifetime
            for i, l_comp in enumerate(model.lifetime_components, 1)
        } | {
            f"intensity_{i}": l_comp.intensity
            for i, l_comp in enumerate(model.lifetime_components, 1)
        }
        print_params(list(params.values()), list(params.keys()))

        print()
        print(f"{' Resolution Components ':=^100}")
        print()

        if not res_norm_okay:
            print("!!! intensities don't add up to 100% !!!")
            print()

        params = {
            f"fwhm_{i}": r_comp.fwhm
            for i, r_comp in enumerate(model.resolution_components, 1)
        } | {
            f"intensity_{i}": r_comp.intensity
            for i, r_comp in enumerate(model.resolution_components, 1)
        } | {
            f"t0_{i}": r_comp.t0
            for i, r_comp in enumerate(model.resolution_components, 1)
        }
        print_params(list(params.values()), list(params.keys()))

        print()
        print("-"*100)
