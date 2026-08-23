from copy import copy, deepcopy
from dataclasses import dataclass
from enum import Enum, auto

import numpy as np
import numpy.typing as npt

from ..constants import M_E_KEV
from ..core.fitting import SimpleFitParam
from ..core.utils import (
    EcalCorrectionOrder, to_min_type_spectrum, Spectrum, EnergyDetector, Unit,
    KeV
)
from .fitting import (
    erf_linear_background, fit_gaussian, fit_erf_linear_1_gauss,
    fit_erf_linear_2_gauss, fit_erf_linear_3_gauss,
)


class PeakModel(Enum):
    GAUSS = auto()
    ERFLINEAR1GAUSS = auto()
    ERFLINEAR2GAUSS = auto()
    ERFLINEAR3GAUSS = auto()


class BinIntegrationScheme(Enum):
    CONST = auto()
    LINEAR = auto()


@dataclass
class PeakArea:
    width: Unit

    def to_bnds_kev(self, eres: None | float) -> tuple[float, float]:
        width_kev = self.width.to_kev(eres)
        return (M_E_KEV - width_kev / 2, M_E_KEV + width_kev / 2)


@dataclass
class PeakAreaWithOffset:
    width: Unit
    offset: Unit

    def to_bnds_kev(self, eres: None | float) -> tuple[float, float]:
        width_kev = self.width.to_kev(eres)
        offset_kev = self.offset.to_kev(eres)
        return (
            M_E_KEV - width_kev / 2 + offset_kev,
            M_E_KEV + width_kev / 2 + offset_kev
        )


@dataclass
class SpectrumArea:
    left_bnd: float
    right_bnd: float

    def to_bnds_kev(self, _: None | float) -> tuple[float, float]:
        return (self.left_bnd, self.right_bnd)


Area = PeakArea | PeakAreaWithOffset | SpectrumArea


@dataclass
class LineshapeParam:
    val: float
    err: float
    num: list[Area]
    denom: list[Area]
    in_corrected_peak: bool


@dataclass
class LineshapeParamDefinition:
    name: str
    num: list[Area]
    denom: list[Area]
    in_corrected_peak: bool


STD_LINESHAPE_PARAMS: list[LineshapeParamDefinition] = [
    LineshapeParamDefinition(
        name = "S",
        num = [PeakArea(KeV(1))],
        denom = [PeakArea(KeV(20))],
        in_corrected_peak = True,
    ),
    LineshapeParamDefinition(
        name = "W",
        num = [
            PeakAreaWithOffset(KeV(1), KeV(3)),
            PeakAreaWithOffset(KeV(1), KeV(-3))
        ],
        denom = [PeakArea(KeV(20))],
        in_corrected_peak = True,
    ),
    LineshapeParamDefinition(
        name = "V/P",
        num = [SpectrumArea(400, 500)],
        denom = [SpectrumArea(501, 521)],
        in_corrected_peak = False,
    ),
    LineshapeParamDefinition(
        name = "P/T",
        num = [SpectrumArea(501, 521)],
        denom = [SpectrumArea(-np.inf, np.inf)],
        in_corrected_peak = False,
    ),
]


@dataclass
class DBSpectrum:
    spectrum: Spectrum
    detector: EnergyDetector
    counts: int
    dcounts: float
    peak: None | npt.NDArray[np.floating]
    peak_bnds: None | tuple[int, int]
    peak_counts: None | float
    dpeak_counts: None | float
    peak_params: dict[str, SimpleFitParam]
    lineshape_params: dict[str, LineshapeParam]

    def __init__(
        self, spectrum: npt.ArrayLike, detector: EnergyDetector,
    ) -> None:
        counts = np.sum(spectrum)
        dcounts = np.sqrt(counts)

        self.spectrum = to_min_type_spectrum(spectrum)
        self.detector = detector
        self.counts = counts
        self.dcounts = dcounts
        self.peak = None
        self.peak_bnds = None
        self.peak_counts = None
        self.dpeak_counts = None
        self.peak_params = {}
        self.lineshape_params = {}

    def default_analyze(self) -> None:
        self.correct_ecal(EcalCorrectionOrder.FIRST)

        self.extract_peak(30, PeakModel.ERFLINEAR1GAUSS)

        for param in STD_LINESHAPE_PARAMS:
            self.calc_lineshape_param(param, BinIntegrationScheme.CONST)

        self.peak = None
        self.peak_bnds = None

    def _get_peak_energies(self) -> npt.NDArray[np.float64]:
        ecal = self.detector.corrected_ecal or self.detector.ecal

        assert (bnds := self.peak_bnds) is not None

        return (
            np.arange(bnds[0], bnds[1], dtype=np.float64) * ecal.scale
            + ecal.offset
        )

    def _fit_peak(self, peak_model: PeakModel) -> None:
        x = self._get_peak_energies()
        assert (y := self.peak) is not None

        m = np.max(y)

        match peak_model:
            case PeakModel.GAUSS:
                self.peak_params = fit_gaussian(x, y, [m, 511, 1])
            case PeakModel.ERFLINEAR1GAUSS:
                self.peak_params = fit_erf_linear_1_gauss(
                    x, y, [m, 511, 1, -10, 0, 0]
                )
            case PeakModel.ERFLINEAR2GAUSS:
                self.peak_params = fit_erf_linear_2_gauss(
                    x, y, [m / 2, 511, 1, m / 2, 511, 1.5, -10, 0, 0]
                )
            case PeakModel.ERFLINEAR3GAUSS:
                self.peak_params = fit_erf_linear_3_gauss(
                    x, y,
                    [m / 3, 511, 1, m / 3, 511, 1.5, m / 3, 511, 2, -10, 0, 0]
                )

    def _subtract_bg(self, peak_model: PeakModel) -> None:
        if peak_model == PeakModel.GAUSS:
            return

        self._fit_peak(peak_model)

        x = self._get_peak_energies()
        assert (y := self.peak) is not None

        corr = erf_linear_background(
            x,
            self.peak_params["erf_amp"].val,
            self.peak_params["x0_1"].val,
            self.peak_params["sig_1"].val,
            self.peak_params["lin"].val,
            self.peak_params["const"].val,
        )

        y -= corr

    def extract_peak(
        self, peak_width: float, peak_model: PeakModel
    ) -> npt.NDArray[np.float64]:
        ecal = self.detector.corrected_ecal or self.detector.ecal

        left_peak_idx = ecal.to_index_rounded(M_E_KEV - peak_width / 2)
        right_peak_idx = ecal.to_index_rounded(M_E_KEV + peak_width / 2)

        self.peak_bnds = (left_peak_idx, right_peak_idx)

        self.peak = (
            self.spectrum[left_peak_idx:right_peak_idx].astype(np.float64)
        )

        self._subtract_bg(peak_model)

        return self.peak

    def correct_ecal(self, order: EcalCorrectionOrder) -> None:
        if order == EcalCorrectionOrder.NONE:
            return

        if self.peak is None:
            self.extract_peak(20, PeakModel.GAUSS)

        self._fit_peak(PeakModel.GAUSS)

        c = self.peak_params["x0_1"].val

        ecal = copy(self.detector.ecal)

        match order:
            case EcalCorrectionOrder.ZEROTH:
                ecal.offset += M_E_KEV - c
            case EcalCorrectionOrder.FIRST:
                ecal.scale *= (M_E_KEV - ecal.offset) / (c - ecal.offset);

        self.detector.corrected_ecal = ecal

    @staticmethod
    def _integrate_const(
        spectrum: npt.NDArray, lower_ch_bnd: float, upper_ch_bnd: float
    ) -> float:
        counts = 0.

        lowest_ch = int(round(lower_ch_bnd))
        highest_ch = int(round(upper_ch_bnd))

        counts += np.sum(spectrum[lowest_ch:highest_ch+1], dtype=np.float64)

        # Left edge counts

        lowest_ch_counts = float(spectrum[lowest_ch])
        x1 = lowest_ch - 0.5
        counts -= lowest_ch_counts * (lower_ch_bnd - x1)

        # Right edge counts

        highest_ch_counts = float(spectrum[highest_ch])
        x2 = highest_ch + 0.5
        counts -= highest_ch_counts * (x2 - upper_ch_bnd)

        return counts

    @staticmethod
    def _integrate_linear(
        spectrum: npt.NDArray, lower_ch_bnd: float, upper_ch_bnd: float
    ) -> float:
        counts = 0.

        lowest_ch = int(np.floor(lower_ch_bnd))
        highest_ch = int(np.ceil(upper_ch_bnd))

        counts += np.sum(spectrum[lowest_ch+1:highest_ch], dtype=np.float64)

        counts += 0.5 * spectrum[lowest_ch]
        counts += 0.5 * spectrum[highest_ch]

        # Left edge counts

        x1 = float(lowest_ch)

        y1 = float(spectrum[lowest_ch])
        y2 = float(spectrum[lowest_ch + 1])

        y_l = (y2 - y1) * (lower_ch_bnd - x1) + y1

        counts -= 0.5 * (y1 + y_l) * (lower_ch_bnd - x1)

        # Left edge counts

        x1 = float(highest_ch - 1)
        x2 = float(highest_ch)

        y1 = float(spectrum[highest_ch - 1])
        y2 = float(spectrum[highest_ch])

        y_u = (y2 - y1) * (upper_ch_bnd - x1) + y1

        counts -= 0.5 * (y_u + y2) * (x2 - upper_ch_bnd)

        return counts

    @staticmethod
    def _integrate_generic(
        spectrum: npt.NDArray,
        lower_ch_bnd: float,
        upper_ch_bnd: float,
        bin_integration_scheme: BinIntegrationScheme,
    ) -> float:
        if np.isinf(lower_ch_bnd):
            lower_ch_bnd = 0.

        if np.isinf(upper_ch_bnd):
            upper_ch_bnd = float(len(spectrum) - 1)

        if (
            lower_ch_bnd < 0 or lower_ch_bnd > len(spectrum) - 1 or
            upper_ch_bnd < 0 or upper_ch_bnd > len(spectrum) - 1
        ):
            raise ValueError(
                "Attempting to integrate an area larger than the spectrum or "
                "extracted peak"
            )

        match bin_integration_scheme:
            case BinIntegrationScheme.CONST:
                return DBSpectrum._integrate_const(
                    spectrum, lower_ch_bnd, upper_ch_bnd
                )
            case BinIntegrationScheme.LINEAR:
                return DBSpectrum._integrate_linear(
                    spectrum, lower_ch_bnd, upper_ch_bnd
                )

    def integrate(
        self,
        area: Area,
        in_corrected_peak: bool,
        bin_integration_scheme: BinIntegrationScheme,
    ) -> float:
        lower_bnd, upper_bnd = area.to_bnds_kev(self.detector.eres)

        if in_corrected_peak and self.peak is None:
            raise RuntimeError(
                "Attempting analysis on background subtracted peak, but no "
                "subtraction has been performed yet"
            )

        if lower_bnd >= upper_bnd:
            raise ValueError(
                "Input parameter violates constraint "
                f"Lower bound < upper bound (got {lower_bnd} < {upper_bnd})"
            )

        ecal = self.detector.corrected_ecal or self.detector.ecal

        lower_ch_bnd = ecal.to_index_float(lower_bnd)
        upper_ch_bnd = ecal.to_index_float(upper_bnd)

        if in_corrected_peak:
            assert (peak_bnds := self.peak_bnds) is not None
            lower_peak_bnd, _ = peak_bnds

            lower_ch_bnd -= lower_peak_bnd
            upper_ch_bnd -= lower_peak_bnd

            assert (peak := self.peak) is not None

            return DBSpectrum._integrate_generic(
                peak, lower_ch_bnd, upper_ch_bnd, bin_integration_scheme,
            )
        else:
            return DBSpectrum._integrate_generic(
                self.spectrum,
                lower_ch_bnd,
                upper_ch_bnd,
                bin_integration_scheme,
            )

    def _get_area_bnds(self, areas: list[Area]) -> list[tuple[float, float]]:
        bnds = [area.to_bnds_kev(self.detector.eres) for area in areas]

        for i in range(len(areas) - 1):
            for j in range(i + 1, len(areas)):
                if not (bnds[i][1] < bnds[j][0] or bnds[j][1] < bnds[i][0]):
                    raise ValueError(
                        "Lineshape numerator or denominator areas overlap"
                    )

        return bnds

    def calc_lineshape_param(
        self,
        definition: LineshapeParamDefinition,
        bin_integration_scheme: BinIntegrationScheme,
    ) -> LineshapeParam:
        num_bnds = self._get_area_bnds(definition.num)
        denom_bnds = self._get_area_bnds(definition.denom)

        num = sum([
            self.integrate(
                area, definition.in_corrected_peak, bin_integration_scheme
            ) for area in definition.num
        ])
        denom = sum([
            self.integrate(
                area, definition.in_corrected_peak, bin_integration_scheme
            ) for area in definition.denom
        ])

        common = 0.

        for num_bnd in num_bnds:
            for denom_bnd in denom_bnds:
                left_bnd = max(num_bnd[0], denom_bnd[0])
                right_bnd = min(num_bnd[1], denom_bnd[1])

                if left_bnd < right_bnd:
                    common += self.integrate(
                        SpectrumArea(left_bnd, right_bnd),
                        definition.in_corrected_peak,
                        bin_integration_scheme,
                    )

        val = num / denom
        err = val * np.sqrt(
            (num - common) / num**2 +
            (denom - common) / denom**2 +
            common * ((denom - num) / (denom * num))**2
        )

        ls_param = LineshapeParam(
            val,
            err,
            deepcopy(definition.num),
            deepcopy(definition.denom),
            definition.in_corrected_peak,
        )

        self.lineshape_params[definition.name] = ls_param

        return ls_param
