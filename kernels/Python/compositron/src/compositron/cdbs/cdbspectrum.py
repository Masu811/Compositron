from copy import copy, deepcopy
from dataclasses import dataclass
from enum import Enum, auto

import numpy as np
import numpy.typing as npt

from .anti_aliasing import (
    Polygon, Rectangle, Parallelogram, Ellipse, Vertex, ConvexToVerticesCounterClockwise,
    get_bounding_box, agg_aa, intersect_convex_polygons
)
from ..constants import M_E_KEV
from ..core.fitting import SimpleFitParam
from ..core.utils import (
    EcalCorrectionOrder, EnergyDetector, Spectrum, to_min_type_spectrum,
    LinearCalibration, Unit, EnergyDetectorPair, KeV, Eres
)
from ..dbs import DBSpectrum
from .fitting import fit_gauss2d


class BackgroundModel(Enum):
    NONE = auto()


class Axis(Enum):
    FIRST_DET = auto()
    SECOND_DET = auto()


@dataclass
class DiagonalArea:
    width_cel: Unit
    width_cml: Unit
    offset_cel: Unit
    offset_cml: Unit


@dataclass
class AxisAlignedArea:
    first_det_bnds: tuple[float, float]
    second_det_bnds: tuple[float, float]


@dataclass
class EllipseArea:
    radius_cel: Unit
    radius_cml: Unit


Area = DiagonalArea | AxisAlignedArea | EllipseArea


@dataclass
class LinearProjectionBins:
    bin_width: Unit


@dataclass
class CustomProjectionBins:
    bins: npt.NDArray

    def __init__(self, bins: npt.ArrayLike) -> None:
        self.bins = np.asarray(bins)


ProjectionBins = LinearProjectionBins | CustomProjectionBins


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
        num = [
            DiagonalArea(KeV(2), Eres(1), KeV(0), KeV(0))
        ],
        denom = [
            DiagonalArea(KeV(np.inf), Eres(1), KeV(0), KeV(0))
        ],
        in_corrected_peak = True,
    ),
    LineshapeParamDefinition(
        name = "W",
        num = [
            DiagonalArea(KeV(1), Eres(1), KeV(3), KeV(0)),
            DiagonalArea(KeV(1), Eres(1), KeV(-3), KeV(0)),
        ],
        denom = [
            DiagonalArea(KeV(np.inf), Eres(1), KeV(0), KeV(0))
        ],
        in_corrected_peak = True,
    ),
    LineshapeParamDefinition(
        name = "P/T",
        num = [
            DiagonalArea(KeV(np.inf), Eres(1), KeV(0), KeV(0))
        ],
        denom = [
            AxisAlignedArea(
                first_det_bnds = (-np.inf, np.inf),
                second_det_bnds = (-np.inf, np.inf),
            )
        ],
        in_corrected_peak = False,
    ),
]


def ij_to_xy(i: float, ecal: LinearCalibration) -> float:
    return ecal.from_index(i)


def xy_to_ij(x: float, ecal: LinearCalibration) -> float:
    return ecal.to_index_float(x)


def xy_to_uv(x: float, y: float) -> tuple[float, float]:
    return (0.5 * (x - y), 0.5 * (x + y) - M_E_KEV)


def uv_to_xy(u: float, v: float) -> tuple[float, float]:
    return (M_E_KEV + u + v, M_E_KEV - u + v)


def ij_to_uv(
    i: float, j: float, ecal: tuple[LinearCalibration, LinearCalibration]
) -> tuple[float, float]:
    x, y = (ij_to_xy(j, ecal[1]), ij_to_xy(i, ecal[0]))
    return xy_to_uv(x, y)


def uv_to_ij(
    u: float, v: float, ecal: tuple[LinearCalibration, LinearCalibration]
) -> tuple[float, float]:
    x, y = uv_to_xy(u, v)
    return (xy_to_ij(y, ecal[0]), xy_to_ij(x, ecal[1]))


def convert_rectangle(
    first_det_bnds: tuple[float, float],
    second_det_bnds: tuple[float, float],
    ecal: tuple[LinearCalibration, LinearCalibration],
    peak_bnds: tuple[tuple[int, int], tuple[int, int]],
) -> Rectangle:
    return Rectangle(
        ecal[0].to_index_float(first_det_bnds[1]) - peak_bnds[0][0],
        ecal[0].to_index_float(first_det_bnds[0]) - peak_bnds[0][0],
        ecal[1].to_index_float(second_det_bnds[0]) - peak_bnds[1][0],
        ecal[1].to_index_float(second_det_bnds[1]) - peak_bnds[1][0],
    )


def convert_parallelogram(
    width_cel: float,
    width_cml: float,
    offset_cel: float,
    offset_cml: float,
    ecal: tuple[LinearCalibration, LinearCalibration],
    peak_bnds: tuple[tuple[int, int], tuple[int, int]],
) -> Parallelogram:
    i, j = uv_to_ij(offset_cel, offset_cml, ecal)
    w1 = np.sqrt(0.5 * (
        (width_cel / ecal[1].scale)**2 + (width_cel / ecal[0].scale)**2
    ))
    w2 = np.sqrt(0.5 * (
        (width_cml / ecal[1].scale)**2 + (-width_cml / ecal[0].scale)**2
    ))
    return Parallelogram(
        i - peak_bnds[0][0],
        j - peak_bnds[1][0],
        -ecal[1].scale / ecal[0].scale,
        ecal[1].scale / ecal[0].scale,
        w1,
        w2,
    )


def convert_ellipse(
    radius_cel: float,
    radius_cml: float,
    ecal: tuple[LinearCalibration, LinearCalibration],
    peak_bnds: tuple[tuple[int, int], tuple[int, int]],
) -> Ellipse:
    w1 = np.sqrt(0.5 * (
        (radius_cel / ecal[1].scale)**2 + (radius_cel / ecal[0].scale)**2
    ))
    w2 = np.sqrt(0.5 * (
        (radius_cml / ecal[1].scale)**2 + (-radius_cml / ecal[0].scale)**2
    ))
    return Ellipse(
        ecal[0].to_index_float(M_E_KEV) - peak_bnds[0][0],
        ecal[1].to_index_float(M_E_KEV) - peak_bnds[1][0],
        w1,
        w2,
        -np.atan(ecal[1].scale / ecal[0].scale),
    )


@dataclass
class FoldedProjection:
    energies: npt.NDArray[np.floating]
    spectrum: npt.NDArray[np.floating]


@dataclass
class Projection:
    spectrum: npt.NDArray[np.floating]
    parent_detector_name: str
    ecal: None | LinearCalibration
    bins: None | npt.NDArray[np.floating]
    counts: float
    dcounts: float

    def __init__(
        self,
        spectrum: npt.NDArray[np.floating],
        parent_detector_name: str,
        ecal: None | LinearCalibration,
        bins: None | npt.NDArray[np.floating],
    ) -> None:
        counts = np.sum(spectrum)
        dcounts = np.sqrt(counts)

        self.spectrum = spectrum
        self.parent_detector_name = parent_detector_name
        self.ecal = ecal
        self.bins = bins
        self.counts = float(counts)
        self.dcounts = dcounts

    def get_energies(self) -> npt.NDArray[np.floating]:
        if self.ecal is not None:
            return (
                np.arange(len(self.spectrum), dtype=np.float64) * self.ecal.scale
                + self.ecal.offset
            )
        elif self.bins is not None:
            return 0.5 * (self.bins[:-1] + self.bins[1:])

        raise ValueError("Projection is missing ecal and bins")

    def fold(self) -> FoldedProjection:
        if self.ecal is not None:
            center_idx = self.ecal.to_index_float(0)

            if abs(center_idx % 1 - 0.5) > 1e-2:
                raise ValueError(
                    "Could not fold projection as bins are not symmetric around 0"
                )

            last_left_idx = self.ecal.to_index_rounded(0)
            first_right_idx = last_left_idx + 1

            len_left = last_left_idx + 1
            len_right = len(self.spectrum) - len_left

            len_new = min(len_left, len_right)

            spectrum = np.zeros((len_new,), dtype=np.float64)

            spectrum += self.spectrum[first_right_idx:first_right_idx + len_new]
            spectrum += self.spectrum[last_left_idx - len_new + 1:last_left_idx + 1][::-1]

            energies = (
                np.arange(
                    first_right_idx,
                    first_right_idx + len_new,
                    dtype=np.float64
                ) * self.ecal.scale + self.ecal.offset
            )

            return FoldedProjection(energies, spectrum)
        elif self.bins is not None:
            center_idx = np.argmin(np.abs(self.bins))

            len_left = center_idx
            len_right = len(self.spectrum) - len_left

            len_new = min(len_left, len_right)

            spectrum, energies = [], []

            for i in range(len_new):
                left_bin_left_edge = self.bins[center_idx - i - 1]
                left_bin_right_edge = self.bins[center_idx - i]

                right_bin_left_edge = self.bins[center_idx + i]
                right_bin_right_edge = self.bins[center_idx + i + 1]

                if (
                    left_bin_left_edge + right_bin_right_edge > 1e-2 or
                    left_bin_right_edge + right_bin_left_edge > 1e-2
                ):
                    raise ValueError(
                        "Could not fold projection as bins are not symmetric around 0"
                    )

                spectrum.append(
                    self.spectrum[center_idx - i - 1] + self.spectrum[center_idx + i]
                )
                energies.append(
                    0.5 * (right_bin_left_edge + right_bin_right_edge)
                )

            return FoldedProjection(np.asarray(energies), np.asarray(spectrum))

        raise ValueError("Projection is missing ecal and bins")

    def to_dbspectrum(self) -> DBSpectrum:
        if self.ecal is None:
            raise ValueError(
                "Could not convert projection to DBSpectrum due to missing "
                "energy calibration"
            )

        return DBSpectrum(
            self.spectrum.astype(np.uint64),
            EnergyDetector(
                name = self.parent_detector_name,
                ecal = self.ecal,
                corrected_ecal = None,
                eres = None,
            )
        )


@dataclass
class CDBSpectrum:
    spectrum: Spectrum
    detpair: EnergyDetectorPair
    counts: int
    dcounts: float
    peak: None | npt.NDArray[np.floating]
    peak_bnds: None | tuple[tuple[int, int], tuple[int, int]]
    peak_counts: None | float
    dpeak_counts: None | float
    peak_params: dict[str, SimpleFitParam]
    lineshape_params: dict[str, LineshapeParam]

    def __init__(
        self, spectrum: npt.ArrayLike, detpair: EnergyDetectorPair,
    ) -> None:
        counts = np.sum(spectrum)
        dcounts = np.sqrt(counts)

        self.spectrum = to_min_type_spectrum(spectrum)
        self.detpair = detpair
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

        self.extract_peak((10., 10.), BackgroundModel.NONE)

        for param in STD_LINESHAPE_PARAMS:
            self.calc_lineshape_param(param)

        self.peak = None
        self.peak_bnds = None

    def project_axes(self, onto_axis: Axis) -> Projection:
        match onto_axis:
            case Axis.FIRST_DET:
                return Projection(
                    spectrum = np.atleast_1d(
                        np.sum(self.spectrum, axis=1, dtype=np.float64)
                    ),
                    parent_detector_name = self.detpair.first_det.name,
                    ecal = self.detpair.first_det.ecal,
                    bins = None,
                )
            case Axis.SECOND_DET:
                return Projection(
                    spectrum = np.atleast_1d(
                        np.sum(self.spectrum, axis=0, dtype=np.float64)
                    ),
                    parent_detector_name = self.detpair.second_det.name,
                    ecal = self.detpair.second_det.ecal,
                    bins = None,
                )

    def _fit_2d_peak(self, bg_model: BackgroundModel) -> None:
        if self.peak is None or self.peak_bnds is None:
            raise ValueError(
                "Attempting analysis on background subtracted peak, but no "
                "subtraction has been performed yet"
            )

        ecal_1 = (
            self.detpair.first_det.corrected_ecal or
            self.detpair.first_det.ecal
        )
        ecal_2 = (
            self.detpair.second_det.corrected_ecal or
            self.detpair.second_det.ecal
        )

        x = np.arange(*self.peak_bnds[1]) * ecal_2.scale + ecal_2.offset
        y = np.arange(*self.peak_bnds[0]) * ecal_1.scale + ecal_1.offset

        x, y = np.meshgrid(x, y)

        match bg_model:
            case BackgroundModel.NONE:
                argmax = np.argmax(self.peak)

                j0 = argmax % self.peak.shape[1]
                i0 = int(argmax // self.peak.shape[1])

                print(i0)
                print(j0)

                max_val = self.peak[i0, j0]

                print(max_val)

                y0 = ecal_1.from_index(i0 + self.peak_bnds[0][0])
                x0 = ecal_2.from_index(j0 + self.peak_bnds[1][0])

                print(x0)
                print(y0)

                self.peak_params = fit_gauss2d(
                    x, y, self.peak, [max_val, x0, y0, 1, 2, -0.78],
                )

    def _subtract_bg(self, bg_model: BackgroundModel) -> None:
        if bg_model == BackgroundModel.NONE:
            return

    def extract_peak(
        self, peak_window_kev: tuple[float, float], bg_model: BackgroundModel,
    ) -> None:
        ecal_1 = (
            self.detpair.first_det.corrected_ecal or
            self.detpair.first_det.ecal
        )
        ecal_2 = (
            self.detpair.second_det.corrected_ecal or
            self.detpair.second_det.ecal
        )

        first_row_e = ecal_1.to_index_float(M_E_KEV - peak_window_kev[0] / 2)
        last_row_e = ecal_1.to_index_float(M_E_KEV + peak_window_kev[0] / 2)

        first_col_e = ecal_2.to_index_float(M_E_KEV - peak_window_kev[1] / 2)
        last_col_e = ecal_2.to_index_float(M_E_KEV + peak_window_kev[1] / 2)

        first_row = min(max(int(first_row_e), 0), self.spectrum.shape[0] - 1)
        last_row = min(max(int(last_row_e) + 1, 0), self.spectrum.shape[0] - 1)

        first_col = min(max(int(first_col_e), 0), self.spectrum.shape[1] - 1)
        last_col = min(max(int(last_col_e) + 1, 0), self.spectrum.shape[1] - 1)

        self.peak = (
            self.spectrum[first_row:last_row, first_col:last_col]
        ).astype(np.float64)

        self.peak_bnds = ((first_row, last_row), (first_col, last_col))

        self._subtract_bg(bg_model)

        peak_counts = np.sum(self.peak)

        self.peak_counts = peak_counts
        self.dpeak_counts = np.sqrt(peak_counts)

    def correct_ecal(self, order: EcalCorrectionOrder) -> None:
        if order == EcalCorrectionOrder.NONE:
            return

        if self.peak is None:
            self.extract_peak((4., 4.), BackgroundModel.NONE)

        if "x0" not in self.peak_params or "y0" not in self.peak_params:
            self._fit_2d_peak(BackgroundModel.NONE)

        x0 = self.peak_params["x0"].val
        y0 = self.peak_params["y0"].val

        ecal_1 = copy(self.detpair.first_det.ecal)
        ecal_2 = copy(self.detpair.second_det.ecal)

        match order:
            case EcalCorrectionOrder.ZEROTH:
                ecal_1.offset += M_E_KEV - y0
                ecal_2.offset += M_E_KEV - x0
            case EcalCorrectionOrder.FIRST:
                ecal_1.scale *= (M_E_KEV - ecal_1.offset) / (y0 - ecal_1.offset)
                ecal_2.scale *= (M_E_KEV - ecal_2.offset) / (x0 - ecal_2.offset)

        self.detpair.first_det.corrected_ecal = ecal_1
        self.detpair.second_det.corrected_ecal = ecal_2

    def _get_max_projection_length(
        self,
        v: float,
        ecal: tuple[LinearCalibration, LinearCalibration],
        peak_bnds: tuple[tuple[int, int], tuple[int, int]],
    ) -> float:
        i1, j1 = uv_to_ij(0., v, ecal)
        i2, j2 = uv_to_ij(0., -v, ecal)
        m = -ecal[1].scale / ecal[0].scale
        t1 = i1 - m * j1
        t2 = i2 - m * j2

        ((i_min, i_max), (j_min, j_max)) = peak_bnds

        i_min = i_min - 0.49
        i_max = i_max + 0.49
        j_min = j_min - 0.49
        j_max = j_max + 0.49

        x1 = ij_to_uv(i_min, (i_min - t1) / m, ecal)
        x2 = ij_to_uv(i_min, (i_min - t2) / m, ecal)
        x3 = ij_to_uv(m * j_min + t1, j_min, ecal)
        x4 = ij_to_uv(m * j_min + t2, j_min, ecal)
        x5 = ij_to_uv(i_max, (i_max - t1) / m, ecal)
        x6 = ij_to_uv(i_max, (i_max - t2) / m, ecal)
        x7 = ij_to_uv(m * j_max + t1, j_max, ecal)
        x8 = ij_to_uv(m * j_max + t2, j_max, ecal)

        return min(abs(v[0]) for v in [x1, x2, x3, x4, x5, x6, x7, x8])

    def _calculate_polygon_boundaries(
        self,
        area: Area,
        ecal: tuple[LinearCalibration, LinearCalibration],
        peak_bnds: tuple[tuple[int, int], tuple[int, int]],
    ) -> ConvexToVerticesCounterClockwise:
        if isinstance(area, DiagonalArea):
            width_cel = area.width_cel.to_kev(self.detpair.eres)
            width_cml = area.width_cml.to_kev(self.detpair.eres)
            offset_cel = area.offset_cel.to_kev(self.detpair.eres)
            offset_cml = area.offset_cml.to_kev(self.detpair.eres)

            if np.isinf(width_cel):
                width_cel = self._get_max_projection_length(
                    0.5 * width_cml, ecal, peak_bnds
                )

            if any(
                not np.isfinite(x) for x in [
                    width_cel, width_cml, offset_cel, offset_cml
                ]
            ):
                raise ValueError("Invalid boundaries for area")

            return convert_parallelogram(
                width_cel, width_cml, offset_cel, offset_cml, ecal, peak_bnds,
            )
        elif isinstance(area, AxisAlignedArea):
            y_min, y_max = area.first_det_bnds
            x_min, x_max = area.second_det_bnds

            if np.isinf(y_min):
                y_min = ecal[0].from_index(peak_bnds[0][0] - 0.49)
            if np.isinf(y_max):
                y_max = ecal[0].from_index(peak_bnds[0][1] + 0.49)
            if np.isinf(x_min):
                x_min = ecal[1].from_index(peak_bnds[1][0] - 0.49)
            if np.isinf(x_max):
                x_max = ecal[1].from_index(peak_bnds[1][1] + 0.49)

            if any(
                not np.isfinite(x) for x in [y_min, y_max, x_min, x_max]
            ):
                raise ValueError("Invalid boundaries for area")

            return convert_rectangle(
                (y_min, y_max), (x_min, x_max), ecal, peak_bnds
            )
        elif isinstance(area, EllipseArea):
            radius_cel = area.radius_cel.to_kev(self.detpair.eres)
            radius_cml = area.radius_cml.to_kev(self.detpair.eres)

            if not np.isfinite(radius_cel) or not np.isfinite(radius_cml):
                raise ValueError("Invalid boundaries for area")

            return convert_ellipse(radius_cel, radius_cml, ecal, peak_bnds)

    def _blend_and_sum(
        self,
        weights: npt.NDArray[np.uint8],
        upper_row: int,
        left_col: int,
        nrows: int,
        ncols: int,
        in_corrected_peak: bool,
    ) -> float:
        if in_corrected_peak:
            if self.peak is None:
                raise ValueError(
                    "Attempting analysis on background subtracted peak, but no "
                    "subtraction has been performed yet"
                )

            view = self.peak[upper_row:upper_row+nrows, left_col:left_col+ncols]

            return np.sum(view * weights, dtype=np.float64) / 255.
        else:
            view = self.spectrum[upper_row:upper_row+nrows, left_col:left_col+ncols]

            import matplotlib.pyplot as plt
            from matplotlib.colors import LogNorm

            plt.imshow(view, origin="lower", norm=LogNorm(vmin=1), cmap="jet")
            plt.show()

            return np.sum(view.astype(np.uint64) * weights, dtype=np.uint64) / 255.

    def integrate(self, area: Area, in_corrected_peak: bool) -> float:
        ecal_1 = (
            self.detpair.first_det.corrected_ecal or
            self.detpair.first_det.ecal
        )
        ecal_2 = (
            self.detpair.second_det.corrected_ecal or
            self.detpair.second_det.ecal
        )
        ecal = (ecal_1, ecal_2)

        if in_corrected_peak:
            if self.peak_bnds is None:
                raise ValueError(
                    "Attempting analysis on background subtracted peak, but no "
                    "subtraction has been performed yet"
                )
            integral_bnds = self.peak_bnds
        else:
            integral_bnds = (
                (0, self.spectrum.shape[0] - 1), (0, self.spectrum.shape[1] - 1)
            )

        nrows = integral_bnds[0][1] - integral_bnds[0][0] + 1
        ncols = integral_bnds[1][1] - integral_bnds[1][0] + 1

        polygon = self._calculate_polygon_boundaries(area, ecal, integral_bnds)

        vertices = polygon.to_vertices()

        left_col, right_col, lower_row, upper_row = get_bounding_box(vertices)

        if right_col > ncols - 1 or lower_row > nrows - 1:
            raise ValueError("Polygon vertices are out of matrix bounds")

        vertices = [Vertex(v.x - left_col, v.y - upper_row) for v in vertices]

        nrows_view = lower_row - upper_row + 1
        ncols_view = right_col - left_col + 1

        weights = agg_aa(nrows_view, ncols_view, Polygon(vertices))

        return self._blend_and_sum(
            weights,
            upper_row,
            left_col,
            nrows_view,
            ncols_view,
            in_corrected_peak
        )

    def calc_lineshape_param(
        self, definition: LineshapeParamDefinition
    ) -> LineshapeParam:
        ecal_1 = (
            self.detpair.first_det.corrected_ecal or
            self.detpair.first_det.ecal
        )
        ecal_2 = (
            self.detpair.second_det.corrected_ecal or
            self.detpair.second_det.ecal
        )
        ecal = (ecal_1, ecal_2)

        if definition.in_corrected_peak:
            if self.peak_bnds is None:
                raise ValueError(
                    "Attempting analysis on background subtracted peak, but no "
                    "subtraction has been performed yet"
                )
            integral_bnds = self.peak_bnds
        else:
            integral_bnds = (
                (0, self.spectrum.shape[0] - 1), (0, self.spectrum.shape[1] - 1)
            )

        num_polygons = [
            self._calculate_polygon_boundaries(area, ecal, integral_bnds)
            for area in definition.num
        ]
        denom_polygons = [
            self._calculate_polygon_boundaries(area, ecal, integral_bnds)
            for area in definition.denom
        ]

        for field in [num_polygons, denom_polygons]:
            for i in range(len(field) - 1):
                for j in range(i + 1, len(field)):
                    intersection = intersect_convex_polygons(
                        field[i], field[j]
                    )

                    if len(intersection.vertices) > 0:
                        raise ValueError(
                            "Lineshape numerator or denominator areas overlap"
                        )

        num_areas = [
            self.integrate(area, definition.in_corrected_peak)
            for area in definition.num
        ]
        denom_areas = [
            self.integrate(area, definition.in_corrected_peak)
            for area in definition.denom
        ]

        common_areas = []

        for num_area in num_polygons:
            for denom_area in denom_polygons:
                intersection = intersect_convex_polygons(
                    num_area, denom_area
                )

                if len(intersection.vertices) == 0:
                    continue

                left_col, right_col, lower_row, upper_row = get_bounding_box(
                    intersection.to_vertices()
                )

                intersection.vertices = [
                    Vertex(v.x - left_col, v.y - upper_row)
                    for v in intersection.vertices
                ]

                nrows_view = lower_row - upper_row + 1
                ncols_view = right_col - left_col + 1

                weights = agg_aa(nrows_view, ncols_view, intersection)

                integral = self._blend_and_sum(
                    weights,
                    upper_row,
                    left_col,
                    nrows_view,
                    ncols_view,
                    definition.in_corrected_peak,
                )

                common_areas.append(integral)

        num = sum(num_areas)
        denom = sum(denom_areas)
        common = sum(common_areas)

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

    def project_digonal(
        self,
        bins: ProjectionBins,
        width: Unit,
        max_length: Unit,
        in_corrected_peak: bool,
    ) -> Projection:
        ecal_1 = (
            self.detpair.first_det.corrected_ecal or
            self.detpair.first_det.ecal
        )
        ecal_2 = (
            self.detpair.second_det.corrected_ecal or
            self.detpair.second_det.ecal
        )
        ecal = (ecal_1, ecal_2)

        if in_corrected_peak:
            if self.peak_bnds is None:
                raise ValueError(
                    "Attempting analysis on background subtracted peak, but no "
                    "subtraction has been performed yet"
                )
            integral_bnds = self.peak_bnds
        else:
            integral_bnds = (
                (0, self.spectrum.shape[0] - 1), (0, self.spectrum.shape[1] - 1)
            )

        width_kev = width.to_kev(self.detpair.eres)
        max_length_kev = max_length.to_kev(self.detpair.eres)

        if isinstance(bins, LinearProjectionBins):
            u_bin = bins.bin_width.to_kev(self.detpair.eres)
            u_max = self._get_max_projection_length(
                width_kev / 2, ecal, integral_bnds
            )

            u_max = min(u_max, max_length_kev)

            num_bins = int(2 * u_max / u_bin)

            if num_bins < 2:
                raise ValueError("Invalid bins for projection")

            max_offset = (num_bins / 2 - 0.5) * u_bin

            ecal = LinearCalibration(
                offset = -max_offset, scale = u_bin
            )

            spectrum = [
                self.integrate(
                    DiagonalArea(
                        KeV(u_bin), width, KeV(ecal.from_index(i)), KeV(0)
                    ),
                    in_corrected_peak
                ) for i in range(num_bins)
            ]

            return Projection(
                np.asarray(spectrum),
                parent_detector_name = self.detpair.name,
                ecal = ecal,
                bins = None,
            )

        elif isinstance(bins, CustomProjectionBins):
            if len(bins.bins) < 2:
                raise ValueError("Invalid bins for projection")

            bins_kev = [x.to_kev(self.detpair.eres) for x in bins.bins]

            spectrum = np.zeros(len(bins_kev) - 1, dtype=np.float64)

            for i in range(len(bins_kev) - 1):
                left_edge = bins_kev[i]
                right_edge = bins_kev[i + 1]

                width_cel = right_edge - left_edge
                offset = 0.5 * (left_edge + right_edge)

                spectrum[i] = self.integrate(
                    DiagonalArea(
                        KeV(width_cel), width, KeV(offset), KeV(0)
                    ),
                    in_corrected_peak
                )

            return Projection(
                np.asarray(spectrum),
                parent_detector_name = self.detpair.name,
                ecal = None,
                bins = np.asarray(bins_kev),
            )
