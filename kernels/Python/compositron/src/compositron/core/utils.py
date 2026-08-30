from dataclasses import dataclass
from enum import Enum, auto
from typing import overload

import numpy as np
import numpy.typing as npt

from ..constants import KEV_PER_M0C


@dataclass
class LinearCalibration:
    offset: float
    scale: float

    @overload
    def from_index(self, index: float) -> float: ...

    @overload
    def from_index(self, index: np.ndarray) -> np.ndarray: ...

    def from_index(self, index: float | np.ndarray) -> float | np.ndarray:
        return self.scale * index + self.offset

    def to_index_rounded(self, value) -> int:
        return int(round((value - self.offset) / self.scale))

    def to_index_trunc(self, value: float) -> int:
        return int((value - self.offset) / self.scale)

    def to_index_float(self, value: float) -> float:
        return (value - self.offset) / self.scale


@dataclass
class EnergyDetector:
    name: str
    ecal: LinearCalibration
    corrected_ecal: None | LinearCalibration
    eres: None | float


@dataclass
class EnergyDetectorPair:
    name: str
    first_det: EnergyDetector
    second_det: EnergyDetector
    eres: None | float


@dataclass
class TimingDetectorPair:
    name: str
    tcal: LinearCalibration
    corrected_tcal: None | LinearCalibration
    tres: None | float


class EcalCorrectionOrder(Enum):
    ZEROTH = auto()
    FIRST = auto()
    NONE = auto()


@dataclass
class KeV:
    value: float

    def to_kev(self, _: None | float) -> float:
        return self.value


@dataclass
class M0C:
    value: float

    def to_kev(self, _: None | float) -> float:
        return self.value * KEV_PER_M0C


@dataclass
class Eres:
    value: float

    def to_kev(self, eres: None | float) -> float:
        if eres is None:
            raise ValueError(
                "Could not convert units of eres to keV due to missing eres"
            )

        return self.value * eres


Unit = KeV | M0C | Eres


Spectrum = npt.NDArray[np.unsignedinteger]


def to_min_type_spectrum(spectrum: npt.ArrayLike) -> Spectrum:
    max = np.max(spectrum)

    if max <= np.iinfo(np.uint8).max:
        return np.asarray(spectrum, dtype=np.uint8)
    elif max <= np.iinfo(np.uint16).max:
        return np.asarray(spectrum, dtype=np.uint16)
    elif max <= np.iinfo(np.uint32).max:
        return np.asarray(spectrum, dtype=np.uint32)
    elif max <= np.iinfo(np.uint64).max:
        return np.asarray(spectrum, dtype=np.uint64)

    raise NotImplementedError(
        "Greater than 64 bit uint values are not supported"
    )
