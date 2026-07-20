import numpy as np
import numpy.typing as npt

from ..core.utils import to_min_type_spectrum


Ecal = tuple[tuple[np.float64, np.float64], tuple[np.float64, np.float64]]

NAN = np.float64(np.nan)


class CDBSpectrum:
    def __init__(
        self,
        spectrum: npt.ArrayLike,
        detpair: str,
        ecal: Ecal,
    ):
        self.spectrum = to_min_type_spectrum(spectrum)
        self.detpair = detpair
        self.ecal = ecal
        self.s = NAN
        self.ds = NAN
        self.w = NAN
        self.ds = NAN
        self.counts = np.sum(self.spectrum, dtype=np.uint64)
        self.dcounts = np.sqrt(self.counts)
        self.roi_counts = NAN
        self.droi_counts = NAN
