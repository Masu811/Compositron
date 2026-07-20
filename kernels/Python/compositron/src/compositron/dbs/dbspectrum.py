import numpy as np
import numpy.typing as npt

from ..core.utils import to_min_type_spectrum


Ecal = tuple[np.float64, np.float64]

NAN = np.float64(np.nan)


class DBSpectrum:
    def __init__(
        self,
        spectrum: npt.ArrayLike,
        detname: str,
        ecal: Ecal,
    ):
        self.spectrum = to_min_type_spectrum(spectrum)
        self.detname = detname
        self.ecal = ecal
        self.s = NAN
        self.ds = NAN
        self.w = NAN
        self.ds = NAN
        self.v2p = NAN
        self.dv2p = NAN
        self.counts = np.sum(self.spectrum, dtype=np.uint64)
        self.dcounts = np.sqrt(self.counts)
        self.peak_counts = NAN
        self.dpeak_counts = NAN
