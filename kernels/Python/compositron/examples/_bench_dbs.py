from time import perf_counter_ns
import numpy as np

from compositron.core.measurement import Measurement
from compositron.core.utils import EcalCorrectionOrder, KeV
from compositron.dbs.dbspectrum import BinIntegrationScheme, LineshapeParamDefinition, PeakArea, PeakModel


def print_mean_and_std(name: str, data: np.ndarray) -> None:
    mean = np.mean(data)
    std = np.std(data)

    print(f"{name}: {mean} ns +/- {std} ns")


def main():
    n_iter = 10000

    imports = np.empty(n_iter)
    corr = np.empty(n_iter)
    bg_sub = np.empty(n_iter)
    s_calc = np.empty(n_iter)

    for i in range(n_iter):
        start = perf_counter_ns()
        m = Measurement(
            "../../../../testdata/depth-profile_Copper_0000.n42"
        )
        stop = perf_counter_ns()
        imports[i] = stop - start

        d = m.dbs["OAB"]

        start = perf_counter_ns()
        d.correct_ecal(EcalCorrectionOrder.FIRST)
        stop = perf_counter_ns()
        corr[i] = stop - start

        start = perf_counter_ns()
        d.extract_peak(30., PeakModel.ERFLINEAR2GAUSS)
        stop = perf_counter_ns()
        bg_sub[i] = stop - start

        param = LineshapeParamDefinition(
            "S",
            [PeakArea(KeV(2.))],
            [PeakArea(KeV(20.))],
            True
        )

        start = perf_counter_ns()
        d.calc_lineshape_param(param, BinIntegrationScheme.CONST)
        stop = perf_counter_ns()
        s_calc[i] = stop - start

    print_mean_and_std("Import", imports)
    print_mean_and_std("Ecal Corr", corr)
    print_mean_and_std("BG Sub", bg_sub)
    print_mean_and_std("S Param", s_calc)


if __name__ == "__main__":
    main()
