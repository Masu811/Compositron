from time import perf_counter_ns
import numpy as np

from compositron.cdbs.cdbspectrum import Axis, DiagonalArea, LinearProjectionBins, LineshapeParamDefinition
from compositron.core.measurement import Measurement
from compositron.core.utils import EcalCorrectionOrder, KeV


def print_mean_and_std(name: str, data: np.ndarray) -> None:
    mean = np.mean(data)
    std = np.std(data)

    print(f"{name}: {mean} ns +/- {std} ns")


def main():
    n_iter = 1000

    imports = np.empty(n_iter)
    corr = np.empty(n_iter)
    s_calc = np.empty(n_iter)
    proj_axes = np.empty(n_iter)
    proj_diag = np.empty(n_iter)

    for i in range(n_iter):
        start = perf_counter_ns()
        m = Measurement(
            "../../../../testdata/depth-profile_Copper_0000.n42"
        )
        stop = perf_counter_ns()
        imports[i] = stop - start

        c = m.cdbs["OAA x OAB"]

        start = perf_counter_ns()
        c.correct_ecal(EcalCorrectionOrder.FIRST)
        stop = perf_counter_ns()
        corr[i] = stop - start

        param = LineshapeParamDefinition(
            name = "S",
            num = [DiagonalArea(KeV(2), KeV(1))],
            denom = [DiagonalArea(KeV(10), KeV(1))],
            in_corrected_peak = False,
        )

        start = perf_counter_ns()
        c.calc_lineshape_param(param)
        stop = perf_counter_ns()
        s_calc[i] = stop - start

        start = perf_counter_ns()
        c.project_digonal(
            LinearProjectionBins(KeV(0.1)),
            KeV(2.),
            KeV(100.),
            False
        )
        stop = perf_counter_ns()
        proj_diag[i] = stop - start

        start = perf_counter_ns()
        c.project_axes(Axis.FIRST_DET)
        stop = perf_counter_ns()
        proj_axes[i] = stop - start

    print_mean_and_std("Import", imports)
    print_mean_and_std("Ecal Corr", corr)
    print_mean_and_std("S Param", s_calc)
    print_mean_and_std("Proj Diag", proj_diag)
    print_mean_and_std("Proj Axes", proj_axes)


if __name__ == "__main__":
    main()
