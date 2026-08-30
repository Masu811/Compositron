import time

import numpy as np

from compositron.cdbs.cdbspectrum import DiagonalArea, LinearProjectionBins, LineshapeParamDefinition
from compositron.core.measurement import Measurement
from compositron.core.utils import EcalCorrectionOrder, KeV

# Import

m = Measurement(
    "../../../../testdata/depth-profile_Copper_0000.n42"
)

c = m.cdbs["OAA x OAB"]

c.correct_ecal(EcalCorrectionOrder.FIRST)

# S parameter calculation

ls_param = LineshapeParamDefinition(
    name = "S",
    num = [DiagonalArea(KeV(2), KeV(1))],
    denom = [DiagonalArea(KeV(10), KeV(1))],
    in_corrected_peak = False,
)

start = time.time()
s = c.calc_lineshape_param(ls_param)
stop = time.time()

print(f"Calculated S parameter: {s.val} +/- {s.err}");
print(f"Calculation duration: {(stop - start) * 1e6:.1f} µs");

# Projection onto diagonal

start = time.time()
p = c.project_digonal(
    LinearProjectionBins(KeV(0.1)),
    KeV(2.),
    KeV(10.),
    False
)
stop = time.time()

print(f"Projection duration: {(stop - start) * 1e6:.1f} µs");

# Write projection to file for inspection

with open("projection.csv", "w") as f:
    for elem in p.spectrum:
        f.write(f"{elem}\n")
