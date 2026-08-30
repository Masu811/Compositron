using DelimitedFiles

using Compositron.Utils
using Compositron.Core
using Compositron.CDBS


m = Measurement(
    "../../../../testdata/depth-profile_Copper_0000.n42"
)

c = m.cdbs["OAA x OAB"]

correct_ecal!(c, EcalCorrOrder_First)

# S parameter calculation

ls_param = LineshapeParamDefinition(
    "S",
    [DiagonalArea(KeV(2), KeV(1))],
    [DiagonalArea(KeV(10), KeV(1))],
    false
)

calc_lineshape_param!(c, ls_param)

println("S calculation duration")

@time calc_lineshape_param!(c, ls_param)

s = c.lineshape_params["S"]

println("Calculated S parameter: $(s.val) +/- $(s.err)")

# Projection onto diagonal

project_diagonal(c, LinearProjectionBins(KeV(0.1)), KeV(2), KeV(10), false)

println("Projection duration: ")

p = @time project_diagonal(c, LinearProjectionBins(KeV(0.1)), KeV(2), KeV(10), false)

# Write projection to file for inspection

writedlm("projection.csv", p.spectrum, ",")
