using DelimitedFiles

using Compositron.Utils: first_order, KeV
using Compositron.Core: Measurement
using Compositron.CDBS: correct_ecal!, LineshapeParamDefinition, DiagonalArea, calc_lineshape_param!, project_diagonal, LinearProjectionBins

m = Measurement(
    "../../../../testdata/depth-profile_Copper_0000.n42"
)

c = m.cdbs["OAA x OAB"]

correct_ecal!(c, first_order)

# S parameter calculation

ls_param = LineshapeParamDefinition(
    "S",
    [DiagonalArea(KeV(2), KeV(1), KeV(0), KeV(0))],
    [DiagonalArea(KeV(Inf), KeV(1), KeV(0), KeV(0))],
    false
)

println("S calculation duration (cold-start): ")

@time calc_lineshape_param!(c, ls_param)

println("S calculation duration (warmed-up): ")

@time calc_lineshape_param!(c, ls_param)

s = c.lineshape_params["S"]

println("Calculated S parameter: $(s.val) +/- $(s.err)")

# Projection onto diagonal

println("Projection duration (cold-start): ")

@time project_diagonal(c, LinearProjectionBins(KeV(0.1)), KeV(2), KeV(10), false)

println("Projection duration (warmed-up): ")

p = @time project_diagonal(c, LinearProjectionBins(KeV(0.1)), KeV(2), KeV(10), false)

# Write projection to file for inspection

writedlm("projection.csv", p.spectrum, ",")
