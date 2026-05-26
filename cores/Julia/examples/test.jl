include("../src/stacs.jl")

using .STACS

m = STACS.DopplerMeasurement("../../../testdata/depth-profile_Copper_0000.n42")

println(STACS.Measurement.shape(m))
