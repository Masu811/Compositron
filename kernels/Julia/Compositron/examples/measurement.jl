include("../src/Compositron.jl")

import .Compositron: Measurement
import .Compositron.Core: shape

m = Measurement(
    "../../../../testdata/depth-profile_Copper_0000.n42"
)

println(shape(m))
