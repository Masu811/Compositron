using Compositron.Core: Measurement, shape

m = Measurement(
    "../../../../testdata/depth-profile_Copper_0000.n42"
)

println(shape(m))
