using Compositron.Core: Measurement, shape
using Compositron.DBS: default_analyze!

m = Measurement(
    "../../../../testdata/depth-profile_Copper_0000.n42"
)

println(shape(m))

default_analyze!(m.dbs["OAB"])

display(m.dbs["OAB"].lineshape_params)
