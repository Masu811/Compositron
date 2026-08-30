using Printf

using Compositron.Core
using Compositron.DBS


m = Measurement(
    "../../../../testdata/depth-profile_Copper_0000.n42"
)

println(shape(m))

d = m.dbs["OAB"]

default_analyze!(d)

s = d.lineshape_params["S"]
w = d.lineshape_params["W"]
v = d.lineshape_params["V/P"]
p = d.lineshape_params["P/T"]

@printf("  S   = %.6f +/- %.6f\n", s.val, s.err)
@printf("  W   = %.6f +/- %.6f\n", w.val, w.err)
@printf("  V/P = %.6f +/- %.6f\n", v.val, v.err)
@printf("  P/T = %.6f +/- %.6f\n", p.val, p.err)
