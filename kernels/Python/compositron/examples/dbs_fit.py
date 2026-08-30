from compositron.core.measurement import Measurement


m = Measurement(
    "../../../../testdata/depth-profile_Copper_0000.n42"
)

print(m.shape)

d = m.dbs["OAB"]

d.default_analyze()

s = d.lineshape_params["S"]
w = d.lineshape_params["W"]
v = d.lineshape_params["V/P"]
p = d.lineshape_params["P/T"]

print(f"  S     = {s.val:.6f} +/- {s.err:.6f}")
print(f"  W     = {w.val:.6f} +/- {w.err:.6f}")
print(f"  V/P   = {v.val:.6f} +/- {v.err:.6f}")
print(f"  P/T   = {p.val:.6f} +/- {p.err:.6f}")
