using Compositron.Core: Measurement, shape

m = Measurement(
    "../../../../testdata/1_W3Re6_2.00keV_301.0K_10.0Mio_PALS_01.05.26_08.42.22.dat",
)

println(shape(m))
