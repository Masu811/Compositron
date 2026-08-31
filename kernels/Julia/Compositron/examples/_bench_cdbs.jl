using Compositron.Core
using Compositron.Utils
using Compositron.CDBS

using Statistics


function work()
    import_bench = @timed Measurement(
        "../../../../testdata/depth-profile_Copper_0000.n42"
    )
    m = import_bench.value
    import_time = import_bench.time * 1e6

    c = m.cdbs["OAA x OAB"]

    corr_bench = @timed correct_ecal!(c, EcalCorrOrder_First)
    corr_time = corr_bench.time * 1e6

    param = LineshapeParamDefinition(
        "S",
        [DiagonalArea(KeV(2), KeV(1))],
        [DiagonalArea(KeV(10), KeV(1))],
        false
    )

    s_calc_bench = @timed calc_lineshape_param!(c, param)
    s_calc_time = s_calc_bench.time * 1e6


    proj_diag_bench = @timed project_diagonal(
        c, LinearProjectionBins(KeV(0.1)), KeV(2), KeV(100), false
    )
    proj_diag_time = proj_diag_bench.time * 1e6

    proj_axes_bench = @timed project_axes(c, Axis_FirstDet)
    proj_axes_time = proj_axes_bench.time * 1e6

    (import_time, corr_time, s_calc_time, proj_diag_time, proj_axes_time)
end


function print_mean_and_std(name::String, data::AbstractVector)
    m = mean(data)
    s = std(data)

    println("$name: $m µs +/- $s µs")
end


n_iter = 10000

imports = zeros(Float64, n_iter)
corr = zeros(Float64, n_iter)
s_calc = zeros(Float64, n_iter)
proj_diag = zeros(Float64, n_iter)
proj_axes = zeros(Float64, n_iter)

work()

for i in 1:n_iter
    (import_time, corr_time, s_calc_time, proj_diag_time, proj_axes_time) = work()

    imports[i] = import_time
    corr[i] = corr_time
    s_calc[i] = s_calc_time
    proj_diag[i] = proj_diag_time
    proj_axes[i] = proj_axes_time
end

print_mean_and_std("Import", imports)
print_mean_and_std("Ecal Corr", corr)
print_mean_and_std("S Param", s_calc)
print_mean_and_std("Proj Diag", proj_diag)
print_mean_and_std("Proj Axes", proj_axes)
