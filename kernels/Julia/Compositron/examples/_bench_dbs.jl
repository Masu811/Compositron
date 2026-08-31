using Compositron.Core
using Compositron.Utils
using Compositron.DBS

using Statistics


function work()
    import_bench = @timed Measurement(
        "../../../../testdata/depth-profile_Copper_0000.n42"
    )
    m = import_bench.value
    import_time = import_bench.time * 1e6

    d = m.dbs["OAB"]

    corr_bench = @timed correct_ecal!(d, EcalCorrOrder_First)
    corr_time = corr_bench.time * 1e6

    bg_sub_bench = @timed extract_peak!(d, 30., PeakModel_ErfLinear2Gauss)
    bg_sub_time = bg_sub_bench.time * 1e6

    param = LineshapeParamDefinition(
        "S",
        [PeakArea(KeV(2))],
        [PeakArea(KeV(20))],
        true,
    )

    s_calc_bench = @timed calc_lineshape_param!(d, param, BinIntScheme_Const)
    s_calc_time = s_calc_bench.time * 1e6

    (import_time, corr_time, bg_sub_time, s_calc_time)
end


function print_mean_and_std(name::String, data::AbstractVector)
    m = mean(data)
    s = std(data)

    println("$name: $m µs +/- $s µs")
end


n_iter = 10000

imports = zeros(Float64, n_iter)
corr = zeros(Float64, n_iter)
bg_sub = zeros(Float64, n_iter)
s_calc = zeros(Float64, n_iter)

work()

for i in 1:n_iter
    (import_time, corr_time, bg_sub_time, s_calc_time) = work()

    imports[i] = import_time
    corr[i] = corr_time
    bg_sub[i] = bg_sub_time
    s_calc[i] = s_calc_time
end

print_mean_and_std("Import", imports)
print_mean_and_std("Ecal Corr", corr)
print_mean_and_std("BG Sub", bg_sub)
print_mean_and_std("S Param", s_calc)
