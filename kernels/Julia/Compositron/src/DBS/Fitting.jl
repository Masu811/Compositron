using LsqFit
using SpecialFunctions: erf

using ..Utils: SimpleFitParam
using ..Constants: TWO_OVER_SQRT_PI

function gauss(x, amp, x0, sig)
    amp * exp(-0.5 * ((x - x0) / sig)^2)
end

function fit_gaussian(x::Vector, y::Vector, init::Vector)
    m(t, p) = gauss.(x, p[1], p[2], p[3])

    function j_m(t, p)
        J = Array{Float64}(undef, length(t), length(p))
        e = gauss.(x, p[1], p[2], p[3])
        u = @. (t - p[2]) / p[3]
        z = @. u / p[3]
        @. J[:,1] = e / p[1]
        @. J[:,2] = e * z
        @. J[:,3] = e * u * z
        J
    end

    w = @. 1 / max(y, 1)

    fit = curve_fit(m, j_m, x, y, w, init)

    opt = coef(fit)

    err = try
        stderror(fit)
    catch
        fill(NaN, length(opt))
    end

    Dict(
        "amp_1" => SimpleFitParam(opt[1], err[1]),
        "x0_1" => SimpleFitParam(opt[2], err[2]),
        "sig_1" => SimpleFitParam(opt[3], err[3]),
    )
end

function erf_linear_background(x, amp, x0, sig, lin, off)
    amp * erf((x - x0) / sig) + lin * x + off
end

function fit_erf_linear_1_gauss(x::Vector, y::Vector, init::Vector)
    function m(t, p)
        u = @. (t - p[2]) / p[3]

        @. p[1] * exp(-0.5 * u^2) + p[4] * erf(u) + p[5] * t + p[6]
    end

    function j_m(t, p)
        J = Array{Float64}(undef, length(t), length(p))
        e = gauss.(x, p[1], p[2], p[3])
        u = @. (t - p[2]) / p[3]
        z = @. u / p[3]
        @. J[:,1] = e / p[1]
        @. J[:,2] = e * z - TWO_OVER_SQRT_PI * exp(-u^2) / p[3]
        @. J[:,3] = e * u * z + TWO_OVER_SQRT_PI * exp(-u^2) * u / p[3]
        @. J[:,4] = erf(u)
        @. J[:,5] = 1
        @. J[:,6] = 0
        J
    end

    w = @. 1 / max(y, 1)

    fit = curve_fit(m, j_m, x, y, w, init)

    opt = coef(fit)

    err = try
        stderror(fit)
    catch
        fill(NaN, length(opt))
    end

    Dict(
        "amp_1" => SimpleFitParam(opt[1], err[1]),
        "x0_1" => SimpleFitParam(opt[2], err[2]),
        "sig_1" => SimpleFitParam(opt[3], err[3]),
        "erf_amp" => SimpleFitParam(opt[4], err[4]),
        "lin" => SimpleFitParam(opt[5], err[5]),
        "const" => SimpleFitParam(opt[6], err[6]),
    )
end
