using LsqFit
using SpecialFunctions: erf

using ..Utils: SimpleFitParam
using ..Constants: TWO_OVER_SQRT_PI

function gauss(x, amp, x0, sig)
    amp * exp(-0.5 * ((x - x0) / sig)^2)
end

function fit_gaussian(x::Vector, y::Vector, init::Vector)
    m(t, p) = gauss.(t, p[1], p[2], p[3])

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
        amp_gauss_1 = p[1];
        x0_1 = p[2];
        sig_1 = p[3];
        amp_erf = p[4];
        lin = p[5];
        off = p[5];

        u_1 = @. (t - x0_1) / sig_1;

        gauss_1 = @. amp_gauss_1 * exp(-0.5 * u_1 * u_1);

        er = @. amp_erf * erf(u_1);

        @. gauss_1 + er + lin * x + off;
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

function fit_erf_linear_2_gauss(x::Vector, y::Vector, init::Vector)
    function m(t, p)
        amp_gauss_1 = p[1];
        x0_1 = p[2];
        sig_1 = p[3];
        amp_gauss_2 = p[4];
        x0_2 = p[5];
        sig_2 = p[6];
        amp_erf = p[7];
        lin = p[8];
        off = p[9];

        u_1 = @. (t - x0_1) / sig_1;
        u_2 = @. (t - x0_2) / sig_2;

        gauss_1 = @. amp_gauss_1 * exp(-0.5 * u_1 * u_1);
        gauss_2 = @. amp_gauss_2 * exp(-0.5 * u_2 * u_2);

        er = @. amp_erf * erf(u_1);

        @. gauss_1 + gauss_2 + erf + lin * x + off;
    end

    function j_m(t, p)
        J = Array{Float64}(undef, length(t), length(p))
        e_1 = gauss.(x, p[1], p[2], p[3])
        e_2 = gauss.(x, p[4], p[5], p[6])
        u_1 = @. (t - p[2]) / p[3]
        u_2 = @. (t - p[5]) / p[6]
        z_1 = @. u / p[3]
        z_2 = @. u / p[6]
        @. J[:,1] = e_1 / p[1]
        @. J[:,2] = e_1 * z_1 - TWO_OVER_SQRT_PI * exp(-u_1^2) / p[3]
        @. J[:,3] = e_1 * u_1 * z_1 + TWO_OVER_SQRT_PI * exp(-u_1^2) * u_1 / p[3]
        @. J[:,4] = e_2 / p[4]
        @. J[:,5] = e_2 * z_2
        @. J[:,6] = e_2 * u_2 * z_2
        @. J[:,7] = erf(u_1)
        @. J[:,8] = 1
        @. J[:,9] = 0
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
        "amp_2" => SimpleFitParam(opt[4], err[4]),
        "x0_2" => SimpleFitParam(opt[5], err[5]),
        "sig_2" => SimpleFitParam(opt[6], err[6]),
        "erf_amp" => SimpleFitParam(opt[7], err[7]),
        "lin" => SimpleFitParam(opt[8], err[8]),
        "const" => SimpleFitParam(opt[9], err[9]),
    )
end

function fit_erf_linear_3_gauss(x::Vector, y::Vector, init::Vector)
    function m(t, p)
        amp_gauss_1 = p[1];
        x0_1 = p[2];
        sig_1 = p[3];
        amp_gauss_2 = p[4];
        x0_2 = p[5];
        sig_2 = p[6];
        amp_gauss_3 = p[7];
        x0_3 = p[8];
        sig_3 = p[9];
        amp_erf = p[10];
        lin = p[11];
        off = p[12];

        u_1 = @. (t - x0_1) / sig_1;
        u_2 = @. (t - x0_2) / sig_2;
        u_3 = @. (t - x0_3) / sig_3;

        gauss_1 = @. amp_gauss_1 * exp(-0.5 * u_1 * u_1);
        gauss_2 = @. amp_gauss_2 * exp(-0.5 * u_2 * u_2);
        gauss_3 = @. amp_gauss_3 * exp(-0.5 * u_3 * u_3);

        er = @. amp_erf * erf(u_1);

        @. gauss_1 + gauss_2 + gauss_3 + erf + lin * x + off;
    end

    function j_m(t, p)
        J = Array{Float64}(undef, length(t), length(p))
        e_1 = gauss.(x, p[1], p[2], p[3])
        e_2 = gauss.(x, p[4], p[5], p[6])
        e_3 = gauss.(x, p[7], p[8], p[9])
        u_1 = @. (t - p[2]) / p[3]
        u_2 = @. (t - p[5]) / p[6]
        u_3 = @. (t - p[8]) / p[9]
        z_1 = @. u / p[3]
        z_2 = @. u / p[6]
        z_3 = @. u / p[9]
        @. J[:,1] = e_1 / p[1]
        @. J[:,2] = e_1 * z_1 - TWO_OVER_SQRT_PI * exp(-u_1^2) / p[3]
        @. J[:,3] = e_1 * u_1 * z_1 + TWO_OVER_SQRT_PI * exp(-u_1^2) * u_1 / p[3]
        @. J[:,4] = e_2 / p[4]
        @. J[:,5] = e_2 * z_2
        @. J[:,6] = e_2 * u_2 * z_2
        @. J[:,7] = e_3 / p[7]
        @. J[:,8] = e_3 * z_3
        @. J[:,9] = e_3 * u_3 * z_3
        @. J[:,10] = erf(u_1)
        @. J[:,11] = 1
        @. J[:,12] = 0
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
        "amp_2" => SimpleFitParam(opt[4], err[4]),
        "x0_2" => SimpleFitParam(opt[5], err[5]),
        "sig_2" => SimpleFitParam(opt[6], err[6]),
        "amp_3" => SimpleFitParam(opt[7], err[7]),
        "x0_3" => SimpleFitParam(opt[8], err[8]),
        "sig_3" => SimpleFitParam(opt[9], err[9]),
        "erf_amp" => SimpleFitParam(opt[10], err[10]),
        "lin" => SimpleFitParam(opt[11], err[11]),
        "const" => SimpleFitParam(opt[12], err[12]),
    )
end
