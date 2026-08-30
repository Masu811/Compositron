using NonlinearSolve
using SpecialFunctions: erfc

using ..Utils
using ..Constants: TWO_OVER_SQRT_PI

export gauss, fit_gauss, erf_linear_background, fit_erf_linear_1_gauss,
    fit_erf_linear_2_gauss, fit_erf_linear_3_gauss


function gauss(x, amp, x0, sig)
    amp * exp(-0.5 * ((x - x0) / sig)^2)
end


function fit_gauss(x::Vector, y::Vector, init::Vector)
    m(u, p) = gauss.(p, u[1], u[2], u[3]) .- y

    prob = NonlinearLeastSquaresProblem(m, init, x)

    sol = solve(prob, TrustRegion(; autodiff = AutoForwardDiff()))

    opt = sol.u

    err = try
        Q, R = qr(J)
        Rinv = inv(R)

        cov = red_chi_2 .* (Rinv * Rinv')

        fit_result.cov = cov

        sqrt.(diag(cov))
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
    amp * erfc((x - x0) / sig) + lin * x + off
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

        erf = @. amp_erf * erfc(u_1);

        @. gauss_1 + erf + lin * t + off;
    end

    w = @. 1 / sqrt(max(y, 1))

    prob = NonlinearLeastSquaresProblem((u, p) -> w .* (m(p, u) .- y), init, x)

    sol = solve(prob, TrustRegion(; autodiff = AutoForwardDiff()))

    opt = sol.u

    err = try
        Q, R = qr(J)
        Rinv = inv(R)

        cov = red_chi_2 .* (Rinv * Rinv')

        fit_result.cov = cov

        sqrt.(diag(cov))
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

        erf = @. amp_erf * erfc(u_1);

        @. gauss_1 + gauss_2 + erf + lin * x + off;
    end

    w = @. 1 / sqrt(max(y, 1))

    prob = NonlinearLeastSquaresProblem((u, p) -> w .* (m(p, u) .- y), init, x)

    sol = solve(prob, TrustRegion(; autodiff = AutoForwardDiff()))

    opt = sol.u

    err = try
        Q, R = qr(J)
        Rinv = inv(R)

        cov = red_chi_2 .* (Rinv * Rinv')

        fit_result.cov = cov

        sqrt.(diag(cov))
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

        erf = @. amp_erf * erfc(u_1);

        @. gauss_1 + gauss_2 + gauss_3 + erf + lin * x + off;
    end

    w = @. 1 / sqrt(max(y, 1))

    prob = NonlinearLeastSquaresProblem((u, p) -> w .* (m(p, u) .- y), init, x)

    sol = solve(prob, TrustRegion(; autodiff = AutoForwardDiff()))

    opt = sol.u

    err = try
        Q, R = qr(J)
        Rinv = inv(R)

        cov = red_chi_2 .* (Rinv * Rinv')

        fit_result.cov = cov

        sqrt.(diag(cov))
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
