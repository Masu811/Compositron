using SpecialFunctions: erfc
using NonlinearSolve
using LinearAlgebra

using ..Constants: SQRT_2, SQRT_2PI
using ..Utils: FitParam

function lifetime_spectrum_component(
    t::Vector{Float64},
    n,
    t0,
    tau,
    l_int,
    sig,
    r_int,
    r_t0,
)
    # TODO: rewrite with unicode symbols

    sig_over_tau_sqrt_2 = sig / (tau * SQRT_2)
    sig_over_tau_sqrt_2_sq = sig_over_tau_sqrt_2 ^ 2
    one_over_tau = 1 / tau
    one_over_sig_sqrt_2 = 1 / (sig * SQRT_2)

    t_shifted = t .- (t0 + r_t0)
    e = @. exp(-t_shifted * one_over_tau + sig_over_tau_sqrt_2_sq)
    c = @. erfc(sig_over_tau_sqrt_2 - t_shifted * one_over_sig_sqrt_2)

    @. 0.5 * n * l_int * r_int * one_over_tau * e * c
end

mutable struct FitContext
    n_l::Int
    n_r::Int
    signed_param_idcs::Vector{Int}  # pos: idx in p, neg: idx in fixed_params
    fixed_params::Vector{Float64}
    available_l_intensity::Float64
    available_r_intensity::Float64
    last_l_vary_idx::Int
    last_r_vary_idx::Int
    lower_bounds::Vector{Float64}
    upper_bounds::Vector{Float64}
    init::Vector{Float64}
    varied_param_idx::Int
    fixed_param_idx::Int
end

function _get_param(idx::Int, p, ctx::FitContext)
    signed_idx = ctx.signed_param_idcs[idx]

    if signed_idx > 0
        p[signed_idx]
    else
        ctx.fixed_params[abs(signed_idx)]
    end
end

function _transform_intensities!(p, ctx::FitContext)
    leftover_intensity = ctx.available_l_intensity
    abs_l_intensities = zeros(eltype(p), ctx.n_l)

    for i in 1:ctx.n_l
        signed_idx = ctx.signed_param_idcs[4 + 2 * (i - 1) + 1]

        if signed_idx > 0
            next_intensity = p[signed_idx] * leftover_intensity
            abs_l_intensities[i] = next_intensity
            leftover_intensity = max(leftover_intensity - next_intensity, 0)
        elseif i == ctx.last_l_vary_idx
            abs_l_intensities[i] = leftover_intensity
        else
            abs_l_intensities[i] = ctx.fixed_params[abs(signed_idx)]
        end
    end

    leftover_intensity = ctx.available_r_intensity
    abs_r_intensities = zeros(eltype(p), ctx.n_r)

    for i in 1:ctx.n_r
        signed_idx = ctx.signed_param_idcs[4 + 2 * ctx.n_l + 3 * (i - 1) + 1]

        if signed_idx > 0
            next_intensity = p[signed_idx] * leftover_intensity
            abs_r_intensities[i] = next_intensity
            leftover_intensity = max(leftover_intensity - next_intensity, 0)
        elseif i == ctx.last_r_vary_idx
            abs_r_intensities[i] = leftover_intensity
        else
            abs_r_intensities[i] = ctx.fixed_params[abs(signed_idx)]
        end
    end

    return abs_l_intensities, abs_r_intensities
end

function _get_all_params!(p, ctx::FitContext)
    abs_l_intensities, abs_r_intensities = _transform_intensities!(p, ctx)

    p_buffer = zeros(eltype(p), 3 + 2 * ctx.n_l, 3 * ctx.n_r)

    for i in 1:3
        p_buffer[i] = _get_param(i, p, ctx)
    end

    for i in 1:ctx.n_l
        base_idx = 4 + 2 * (i - 1)
        p_buffer[base_idx] = _get_param(base_idx, p, ctx)
        p_buffer[base_idx + 1] = abs_l_intensities[i]
    end

    for i in 1:ctx.n_r
        base_idx = 4 + 2 * ctx.n_l + 3 * (i - 1)
        p_buffer[base_idx] = _get_param(base_idx, p, ctx)
        p_buffer[base_idx + 1] = abs_r_intensities[i]
        p_buffer[base_idx + 2] = _get_param(base_idx + 2, p, ctx)
    end

    p_buffer
end

function _m(t, p, ctx::FitContext)
    params = _get_all_params!(p, ctx)
    n_l = ctx.n_l
    n_r = ctx.n_r

    n = params[1]
    bg = params[2]
    t0 = params[3]

    y = fill(bg, length(t))

    for i in 1:n_l
        l_base_idx = 4 + 2 * (i - 1)
        tau = params[l_base_idx]
        l_int = params[l_base_idx + 1]

        for j in 1:n_r
            r_base_idx = 4 + 2 * n_l + 3 * (j - 1)
            sig = params[r_base_idx]
            r_int = params[r_base_idx + 1]
            r_t0 = params[r_base_idx + 2]

            y .+= lifetime_spectrum_component(
                t, n, t0, tau, l_int, sig, r_int, r_t0
            )
        end
    end

    y
end

function _j_m(t, p, ctx::FitContext)
    params = _get_all_params!(p, ctx)

    n_l = ctx.n_l
    n_r = ctx.n_r

    n = params[1]
    t0 = params[3]

    n_dpoints = length(t)

    y = fill(0., n_dpoints)

    terms = Matrix{Float64}(undef, n_dpoints, n_l * n_r)
    deriv = [Matrix{Float64}(undef, n_dpoints, 4) for _ in  1:(n_l * n_r)]

    for i in 1:n_l
        l_base_idx = 4 + 2 * (i - 1)
        tau = params[l_base_idx]
        l_int = params[l_base_idx + 1]

        for j in 1:n_r
            col_idx = j + (i - 1) * n_r
            r_base_idx = 4 + 2 * n_l + 3 * (j - 1)

            sig = params[r_base_idx]
            r_int = params[r_base_idx + 1]
            r_t0 = params[r_base_idx + 2]

            sig_over_tau_sqrt_2 = sig / (tau * SQRT_2)
            sig_over_tau_sqrt_2_sq = sig_over_tau_sqrt_2 ^ 2
            one_over_tau = 1 / tau
            one_over_tau_sq = one_over_tau^2
            one_over_sig_sqrt_2 = 1 / (sig * SQRT_2)
            one_over_tau_sqrt_2_pi = 1 / (tau * SQRT_2PI)

            t_shifted = t .- (t0 + r_t0)

            a = @. -t_shifted  * one_over_tau + sig_over_tau_sqrt_2_sq
            b = @. sig_over_tau_sqrt_2 - t_shifted * one_over_sig_sqrt_2

            decay = @. n * 0.5 * l_int * r_int * one_over_tau * exp(a) * erfc(b)

            exp_comb = @. n * l_int * r_int * one_over_tau_sqrt_2_pi * exp(a - b^2)

            @. terms[:, col_idx] = decay

            @. y += decay

            df_dt0 = @. -exp_comb / sig + decay / tau

            @. deriv[col_idx][:, 1] = df_dt0
            @. deriv[col_idx][:, 2] = (
                exp_comb * sig * one_over_tau_sq
                + decay * (
                    t_shifted * one_over_tau_sq - sig^2 / tau^3 - one_over_tau
                )
            )
            @. deriv[col_idx][:, 3] = (
                decay * sig * one_over_tau_sq
                - exp_comb * (t_shifted / sig^2 + one_over_tau)
            )
            @. deriv[col_idx][:, 4] = df_dt0
        end
    end

    df_dl_int = Matrix{Float64}(undef, n_dpoints, n_l)

    for i in 1:n_l
        signed_idx = ctx.signed_param_idcs[4 + 2 * (i - 1) + 1]
        if signed_idx > 0
            tmp = fill(0., n_dpoints)

            for n in i:n_l
                for m in 1:n_r
                    a = if n == i
                        1 / (p[signed_idx] + 1e-5)
                    else
                        -1 / (1 - p[signed_idx] + 1e-5)
                    end

                    @. tmp += a * terms[:, m + n_r * (n - 1)]
                end
            end

            @. df_dl_int[:, i] = tmp
        else
            @. df_dl_int[:, i] = 0.
        end
    end

    df_dr_int = Matrix{Float64}(undef, n_dpoints, n_r)

    for i in 1:n_r
        signed_idx = ctx.signed_param_idcs[4 + 2 * n_l + 3 * (i - 1) + 1]
        if signed_idx > 0
            tmp = fill(0., n_dpoints)

            for m in i:n_r
                for n in 1:n_l
                    a = if m == i
                        1 / (p[signed_idx] + 1e-5)
                    else
                        -1 / (1 - p[signed_idx] + 1e-5)
                    end

                    @. tmp += a * terms[:, m + n_r * (n - 1)]
                end
            end

            @. df_dr_int[:, i] = tmp
        else
            @. df_dr_int[:, i] = 0.
        end
    end

    j = Matrix{Float64}(undef, n_dpoints, length(p))
    j_idx = 1

    if ctx.signed_param_idcs[1] > 0
        @. j[:, j_idx] = y / n
        j_idx += 1
    end

    if ctx.signed_param_idcs[2] > 0
        @. j[:, j_idx] = 1.
        j_idx += 1
    end

    if ctx.signed_param_idcs[3] > 0
        j[:, j_idx] .= sum(mat[:, 1] for mat in deriv)
        j_idx += 1
    end

    for i in 1:n_l
        flat_idx = 4 + 2 * (i - 1)

        if ctx.signed_param_idcs[flat_idx] > 0
            j[:, j_idx] .= sum(
                deriv[k][:, 2] for k in ((i - 1) * n_r + 1):(i * n_r)
            )
            j_idx += 1
        end

        if ctx.signed_param_idcs[flat_idx + 1] > 0
            @. j[:, j_idx] = df_dl_int[:, i]
            j_idx += 1
        end
    end

    for i in 1:n_r
        flat_idx = 4 + 2 * n_l + 3 * (i - 1)

        if ctx.signed_param_idcs[flat_idx] > 0
            j[:, j_idx] .= sum(
                deriv[k][:, 3] for k in i:n_r:(n_l * n_r)
            )
            j_idx += 1
        end

        if ctx.signed_param_idcs[flat_idx + 1] > 0
            @. j[:, j_idx] = df_dr_int[:, i]
            j_idx += 1
        end

        if ctx.signed_param_idcs[flat_idx + 2] > 0
            j[:, j_idx] .= sum(
                deriv[k][:, 4] for k in i:n_r:(n_l * n_r)
            )
            j_idx += 1
        end
    end

    j
end

function _set_model_defaults!(
    model::LifetimeModel, x::Vector{Float64}, counts::Int, peak_center::Float64
)
    if model.scale_component.default_initialized
        model.scale_component.val = 2. * counts
    end

    if model.shift_component.default_initialized
        l = x[1]
        r = x[end]
        d = r - l
        model.shift_component.min = l - d
        model.shift_component.max = r + d
        model.shift_component.val = peak_center
    end

    for comp in model.lifetime_components
        comp.lifetime.min = max(0, comp.lifetime.min)
        comp.intensity.min = max(0, comp.intensity.min)
        comp.intensity.max = min(1, comp.intensity.max)
    end

    for comp in model.resolution_components
        comp.fwhm.min = max(0, comp.fwhm.min)
        comp.intensity.min = max(0, comp.intensity.min)
        comp.intensity.max = min(1, comp.intensity.max)
    end
end

function _normalize_and_fix_params!(model::LifetimeModel, ctx::FitContext)
    l_norm = sum(comp.intensity.val for comp in model.lifetime_components)

    for comp in model.lifetime_components
        comp.intensity.val /= l_norm
    end

    r_norm = sum(comp.intensity.val for comp in model.resolution_components)

    for comp in model.resolution_components
        comp.intensity.val /= r_norm
    end

    for i in ctx.n_l:-1:1
        param = model.lifetime_components[i].intensity
        if param.vary
            ctx.last_l_vary_idx = i
            param.val = 1
            param.vary = false
            break
        end
    end

    for i in ctx.n_r:-1:1
        param = model.resolution_components[i].intensity
        if param.vary
            ctx.last_r_vary_idx = i
            param.val = 1
            param.vary = false
            break
        end
    end

    if all(comp.t0.vary for comp in model.resolution_components)
        model.resolution_components[1].t0.vary = false
    end
end

function _gobble_component!(
    val::Float64, vary::Bool, min::Float64, max::Float64, ctx::FitContext
)
    if vary
        push!(ctx.init, val)
        push!(ctx.lower_bounds, min)
        push!(ctx.upper_bounds, max)
        push!(ctx.signed_param_idcs, ctx.varied_param_idx)
        ctx.varied_param_idx += 1
    else
        push!(ctx.fixed_params, val)
        push!(ctx.signed_param_idcs, ctx.fixed_param_idx)
        ctx.fixed_param_idx -= 1
    end
end

function _gobble_model!(model::LifetimeModel, ctx::FitContext)
    for comp in [
        model.scale_component, model.background_component, model.shift_component
    ]
        _gobble_component!(comp.val, comp.vary, comp.min, comp.max, ctx)
    end

    remaining_l_int = ctx.available_l_intensity

    for (i, comp) in enumerate(model.lifetime_components)
        lt = comp.lifetime
        _gobble_component!(lt.val, lt.vary, lt.min, lt.max, ctx)

        l_int = comp.intensity
        if l_int.vary
            q = l_int.val / remaining_l_int
            _gobble_component!(q, l_int.vary, l_int.min, l_int.max, ctx)
            remaining_l_int -= l_int.val
        else
            _gobble_component!(l_int.val, l_int.vary, l_int.min, l_int.max, ctx)
        end
    end

    remaining_r_int = ctx.available_r_intensity

    for (i, comp) in enumerate(model.resolution_components)
        fwhm = comp.fwhm
        _gobble_component!(
            fwhm.val / FWHM_OVER_SIGMA,
            fwhm.vary,
            fwhm.min / FWHM_OVER_SIGMA,
            fwhm.max / FWHM_OVER_SIGMA,
            ctx
        )

        r_int = comp.intensity
        if r_int.vary
            q = r_int.val / remaining_r_int
            _gobble_component!(q, r_int.vary, r_int.min, r_int.max, ctx)
            remaining_r_int -= r_int.val
        else
            _gobble_component!(r_int.val, r_int.vary, r_int.min, r_int.max, ctx)
        end

        t0 = comp.t0
        _gobble_component!(t0.val, t0.vary, t0.min, t0.max, ctx)
    end
end

function _prepare_fit!(
    model::LifetimeModel,
    x::Vector{Float64},
    counts::Int,
    peak_center::Float64
)::FitContext
    _set_model_defaults!(model, x, counts, peak_center)

    n_l = length(model.lifetime_components)
    n_r = length(model.resolution_components)
    n_max = 3 + 2 * n_l + 3 * n_r

    ctx = FitContext(
        n_l,
        n_r,
        Vector{Int}(),  # signed_param_idcs
        Vector{Float64}(),  # fixed_params
        1. - sum(
            comp.intensity.val
            for comp in model.lifetime_components if !comp.intensity.vary;
            init=0.
        ),  # available_l_intensity
        1. - sum(
            comp.intensity.val
            for comp in model.resolution_components if !comp.intensity.vary;
            init=0.
        ),  # available_r_intensity
        -1,  # last_l_vary_idx
        -1,  # last_r_vary_idx
        Vector{Float64}(),  # lower_bounds
        Vector{Float64}(),  # upper_bounds
        Vector{Float64}(),  # init
        1,  # varied_param_idx
        -1,  # fixed_param_idx
    )

    _normalize_and_fix_params!(model, ctx)

    _gobble_model!(model, ctx)

    ctx
end

function _paste_fit_result_values!(
    opt::AbstractVector{Float64}, model::LifetimeModel, ctx::FitContext
)
    params = _get_all_params!(opt, ctx)

    model.scale_component.val = params[1]
    model.background_component.val = params[2]
    model.shift_component.val = params[3]

    for i in 1:ctx.n_l
        base_idx = 4 + 2 * (i - 1)
        model.lifetime_components[i].lifetime.val = params[base_idx]
        model.lifetime_components[i].intensity.val = params[base_idx + 1]
    end

    for i in 1:ctx.n_r
        base_idx = 4 + 2 * ctx.n_l + 3 * (i - 1)
        model.resolution_components[i].fwhm.val = params[base_idx] * FWHM_OVER_SIGMA
        model.resolution_components[i].intensity.val = params[base_idx + 1]
        model.resolution_components[i].t0.val = params[base_idx + 2]
    end
end

function _get_error(
    idx::Int, err::AbstractVector{Float64}, ctx::FitContext
)::Float64
    signed_idx = ctx.signed_param_idcs[idx]

    if signed_idx > 0
        err[signed_idx]
    else
        0.
    end
end

function _paste_fit_result_errors!(
    err::AbstractVector{Float64},
    opt::AbstractVector{Float64},
    model::LifetimeModel,
    ctx::FitContext
)
    model.scale_component.err = _get_error(1, err, ctx)
    model.background_component.err = _get_error(2, err, ctx)
    model.shift_component.err = _get_error(3, err, ctx)

    l_cumul_sq_err = 0.

    for i in 1:ctx.n_l
        base_idx = 4 + 2 * (i - 1)
        comp = model.lifetime_components[i]

        comp.lifetime.err = _get_error(base_idx, err, ctx)

        l_int = comp.intensity
        if l_int.vary
            q_err = _get_error(base_idx + 1, err, ctx)
            q = opt[ctx.signed_param_idcs[base_idx + 1]]

            l_int.err = l_int.val * sqrt(
                l_cumul_sq_err + (q_err / (q + 1e-12)) ^ 2
            )

            l_cumul_sq_err += (q_err / (1 - q + 1e-12)) ^ 2
        elseif i == ctx.last_l_vary_idx
            l_int.err = l_int.val * sqrt(l_cumul_sq_err)
        else
            l_int.err = 0
        end
    end

    r_cumul_sq_err = 0.

    for i in 1:ctx.n_r
        base_idx = 4 + 2 * ctx.n_l + 3 * (i - 1)
        comp = model.resolution_components[i]

        comp.fwhm.err = _get_error(base_idx, err, ctx) * FWHM_OVER_SIGMA

        r_int = comp.intensity
        if r_int.vary
            q_err = _get_error(base_idx + 1, err, ctx)
            q = opt[ctx.signed_param_idcs[base_idx + 1]]

            r_int.err = r_int.val * sqrt(
                r_cumul_sq_err + (q_err / (q + 1e-12)) ^ 2
            )

            r_cumul_sq_err += (q_err / (1 - q + 1e-12)) ^ 2
        elseif i == ctx.last_r_vary_idx
            r_int.err = r_int.val * sqrt(r_cumul_sq_err)
        else
            r_int.err = 0
        end

        comp.t0.err = _get_error(base_idx + 2, err, ctx)
    end
end

mutable struct FitResult
    model::LifetimeModel
    n_dpoints::Int
    fit_status::ReturnCode.T
    n_eval::Int
    red_chi_2::Float64
    cov::Union{Nothing, Matrix{Float64}}
end

function _post_fit!(
    sol,
    model::LifetimeModel,
    ctx::FitContext,
    x::Vector{Float64},
    w::Vector{Float64},
)::FitResult
    opt = sol.u
    res = sol.resid

    _paste_fit_result_values!(opt, model, ctx)

    J = _j_m(x, opt, ctx) .* w

    dof = length(res) - length(opt)
    red_chi_2 = sum(abs2, res) / dof

    fit_result = FitResult(
        model,
        length(sol.resid),
        sol.retcode,
        sol.stats.nf,
        red_chi_2,
        nothing,
    )

    try
        Q, R = qr(J)
        Rinv = inv(R)

        cov = red_chi_2 .* (Rinv * Rinv')

        fit_result.cov = cov

        err = sqrt.(diag(cov))

        _paste_fit_result_errors!(err, opt, model, ctx)
    catch
    end

    if ctx.last_l_vary_idx > 0
        model.lifetime_components[ctx.last_l_vary_idx].intensity.vary = true
    end
    if ctx.last_r_vary_idx > 0
        model.resolution_components[ctx.last_r_vary_idx].intensity.vary = true
    end

    fit_result
end

function fit_lifetime_spectrum(
    x::Vector{Float64},
    y::Vector{Float64},
    model::LifetimeModel,
    counts::Int,
    peak_center::Float64
)::FitResult
    ctx = _prepare_fit!(model, x, counts, peak_center)

    w = @. 1 / sqrt(max(y, 1))

    m(u, p) = w .* (_m(p, u, ctx) .- y)
    # j_m(u, p) = _j_m(p, u, ctx) .* w

    prob = NonlinearLeastSquaresProblem(
        m, ctx.init, x; lb=ctx.lower_bounds, ub=ctx.upper_bounds,
    )

    sol = solve(prob, TrustRegion(; autodiff = AutoForwardDiff()))

    _post_fit!(sol, model, ctx, x, w)
end
