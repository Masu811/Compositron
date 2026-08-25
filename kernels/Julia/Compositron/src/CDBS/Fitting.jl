using NonlinearSolve

using ..Utils: SimpleFitParam


function gauss2d(
    x,
    y,
    amp,
    x0,
    y0,
    sig_x,
    sig_y,
    phi,
)
    cos_phi = cos(phi)
    sin_phi = sin(phi)

    dx = x - x0
    dy = y - y0

    amp * exp(
        -0.5 * ((dx * cos_phi - dy * sin_phi) / sig_x) ^ 2
        -0.5 * ((dx * sin_phi + dy * cos_phi) / sig_y) ^ 2
    )
end


function fit_gauss2d(
    x::Vector{Float64},
    y::Vector{Float64},
    z::Matrix{Float64},
    init::Vector{Float64},
)::Dict{String, SimpleFitParam}
    xx = vec([xi for _ in y, xi in x])
    yy = vec([yi for yi in y, _ in x])
    zz = vec(z)

    w = @. 1 / sqrt(max(zz, 1))

    m(u, p) = w .* (
        gauss2d.(p[1], p[2], u[1], u[2], u[3], u[4], u[5], u[6]) .- zz
    )

    prob = NonlinearLeastSquaresProblem(m, init, (xx, yy))

    sol = solve(prob, TrustRegion(; autodiff = AutoForwardDiff()))

    opt = sol.u

    Dict(
        "amp" => SimpleFitParam(opt[1], NaN),
        "x0" => SimpleFitParam(opt[2], NaN),
        "y0" => SimpleFitParam(opt[3], NaN),
        "sig_x" => SimpleFitParam(opt[4], NaN),
        "sig_y" => SimpleFitParam(opt[5], NaN),
        "phi" => SimpleFitParam(opt[6], NaN),
    )
end
