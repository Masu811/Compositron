using Format

using ..Utils: TimingDetectorPair, min_spectrum, from_index
using ..DBS: fit_gaussian
using ..Constants: FWHM_OVER_SIGMA

mutable struct PALSpectrum
    spectrum::Vector{<:Unsigned}
    detpair::TimingDetectorPair
    counts::Int64
    dcounts::Float64
    peak_bnds::Union{Nothing, Tuple{Int, Int}}
    peak_center::Union{Nothing, Float64}
    dpeak_center::Union{Nothing, Float64}
    peak_fwhm::Union{Nothing, Float64}
    dpeak_fwhm::Union{Nothing, Float64}
    fit_result::Union{Nothing, FitResult}

    function PALSpectrum(
        spectrum::AbstractVector{<:Integer}, detpair::TimingDetectorPair
    )
        counts = sum(spectrum)
        dcounts = sqrt(counts)

        new(
            min_spectrum(spectrum),
            detpair,
            counts,
            dcounts,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
        )
    end
end

function _get_peak_center!(p::PALSpectrum)::Float64
    max_idx = argmax(p.spectrum)

    peak_center_window = 15

    a = max(max_idx - peak_center_window, 0)
    b = min(max_idx + peak_center_window, length(p.spectrum))

    y = Float64.(p.spectrum[a:b])
    x = Float64.(a:b)

    params = fit_gaussian(x, y, [maximum(y), peak_center_window, 100])

    peak_center = params["x0_1"]
    peak_width = params["sig_1"]

    p.peak_center = peak_center.val + a
    p.dpeak_center = peak_center.err
    p.peak_fwhm = peak_width.val * FWHM_OVER_SIGMA
    p.dpeak_fwhm = peak_width.err * FWHM_OVER_SIGMA

    peak_center.val + a
end

function fit!(
    p::PALSpectrum,
    model::LifetimeModel,
    fit_start_idx::Integer,
    fit_end_idx::Integer,
)
    peak_center = something(p.peak_center, _get_peak_center!(p))

    y = Float64.(p.spectrum[fit_start_idx:fit_end_idx])
    x = fit_start_idx:fit_end_idx .|> i -> from_index(p.detpair.tcal, i)

    p.fit_result = fit_lifetime_spectrum(
        x, y, model, p.counts, from_index(p.detpair.tcal, peak_center)
    )
end

function _print_params(params::Vector{FitParam}, names::Vector{String})::String
    text = ""
    n_rows = length(names)

    values = [
        begin
            if isnan(param.val)
                "NaN"
            elseif startswith(name, "int")
                format("{:.2f}", 100 * param.val)
            elseif abs(param.val) >= 1e4
                format("{:.3e}", param.val)
            else
                format("{:.2f}", param.val)
            end
        end
        for (name, param) in zip(names, params)
    ]

    abs_errors = [
        begin
            if isnan(param.err)
                "NaN"
            elseif startswith(name, "int")
                format("{:.1f}", 100 * param.err)
            elseif abs(param.err) >= 1e4
                format("{:.3e}", param.err)
            else
                format("{:.2f}", param.err)
            end
        end
        for (name, param) in zip(names, params)
    ]

    rel_errors = [
        begin
            x = abs(param.err / param.val * 100)
            if isnan(x)
                "NaN %"
            else
                format("{:.2f} %", x)
            end
        end
        for param in params
    ]

    mins = [
        begin
            if isnan(param.min)
                "NaN"
            elseif isinf(param.min)
                if param.min > 0
                    "Inf"
                else
                    "-Inf"
                end
            elseif startswith(name, "int")
                format("{:.1f}", 100 * param.min)
            else
                format("{:.1f}", param.min)
            end
        end
        for (name, param) in zip(names, params)
    ]

    maxs = [
        begin
            if isnan(param.max)
                "NaN"
            elseif isinf(param.min)
                if param.min > 0
                    "Inf"
                else
                    "-Inf"
                end
            elseif startswith(name, "int")
                format("{:.1f}", 100 * param.max)
            else
                format("{:.1f}", param.max)
            end
        end
        for (name, param) in zip(names, params)
    ]

    fixed = [
        begin
            if param.vary
                ""
            else
                "Fixed"
            end
        end
        for param in params
    ]

    at_bnd = [
        begin
            if abs(param.val - param.min) < 1e-3 || abs(param.val - param.max) < 1e-3
                "At Boundary"
            else
                ""
            end
        end
        for param in params
    ]

    columns = [
        names, values, abs_errors, rel_errors, mins, maxs, fixed, at_bnd
    ]

    col_widths = [maximum(length(row) for row in col) for col in columns]

    headers = [
        "Parameter", "Value", "Stderr (abs)", "Stderr (rel)", "Min", "Max",
        "Fixed", ""
    ]

    head_widths = [length(head) for head in headers]

    widths = [max(c, h) for (c, h) in zip(col_widths, head_widths)]

    formats = ["<", ">", ">", ">", ">", ">", "<", "<"]

    for (i, head) in enumerate(headers)
        text *= format(format(" {{:^{}}} ", widths[i]), head)
    end

    text *= "\n"

    for wid in widths
        text *= format(format(" {{:-^{}}} ", wid), "")
    end

    text *= "\n"

    for row in 1:n_rows
        for (col, column) in enumerate(columns)
            if formats[col] == "<"
                text *= format(format(" {{:<{}}} ", widths[col]), column[row])
            else
                text *= format(format(" {{:>{}}} ", widths[col]), column[row])
            end
        end
        text *= "\n"
    end

    text
end

function fit_report(p::PALSpectrum)
    if isnothing(p.fit_result)
        throw(ErrorException("No fit has been performed yet"))
    end

    text = ""

    result = something(p.fit_result)
    model = result.model

    n_l = length(model.lifetime_components)
    n_r = length(model.resolution_components)

    thick_hline = format("{:#^100}\n", "")
    thin_hline = format("{:-^100}\n", "")

    text *= thick_hline
    text *= format("{:^100}\n", "Fit report")
    text *= thick_hline
    text *= "\n"
    text *= format("Statistics: {}\n", p.counts)
    text *= format("Number of Data Points: {}\n", result.n_dpoints)
    text *= "\n"
    text *= format("Lifetime Components:   {}\n", n_l)
    text *= format("Resolution Components: {}\n", n_r)
    if model.background_component.vary
        text *= "Background Components: 1\n"
    else
        text *= "Background Components: 0\n"
    end
    text *= "\n"
    text *= "Fit Method: Trust-Region\n"
    text *= format("Fit Status: {}\n", result.fit_status)
    text *= format("Number of Function Evaluations: {}\n", result.n_eval)
    if isnothing(result.red_chi_2)
        text *= "Reduced Chi squared: NaN\n"
    else
        text *= format("Reduced Chi squared: {}\n", result.red_chi_2)
    end
    text *= "\n"
    text *= thin_hline
    text *= format("{:^100}\n", "Parameters")
    text *= thin_hline
    text *= "\n"
    text *= format("{:=^100}\n", " General Parameters ")
    text *= "\n"

    text *= _print_params(
        [model.scale_component, model.shift_component, model.background_component],
        ["N", "t0", "background"]
    )

    text *= "\n"
    text *= format("{:=^100}\n", " Lifetime Components ")
    text *= "\n"

    lt_norm = sum(comp.intensity.val for comp in model.lifetime_components)

    if abs(lt_norm  - 1) > 1e-3
        text *= "!!! intensities don't add up to 100% !!!\n\n"
    end

    params = Vector{FitParam}()
    names = Vector{String}()

    for (i, comp) in enumerate(model.lifetime_components)
        push!(params, comp.lifetime)
        push!(names, "lifetime_$i")
    end
    for (i, comp) in enumerate(model.lifetime_components)
        push!(params, comp.intensity)
        push!(names, "intensity_$i")
    end

    text *= _print_params(params, names)

    text *= "\n"
    text *= format("{:=^100}\n", " Resolution Components ")
    text *= "\n"

    res_norm = sum(comp.intensity.val for comp in model.resolution_components)

    if abs(res_norm  - 1) > 1e-3
        text *= "!!! intensities don't add up to 100% !!!\n\n"
    end

    params = Vector{FitParam}()
    names = Vector{String}()

    for (i, comp) in enumerate(model.resolution_components)
        push!(params, comp.fwhm)
        push!(names, "fwhm_$i")
    end
    for (i, comp) in enumerate(model.resolution_components)
        push!(params, comp.intensity)
        push!(names, "intensity_$i")
    end
    for (i, comp) in enumerate(model.resolution_components)
        push!(params, comp.t0)
        push!(names, "t0_$i")
    end

    text *= _print_params(params, names)
    text *= "\n"
    text *= thin_hline

    text
end
