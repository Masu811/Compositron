using ..Utils: EnergyDetector, min_spectrum, SimpleFitParam,
    EcalCorrectionOrder, first_order, zeroth_order, no_corr, Unit, KeV, to_kev,
    to_index_float, to_index_rounded, from_index
using ..Constants: M_E_KEV

@enum PeakModel gauss_peak erf_linear_1_gauss_peak erf_linear_2_gauss_peak erf_linear_3_gauss_peak

@enum BinIntegrationScheme const_int_scheme linear_int_scheme

abstract type Area end

struct PeakArea <: Area
    width::Unit
end

struct PeakAreaWithOffset <: Area
    width::Unit
    offset::Unit
end

struct SpectrumArea <: Area
    left_bnd::Float64
    right_bnd::Float64
end

function to_bnds_kev(
    area::PeakArea, eres::Union{Nothing, Float64}
)::Union{Nothing, Tuple{Float64, Float64}}
    width_kev = to_kev(area.width, eres)

    isnothing(width_kev) && return nothing

    (M_E_KEV - width_kev / 2, M_E_KEV + width_kev / 2)
end

function to_bnds_kev(
    area::PeakAreaWithOffset, eres::Union{Nothing, Float64}
)::Union{Nothing, Tuple{Float64, Float64}}
    width_kev = to_kev(area.width, eres)

    isnothing(width_kev) && return nothing

    offset_kev = to_kev(area.offset, eres)

    isnothing(offset_kev) && return nothing

    (M_E_KEV - width_kev / 2 + offset_kev, M_E_KEV + width_kev / 2 + offset_kev)
end

function to_bnds_kev(
    area::SpectrumArea, eres::Union{Nothing, Float64}
)::Union{Nothing, Tuple{Float64, Float64}}
    (area.left_bnd, area.right_bnd)
end

struct LineshapeParam
    val::Float64
    err::Float64
    num::Vector{Area}
    denom::Vector{Area}
    in_corrected_peak::Bool
end

struct LineshapeParamDefinition
    name::String
    num::Vector{Area}
    denom::Vector{Area}
    in_corrected_peak::Bool
end

const STD_LINESHAPE_PARAMS::Vector{LineshapeParamDefinition} = [
    LineshapeParamDefinition(
        "S",
        [PeakArea(KeV(1))],
        [PeakArea(KeV(20))],
        true,
    ),
    LineshapeParamDefinition(
        "W",
        [
            PeakAreaWithOffset(KeV(1), KeV(3))
            PeakAreaWithOffset(KeV(1), KeV(-3))
        ],
        [PeakArea(KeV(20))],
        true,
    ),
    LineshapeParamDefinition(
        "V/P",
        [SpectrumArea(400, 500)],
        [SpectrumArea(501, 521)],
        false,
    ),
    LineshapeParamDefinition(
        "P/T",
        [SpectrumArea(501, 521)],
        [SpectrumArea(-Inf, Inf)],
        false,
    ),
]

mutable struct DBSpectrum
    spectrum::Vector{<:Unsigned}
    detector::EnergyDetector
    counts::UInt64
    dcounts::Float64
    peak::Union{Nothing, Vector{Float64}}
    peak_bnds::Union{Nothing, Tuple{Int, Int}}
    peak_counts::Union{Nothing, Float64}
    dpeak_counts::Union{Nothing, Float64}
    peak_params::Dict{String, SimpleFitParam}
    lineshape_params::Dict{String, LineshapeParam}

    function DBSpectrum(
        spectrum::AbstractVector{<:Integer}, detector::EnergyDetector
    )
        counts = sum(spectrum)
        dcounts = sqrt(counts)

        new(
            min_spectrum(spectrum),
            detector,
            counts,
            dcounts,
            nothing,
            nothing,
            nothing,
            nothing,
            Dict{String, SimpleFitParam}(),
            Dict{String, LineshapeParam}(),
        )
    end
end

function default_analyze!(d::DBSpectrum)
    correct_ecal!(d, first_order)

    extract_peak!(d, 30, erf_linear_1_gauss_peak)

    for param in STD_LINESHAPE_PARAMS
        calc_lineshape_param!(d, param, const_int_scheme)
    end

    d.peak = nothing
    d.peak_bnds = nothing
end

function _get_peak_energies(d::DBSpectrum)::Vector{Float64}
    ecal = something(d.detector.corrected_ecal, d.detector.ecal)
    bnds = something(d.peak_bnds)

    bnds[1]:bnds[2] .|> i -> from_index(ecal, i)
end

function _fit_peak!(d::DBSpectrum, peak_model::PeakModel)
    x = _get_peak_energies(d)
    y = something(d.peak)

    if peak_model == gauss_peak
        d.peak_params = fit_gaussian(x, y, [maximum(y), 511, 1])
    elseif peak_model == erf_linear_1_gauss_peak
        d.peak_params = fit_erf_linear_1_gauss(
            x, y, [maximum(y), 511, 1, -10, 0, 0]
        )
    elseif peak_model == erf_linear_2_gauss_peak
        m = maximum(y)
        d.peak_params = fit_erf_linear_2_gauss(
            x, y, [m / 2, 511, 1, m / 2, 511, 1.5, -10, 0, 0]
        )
    else
        m = maximum(y)
        d.peak_params = fit_erf_linear_3_gauss(
            x, y, [m / 3, 511, 1, m / 3, 511, 1.5, m / 3, 511, 2.0, -10, 0, 0]
        )
    end
end

function _subtract_bg!(d::DBSpectrum, peak_model::PeakModel)
    peak_model == gauss_peak && return

    _fit_peak!(d, peak_model)

    x = _get_peak_energies(d)

    corr = erf_linear_background.(
        x,
        d.peak_params["erf_amp"].val,
        d.peak_params["x0_1"].val,
        d.peak_params["sig_1"].val,
        d.peak_params["lin"].val,
        d.peak_params["const"].val,
    )

    d.peak .-= corr
end

function extract_peak!(d::DBSpectrum, peak_width::Real, peak_model::PeakModel)
    ecal = something(d.detector.corrected_ecal, d.detector.ecal)

    left_peak_idx = to_index_rounded(ecal, M_E_KEV - peak_width / 2)
    right_peak_idx = to_index_rounded(ecal, M_E_KEV + peak_width / 2)

    if right_peak_idx >= length(d.spectrum)
        throw(ErrorException(
            "Attempting to integrate an area larger than the spectrum or "
            * "extracted peak"
        ))
    end

    d.peak_bnds = (left_peak_idx, right_peak_idx)

    d.peak = Float64.(d.spectrum[left_peak_idx:right_peak_idx])

    _subtract_bg!(d, peak_model)
end

function correct_ecal!(d::DBSpectrum, order::EcalCorrectionOrder)
    order == no_corr && return

    if isnothing(d.peak)
        extract_peak!(d, 20, gauss_peak)
    end

    _fit_peak!(d, gauss_peak)

    c = d.peak_params["x0_1"].val

    ecal = d.detector.ecal

    if order == zeroth_order
        ecal.offset += M_E_KEV - c
    else
        ecal.scale *= (M_E_KEV - ecal.offset) / (c - ecal.offset)
    end

    d.detector.corrected_ecal = ecal
end

function _integrate_const(
    spectrum::Vector, lower_ch_bnd::Float64, upper_ch_bnd::Float64
)::Float64
    counts = 0.

    lowest_ch = Int(round(lower_ch_bnd))
    highest_ch = Int(round(upper_ch_bnd))

    counts += sum(Float64.(spectrum[lowest_ch:highest_ch]))

    # Left edge counts

    lowest_ch_counts = Float64(spectrum[lowest_ch])
    x1 = lowest_ch - 0.5
    counts -= lowest_ch_counts * (lower_ch_bnd - x1)

    # Right edge counts

    highest_ch_counts = Float64(spectrum[highest_ch])
    x2 = highest_ch + 0.5
    counts -= highest_ch_counts * (x2 - upper_ch_bnd)

    counts
end

function _integrate_linear(
    spectrum::Vector, lower_ch_bnd::Float64, upper_ch_bnd::Float64
)::Float64
    counts = 0.

    lowest_ch = Int(floor(lower_ch_bnd))
    highest_ch = Int(ceil(upper_ch_bnd))

    counts += sum(Float64.(spectrum[lowest_ch+1:highest_ch-1]))

    counts += 0.5 * spectrum[lowest_ch]
    counts += 0.5 * spectrum[highest_ch]

    # Left edge counts

    x1 = Float64(lowest_ch)

    y1 = Float64(spectrum[lowest_ch])
    y2 = Float64(spectrum[lowest_ch + 1])

    y_l = (y2 - y1) * (lower_ch_bnd - x1) + y1

    counts -= 0.5 * (y1 + y_l) * (lower_ch_bnd - x1)

    # Right edge counts

    x1 = Float64(highest_ch - 1)
    x2 = Float64(highest_ch)

    y1 = Float64(spectrum[highest_ch - 1])
    y2 = Float64(spectrum[highest_ch])

    y_u = (y2 - y1) * (upper_ch_bnd - x1) + y1

    counts -= 0.5 * (y_u + y2) * (x2 - upper_ch_bnd)

    counts
end

function _integrate_generic(
    spectrum::Vector,
    lower_ch_bnd::Float64,
    upper_ch_bnd::Float64,
    bin_integration_scheme::BinIntegrationScheme,
)::Float64
    if isinf(lower_ch_bnd)
        lower_ch_bnd = 1.
    end
    if isinf(upper_ch_bnd)
        upper_ch_bnd = Float64(length(spectrum))
    end

    if (
        lower_ch_bnd < 1 || lower_ch_bnd > length(spectrum) ||
        upper_ch_bnd < 1 || upper_ch_bnd > length(spectrum)
    )
        throw(ErrorException(
            "Attempting to integrate an area larger than the spectrum or "
            * "extracted peak"
        ))
    end

    if bin_integration_scheme == const_int_scheme
        return _integrate_const(spectrum, lower_ch_bnd, upper_ch_bnd)
    else
        return _integrate_linear(spectrum, lower_ch_bnd, upper_ch_bnd)
    end
end

function integrate(
    d::DBSpectrum,
    area::Area,
    in_corrected_peak::Bool,
    bin_integration_scheme::BinIntegrationScheme,
)::Float64
    bnds = to_bnds_kev(area, d.detector.eres)

    if isnothing(bnds)
        throw(ErrorException(
            "Could not convert units of eres to keV due to missing eres"
        ))
    end

    lower_bnd, upper_bnd = bnds

    if lower_bnd >= upper_bnd
        throw(ErrorException(
            "Input parameter violates constaint lower bound < upper bound"
        ))
    end

    if in_corrected_peak && isnothing(d.peak)
        throw(ErrorException(
            "Attempting analysis on background subtracted peak, but no "
            * "subtraction has been performed yet"
        ))
    end

    ecal = something(d.detector.corrected_ecal, d.detector.ecal)

    lower_ch_bnd = to_index_float(ecal, lower_bnd)
    upper_ch_bnd = to_index_float(ecal, upper_bnd)

    if in_corrected_peak
        lower_peak_bnd, _ = something(d.peak_bnds)
        lower_ch_bnd -= lower_peak_bnd
        upper_ch_bnd -= lower_peak_bnd
        peak = something(d.peak)
        return _integrate_generic(
            peak, lower_ch_bnd, upper_ch_bnd, bin_integration_scheme
        )
    else
        return _integrate_generic(
            d.spectrum, lower_ch_bnd, upper_ch_bnd, bin_integration_scheme
        )
    end
end

function _get_area_bnds(
    d::DBSpectrum, areas::Vector{Area}
)::Vector{Tuple{Float64, Float64}}
    bnds = to_bnds_kev.(areas, d.detector.eres)

    if any(isnothing.(bnds))
        throw(ErrorException(
            "Could not convert units of eres to keV due to missing eres"
        ))
    end

    for i in 1:length(areas) - 1
        for j in i+1:length(areas)
            if !(bnds[i][2] < bnds[j][1] || bnds[j][2] < bnds[i][1])
                throw(ErrorException(
                    "Lineshape numerator or denominator areas overlap"
                ))
            end
        end
    end

    bnds
end

function calc_lineshape_param!(
    d::DBSpectrum,
    definition::LineshapeParamDefinition,
    bin_integration_scheme::BinIntegrationScheme,
)::LineshapeParam
    num_bnds = _get_area_bnds(d, definition.num)
    denom_bnds = _get_area_bnds(d, definition.denom)

    num = sum(integrate(
        d, area, definition.in_corrected_peak, bin_integration_scheme
    ) for area in definition.num)

    denom = sum(integrate(
        d, area, definition.in_corrected_peak, bin_integration_scheme
    ) for area in definition.denom)

    common = 0.

    for num_bnd in num_bnds
        for denom_bnd in denom_bnds
            left_bnd = max(num_bnd[1], denom_bnd[1])
            right_bnd = min(num_bnd[2], denom_bnd[2])

            if left_bnd < right_bnd
                common += integrate(
                    d,
                    SpectrumArea(left_bnd, right_bnd),
                    definition.in_corrected_peak,
                    bin_integration_scheme,
                )
            end
        end
    end

    val = num / denom
    err = val * sqrt(
        (num - common) / num^2 +
        (denom - common) / denom^2 +
        common * ((denom - num) / (denom * num))^2
    )

    ls_param = LineshapeParam(
        val, err, definition.num, definition.denom, definition.in_corrected_peak
    )

    d.lineshape_params[definition.name] = ls_param

    return ls_param
end
