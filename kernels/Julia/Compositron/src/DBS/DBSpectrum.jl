using ..Utils: EnergyDetector, min_spectrum, SimpleFitParam

struct LineshapeParam
end

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
