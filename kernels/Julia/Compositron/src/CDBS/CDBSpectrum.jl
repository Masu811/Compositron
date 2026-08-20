using ..Utils: EnergyDetectorPair, min_spectrum2d, SimpleFitParam

struct LineshapeParam
end

mutable struct CDBSpectrum
    spectrum::Matrix{<:Unsigned}
    detpair::EnergyDetectorPair
    counts::UInt64
    dcounts::Float64
    peak::Union{Nothing, Matrix{Float64}}
    peak_bnds::Union{Nothing, Tuple{Tuple{Int, Int}, Tuple{Int, Int}}}
    peak_counts::Union{Nothing, Float64}
    dpeak_counts::Union{Nothing, Float64}
    peak_params::Dict{String, SimpleFitParam}
    lineshape_params::Dict{String, LineshapeParam}

    function CDBSpectrum(
        spectrum::AbstractMatrix{<:Integer},
        detpair::EnergyDetectorPair,
    )
        counts = sum(spectrum)
        dcounts = sqrt(counts)

        new(
            min_spectrum2d(spectrum),
            detpair,
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
