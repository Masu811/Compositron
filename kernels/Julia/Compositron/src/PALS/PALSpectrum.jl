using ..Utils: TimingDetectorPair, min_spectrum

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
