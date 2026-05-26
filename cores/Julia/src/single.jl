module Single

import ..Utils: min_spectrum

export SingleSpectrum
export analyze!

mutable struct SingleSpectrum
    spectrum::Vector{<:Unsigned}
    detname::String
    ecal::Tuple{Float64, Float64}
    s::Float64
    ds::Float64
    w::Float64
    dw::Float64
    v2p::Float64
    dv2p::Float64
    counts::UInt64
    dcounts::Float64
    peak_counts::Float64
    dpeak_counts::Float64

    function SingleSpectrum(
        spectrum::AbstractVector{<:Integer},
        detname::String,
        ecal::Tuple{Float64, Float64},
    )
        counts = sum(spectrum)
        dcounts = sqrt(counts)
        s = new()
        s.spectrum = min_spectrum(spectrum)
        s.detname = detname
        s.ecal = ecal
        s.s = NaN
        s.ds = NaN
        s.w = NaN
        s.dw = NaN
        s.v2p = NaN
        s.dv2p = NaN
        s.counts = counts
        s.dcounts = dcounts
        s.peak_counts = NaN
        s.dpeak_counts = NaN

        s
    end
end

function analyze!(s::SingleSpectrum)
end

end  # module Single
