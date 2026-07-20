import ..Utils: min_spectrum

const Ecal = Tuple{Float64, Float64}

mutable struct DBSpectrum
    spectrum::Vector{<:Unsigned}
    detname::String
    ecal::Ecal
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

    function DBSpectrum(
        spectrum::AbstractVector{<:Integer},
        detname::String,
        ecal::Ecal,
    )
        counts = sum(spectrum)
        dcounts = sqrt(counts)

        d = new()
        d.spectrum = min_spectrum(spectrum)
        d.detname = detname
        d.ecal = ecal
        d.s = NaN
        d.ds = NaN
        d.w = NaN
        d.dw = NaN
        d.v2p = NaN
        d.dv2p = NaN
        d.counts = counts
        d.dcounts = dcounts
        d.peak_counts = NaN
        d.dpeak_counts = NaN

        d
    end
end
