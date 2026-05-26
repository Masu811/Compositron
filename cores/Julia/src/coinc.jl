module Coinc

import ..Utils: min_spectrum2d

export CoincidenceSpectrum

mutable struct CoincidenceSpectrum
    hist::Matrix{<:Unsigned}
    detpair::String
    ecal::Tuple{Tuple{Float64, Float64}, Tuple{Float64, Float64}}
    s::Float64
    ds::Float64
    w::Float64
    dw::Float64
    v2p::Float64
    dv2p::Float64
    counts::UInt64
    dcounts::Float64
    roi_counts::Float64
    droi_counts::Float64

    function CoincidenceSpectrum(
        hist::AbstractMatrix{<:Integer},
        detpair::String,
        ecal::Tuple{Tuple{Float64, Float64}, Tuple{Float64, Float64}}
    )
        counts = sum(hist)
        dcounts = sqrt(counts)
        c = new()
        c.hist = min_spectrum2d(hist)
        c.detpair = detpair
        c.ecal = ecal
        c.s = NaN
        c.ds = NaN
        c.w = NaN
        c.dw = NaN
        c.v2p = NaN
        c.dv2p = NaN
        c.counts = counts
        c.dcounts = dcounts
        c.roi_counts = NaN
        c.droi_counts = NaN

        c
    end
end

end  # module Coinc
