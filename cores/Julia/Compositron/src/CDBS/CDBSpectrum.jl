import ..Utils: min_spectrum2d

const Ecal = Tuple{Tuple{Float64, Float64}, Tuple{Float64, Float64}}

mutable struct CDBSpectrum
    spectrum::Matrix{<:Unsigned}
    detpair::String
    ecal::Ecal
    s::Float64
    ds::Float64
    w::Float64
    dw::Float64
    counts::UInt64
    dcounts::Float64
    roi_counts::Float64
    droi_counts::Float64

    function CDBSpectrum(
        spectrum::AbstractMatrix{<:Integer},
        detpair::String,
        ecal::Ecal,
    )
        counts = sum(spectrum)
        dcounts = sqrt(counts)

        c = new()
        c.spectrum = min_spectrum2d(spectrum)
        c.detpair = detpair
        c.ecal = ecal
        c.s = NaN
        c.ds = NaN
        c.w = NaN
        c.dw = NaN
        c.counts = counts
        c.dcounts = dcounts
        c.roi_counts = NaN
        c.droi_counts = NaN

        c
    end
end
