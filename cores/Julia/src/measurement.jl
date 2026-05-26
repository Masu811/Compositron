module Measurement

using DataStructures

import ..Single.SingleSpectrum
import ..Coinc.CoincidenceSpectrum

export DopplerMeasurement

mutable struct DopplerMeasurement
    name::Union{String, Nothing}
    filename::Union{String, Nothing}
    singles::OrderedDict{String, SingleSpectrum}
    coinc::OrderedDict{String, CoincidenceSpectrum}
    metadata::Dict{String, String}

    function DopplerMeasurement()
        new(
            nothing,
            nothing,
            OrderedDict(),
            OrderedDict(),
            Dict(),
        )
    end
end

struct DopplerMeasurementShape
    s::Int
    c::Int
end

function shape(m::DopplerMeasurement)
    s = length(m.singles)
    c = length(m.coinc)
    DopplerMeasurementShape(s, c)
end

end  # module Measurement
