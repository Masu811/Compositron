import DataStructures: OrderedDict

import ..DBS: DBSpectrum
import ..CDBS: CDBSpectrum

struct MeasurementShape
    dbs::Int
    cdbs::Int
end

mutable struct Measurement
    filename::Union{String, Nothing}
    name::Union{String, Nothing}
    dbs::OrderedDict{String, DBSpectrum}
    cdbs::OrderedDict{String, CDBSpectrum}
    metadata::Dict{String, String}

    function Measurement()
        new(
            nothing,
            nothing,
            OrderedDict(),
            OrderedDict(),
            Dict(),
        )
    end
end

function shape(m::Measurement)
    MeasurementShape(length(m.dbs), length(m.cdbs))
end
