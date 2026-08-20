using DataStructures: OrderedDict

using ..DBS: DBSpectrum
using ..CDBS: CDBSpectrum
using ..PALS: PALSpectrum

struct MeasurementShape
    dbs::Int
    cdbs::Int
    pals::Int
end

mutable struct Measurement
    filename::Union{String, Nothing}
    name::Union{String, Nothing}
    dbs::OrderedDict{String, DBSpectrum}
    cdbs::OrderedDict{String, CDBSpectrum}
    pals::OrderedDict{String, PALSpectrum}
    metadata::Dict{String, String}

    function Measurement()
        new(
            nothing,
            nothing,
            OrderedDict(),
            OrderedDict(),
            OrderedDict(),
            Dict(),
        )
    end
end

function shape(m::Measurement)
    MeasurementShape(length(m.dbs), length(m.cdbs), length(m.pals))
end
