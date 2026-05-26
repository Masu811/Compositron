import .Measurement: DopplerMeasurement
import .Importer: import_n42

function DopplerMeasurement(filename::String)::DopplerMeasurement
    m = import_n42(filename)

    m.filename = filename
    m.name = filename

    m
end
