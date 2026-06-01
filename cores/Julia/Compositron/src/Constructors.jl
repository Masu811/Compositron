import .Importers: import_SLOPE_n42

function Measurement(filename::String)
    import_SLOPE_n42(filename)
end
