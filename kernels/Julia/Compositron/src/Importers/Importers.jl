using .Importers: import_slope_n42, import_meps_dat

function Measurement(filepath::String)
    _, ext = splitext(filepath)

    println("$ext")

    if ext == ".n42"
        import_slope_n42(filepath)
    elseif ext == ".dat"
        import_meps_dat(filepath)
    else
        throw(ErrorException("Invalid filepath"))
    end
end
