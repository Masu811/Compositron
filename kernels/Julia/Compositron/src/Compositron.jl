module Compositron

module Constants
    include("Constants.jl")
end

module Utils
    include("Core/Utils.jl")
    include("Core/Fitting.jl")
end

module DBS
    include("DBS/Fitting.jl")
    include("DBS/DBSpectrum.jl")
end

module CDBS
    include("CDBS/CDBSpectrum.jl")
end

module PALS
    include("PALS/Model.jl")
    include("PALS/Fitting.jl")
    include("PALS/PALSpectrum.jl")
end

module Core
    include("Core/Measurement.jl")

    module Importers
        include("Importers/PNG_Importer.jl")
        include("Importers/SLOPE_N42_Importer.jl")
        include("Importers/MePS_DAT_Importer.jl")
    end

    include("Importers/Importers.jl")
end

end # module Compositron
