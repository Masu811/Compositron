import PNGFiles: load

import ..Utils: min_spectrum2d

function import_png(filename::String)::Matrix{<:Unsigned}
    img = load(filename)
    spectrum = img .|> pixel -> (
        (UInt32(reinterpret(UInt8, pixel.r)) << 16)
        | (UInt32(reinterpret(UInt8, pixel.g)) << 8)
        | UInt32(reinterpret(UInt8, pixel.b))
    )
    min_spectrum2d(spectrum)
end
