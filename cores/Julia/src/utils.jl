module Utils

export min_spectrum, min_spectrum2d

function min_dtype(arr::AbstractArray{<:Integer})::Type{<:Unsigned}
    max_val = maximum(arr)
    if max_val <= typemax(UInt8)
        return UInt8
    elseif max_val <= typemax(UInt16)
        return UInt16
    elseif max_val <= typemax(UInt32)
        return UInt32
    else
        return UInt64
    end
end

function min_spectrum(arr::AbstractVector{<:Integer})::Vector{<:Unsigned}
    dtype = min_dtype(arr)
    Vector{dtype}(max.(arr, 0))
end

function min_spectrum2d(arr::AbstractMatrix{<:Integer})::Matrix{<:Unsigned}
    dtype = min_dtype(arr)
    Matrix{dtype}(max.(arr, 0))
end

end  # module Utils
