using LinearAlgebra

using ..Constants: PI

export Vertex, ConvexToVerticesCounterClockwise, Polygon, Rectangle,
    Parallelogram, Ellipse, to_vertices, get_bounding_box, agg_aa,
    intersect_convex_polygons

const libagg = joinpath(@__DIR__, "../../../../shared/agg/libagg.so")


struct Vertex
    x::Float64
    y::Float64
end


abstract type ConvexToVerticesCounterClockwise end


mutable struct Polygon <: ConvexToVerticesCounterClockwise
    vertices::Vector{Vertex}
end


function to_vertices(p::Polygon)::Vector{Vertex}
    x0 = sum(v.x for v in p.vertices) / length(p.vertices)
    y0 = sum(v.y for v in p.vertices) / length(p.vertices)

    angles = [atan((v.y - y0), (v.x - x0)) for v in p.vertices]

    p.vertices[sortperm(angles)]
end


struct Rectangle <: ConvexToVerticesCounterClockwise
    i_min::Float64
    i_max::Float64
    j_min::Float64
    j_max::Float64
end


function to_vertices(r::Rectangle)::Vector{Vertex}
    [
        Vertex(r.j_min, r.i_min),
        Vertex(r.j_min, r.i_max),
        Vertex(r.j_max, r.i_max),
        Vertex(r.j_max, r.i_min),
    ]
end


struct Parallelogram <: ConvexToVerticesCounterClockwise
    center_i::Float64
    center_j::Float64
    slope_1::Float64
    slope_2::Float64
    side_length_1::Float64
    side_length_2::Float64
end


function to_vertices(p::Parallelogram)::Vector{Vertex}
    half_width_1 = p.side_length_1 / 2
    half_width_2 = p.side_length_2 / 2

    w = [1., p.slope_1]
    h = [1., p.slope_2]

    w *= half_width_1 / norm(w)
    h *= half_width_2 / norm(h)

    x0 = [p.center_j, p.center_i]

    v1 = x0 + (w - h)
    v2 = x0 + (w + h)
    v3 = x0 + (-w + h)
    v4 = x0 + (-w - h)

    return [
        Vertex(v1[1], v1[2]), Vertex(v2[1], v2[2]),
        Vertex(v3[1], v3[2]), Vertex(v4[1], v4[2]),
    ]
end


struct Ellipse <: ConvexToVerticesCounterClockwise
    center_i::Float64
    center_j::Float64
    radius_i::Float64
    radius_j::Float64
    phi::Float64
end


function to_vertices(e::Ellipse)::Vector{Vertex}
    sin_phi = sin(e.phi)
    cos_phi = cos(e.phi)

    n_vertices = 128

    vertices = []

    for k in 1:n_vertices
        t = 2 * PI * k / n_vertices
        s = sin(t)
        c = cos(t)

        push!(
            vertices,
            Vertex(
                e.center_j + e.radius_i * c * sin_phi
                    + e.radius_j * s * cos_phi,
                e.center_i + e.radius_i * c * cos_phi
                    - e.radius_j * s * sin_phi,
            )
        )
    end

    vertices
end


function get_bounding_box(vertices::Vector{Vertex})::Tuple{Int, Int, Int, Int}
    left_col = Int(round(minimum(v.x for v in vertices)))
    right_col = Int(round(maximum(v.x for v in vertices)))
    lower_row = Int(round(maximum(v.y for v in vertices)))
    upper_row = Int(round(minimum(v.y for v in vertices)))

    (left_col, right_col, lower_row, upper_row)
end


function agg_aa(
    nrows::Int, ncols::Int, polygon::ConvexToVerticesCounterClockwise
)::Matrix{UInt8}
    if nrows < 3 || ncols < 3
        throw(ErrorException("Input data is too small in size"))
    end

    vertices = to_vertices(polygon)

    if length(vertices) < 3
        throw(ErrorException("Insufficient number of vertices"))
    end

    # Transform our coords (pixel centers at integer coords)
    # to agg's coords (pixel boundaries at integer coords)
    vertices = [Vertex(v.x - 0.5, v.y - 0.5) for v in vertices]

    if any(v.x < 0 || v.x > ncols || v.y < 0 || v.y > nrows for v in vertices)
        throw(ErrorException("Polygon vertices are out of matrix bounds"))
    end

    weights = zeros(UInt8, (nrows, ncols))

    vx = [v.x for v in vertices]
    vy = [v.y for v in vertices]

    # Coords are transposed because agg is row-major, Julia is col-major
    ccall(
        (:agg_aa, libagg),
        Cvoid,
        (Ptr{UInt8}, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint),
        weights, ncols, nrows, vy, vx, length(vertices)
    )

    weights
end


function compute_intersection(a::Vertex, b::Vertex, v1::Vertex, v2::Vertex)::Vertex
    m1 = (a.y - b.y) / (a.x - b.x)
    t1 = a.y - m1 * a.x

    m2 = (v1.y - v2.y) / (v1.x - v2.x)
    t2 = v1.y - m2 * v1.x

    x0 = (t2 - t1) / (m1 - m2)
    y0 = m1 * x0 + t1

    Vertex(x0, y0)
end


function inside(v::Vertex, v1::Vertex, v2::Vertex)::Bool
    a = [v1.x - v2.x, v1.y - v2.y]
    b = [v1.x - v.x, v1.y - v.y]

    a[1] * b[2] - a[2] * b[1] >= 0
end


function intersect_convex_polygons(
    a::ConvexToVerticesCounterClockwise,
    b::ConvexToVerticesCounterClockwise,
)::Polygon
    va = to_vertices(a)
    vb = to_vertices(b)

    if length(va) < 3 || length(vb) < 3
        throw(ErrorException("Insufficient number of vertices"))
    end

    la = length(va)

    output = vb

    for i in eachindex(va)
        edge_start = va[i]
        edge_end = va[i % la + 1]

        input = output
        output = []

        for j in eachindex(input)
            curr = input[j]
            prev = input[(j + length(input) - 2) % length(input) + 1]

            curr_inside = inside(curr, edge_start, edge_end)
            prev_inside = inside(prev, edge_start, edge_end)

            if curr_inside
                if !prev_inside
                    push!(
                        output,
                        compute_intersection(prev, curr, edge_start, edge_end)
                    )
                end
                push!(output, curr)
            elseif prev_inside
                push!(
                    output,
                    compute_intersection(prev, curr, edge_start, edge_end)
                )
            end
        end
    end

    Polygon(output)
end
