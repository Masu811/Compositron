using ..Utils
using ..DBS
using ..Constants

export BackgroundModel, BGModel_None, Axis, Axis_FirstDet, Axis_SecondDet,
    Area, DiagonalArea, DiagonalAreaWithOffset, AxisAlignedArea, EllipseArea,
    ProjectionBins, LinearProjectionBins, CustomProjectionBins, LineshapeParam,
    LineshapeParamDefinition, STD_LINESHAPE_PARAMS, FoldedProjection,
    Projection, get_energies, fold, to_dbspectrum, CDBSpectrum, project_axes,
    extract_peak!, correct_ecal!, integrate, calc_lineshape_param!,
    default_analyze!, project_diagonal


@enum BackgroundModel begin
    BGModel_None
end


@enum Axis begin
    Axis_FirstDet
    Axis_SecondDet
end


abstract type Area end


struct DiagonalArea <: Area
    width_cel::Unit
    width_cml::Unit
end


struct DiagonalAreaWithOffset <: Area
    width_cel::Unit
    width_cml::Unit
    offset_cel::Unit
    offset_cml::Unit
end


struct AxisAlignedArea <: Area
    first_det_bnds::Tuple{Float64, Float64}
    second_det_bnds::Tuple{Float64, Float64}
end


struct EllipseArea <: Area
    radius_cel::Unit
    radius_cml::Unit
end


abstract type ProjectionBins end


struct LinearProjectionBins <: ProjectionBins
    bin_width::Unit
end


struct CustomProjectionBins <: ProjectionBins
    bins::AbstractVector{Float64}
end


struct LineshapeParam
    val::Float64
    err::Float64
    num::Vector{Area}
    denom::Vector{Area}
    in_corrected_peak::Bool
end


struct LineshapeParamDefinition
    name::String
    num::Vector{Area}
    denom::Vector{Area}
    in_corrected_peak::Bool
end


const STD_LINESHAPE_PARAMS::Vector{LineshapeParamDefinition} = [
    LineshapeParamDefinition(
        "S",
        [DiagonalArea(KeV(2), Eres(1))],
        [DiagonalArea(KeV(Inf), Eres(1))],
        true
    ),
    LineshapeParamDefinition(
        "W",
        [
            DiagonalAreaWithOffset(KeV(1), Eres(1), KeV(3), KeV(0)),
            DiagonalAreaWithOffset(KeV(1), Eres(1), KeV(-3), KeV(0)),
        ],
        [DiagonalArea(KeV(Inf), Eres(1))],
        true
    ),
    LineshapeParamDefinition(
        "P/T",
        [DiagonalArea(KeV(Inf), Eres(1))],
        [AxisAlignedArea((-Inf, Inf), (-Inf, Inf))],
        false
    ),
]


function ij_to_xy(i, ecal::LinearCalibration)::Float64
    from_index(ecal, i)
end


function xy_to_ij(x, ecal::LinearCalibration)::Float64
    to_index_float(ecal, x)
end


function xy_to_uv(x::Float64, y::Float64)::Tuple{Float64, Float64}
    (0.5 * (x - y), 0.5 * (x + y) - M_E_KEV)
end


function uv_to_xy(u::Float64, v::Float64)::Tuple{Float64, Float64}
    return (M_E_KEV + u + v, M_E_KEV - u + v)
end


function ij_to_uv(
    i, j, ecal::Tuple{LinearCalibration, LinearCalibration}
)::Tuple{Float64, Float64}
    x, y = (ij_to_xy(j, ecal[2]), ij_to_xy(i, ecal[1]))
    xy_to_uv(x, y)
end


function uv_to_ij(
    u, v, ecal::Tuple{LinearCalibration, LinearCalibration}
)::Tuple{Float64, Float64}
    x, y = uv_to_xy(u, v)
    (xy_to_ij(y, ecal[1]), xy_to_ij(x, ecal[2]))
end


function convert_rectangle(
    first_det_bnds::Tuple{Float64, Float64},
    second_det_bnds::Tuple{Float64, Float64},
    ecal::Tuple{LinearCalibration, LinearCalibration},
    peak_bnds::Tuple{Tuple{Int, Int}, Tuple{Int, Int}},
)::Rectangle
    Rectangle(
        to_index_float(ecal[1], first_det_bnds[2]) - peak_bnds[1][1] + 1,
        to_index_float(ecal[1], first_det_bnds[1]) - peak_bnds[1][1] + 1,
        to_index_float(ecal[2], second_det_bnds[1]) - peak_bnds[2][1] + 1,
        to_index_float(ecal[2], second_det_bnds[2]) - peak_bnds[2][1] + 1,
    )
end


function convert_parallelogram(
    width_cel::Float64,
    width_cml::Float64,
    offset_cel::Float64,
    offset_cml::Float64,
    ecal::Tuple{LinearCalibration, LinearCalibration},
    peak_bnds::Tuple{Tuple{Int, Int}, Tuple{Int, Int}},
)::Parallelogram
    i, j = uv_to_ij(offset_cel, offset_cml, ecal)
    w1 = sqrt(0.5 * (
        (width_cel / ecal[2].scale)^2 + (width_cel / ecal[1].scale)^2
    ))
    w2 = sqrt(0.5 * (
        (width_cml / ecal[2].scale)^2 + (-width_cml / ecal[1].scale)^2
    ))
    Parallelogram(
        i - peak_bnds[1][1] + 1,
        j - peak_bnds[2][1] + 1,
        -ecal[2].scale / ecal[1].scale,
        ecal[2].scale / ecal[1].scale,
        w1,
        w2,
    )
end


function convert_ellipse(
    radius_cel::Float64,
    radius_cml::Float64,
    ecal::Tuple{LinearCalibration, LinearCalibration},
    peak_bnds::Tuple{Tuple{Int, Int}, Tuple{Int, Int}},
)::Ellipse
    w1 = sqrt(0.5 * (
        (radius_cel / ecal[2].scale)^2 + (radius_cel / ecal[1].scale)^2
    ))
    w2 = sqrt(0.5 * (
        (radius_cml / ecal[2].scale)^2 + (-radius_cml / ecal[1].scale)^2
    ))
    Ellipse(
        to_index_float(ecal[1], M_E_KEV) - peak_bnds[1][1] + 1,
        to_index_float(ecal[2], M_E_KEV) - peak_bnds[2][1] + 1,
        w1,
        w2,
        -atan(ecal[2].scale / ecal[1].scale),
    )
end


struct FoldedProjection
    energies::Vector{Float64}
    spectrum::Vector{Float64}
end


struct Projection
    spectrum::Vector{Float64}
    parent_detector_name::String
    ecal::Union{Nothing, LinearCalibration}
    bins::Union{Nothing, Vector{Float64}}
    counts::Float64

    function Projection(
        spectrum::Vector{Float64},
        parent_detector_name::String,
        ecal::Union{Nothing, LinearCalibration},
        bins::Union{Nothing, Vector{Float64}},
    )
        new(
            spectrum,
            parent_detector_name,
            ecal,
            bins,
            sum(spectrum),
        )
    end
end


function get_energies(p::Projection)::Vector{Float64}
    if !isnothing(p.ecal)
        1:length(p.spectrum) .|> i -> from_index(p.ecal, i)
    elseif !isnothing(p.bins)
        @. 0.5 * (p.bins[1:end-1] + p.bins[2:end])
    else
        throw(ErrorException("Projection is missing ecal and bins"))
    end
end


function fold(p::Projection)::FoldedProjection
    if !isnothing(p.ecal)
        center_idx = to_index_float(p.ecal, 0)

        if abs(center_idx % 1 - 0.5) > 1e-2
            throw(ErrorException(
                "Could not fold projection as bins are not symmetric around 0"
            ))
        end

        last_left_idx = to_index_rounded(p.ecal, 0)
        first_right_idx = last_left_idx + 1

        len_left = last_left_idx
        len_right = length(p.spectrum) - len_left

        len_new = min(len_left, len_right)

        spectrum = zeros(Float64, len_new)

        spectrum .+= p.spectrum[first_right_idx:first_right_idx + len_new]
        spectrum .+= p.spectrum[last_left_idx:-1:last_left_idx - len_new]

        energies = (
            first_right_idx:first_right_idx + len_new
            .|> i -> from_index(p.ecal, i)
        )

        FoldedProjection(energies, spectrum)
    elseif isnothing(p.bins)
        center_idx = argmin(abs(p.bins))

        len_left = center_idx
        len_right = length(p.spectrum) - len_left

        len_new = min(len_left, len_right)

        spectrum = zeros(Float64, len_new)
        energies = zeros(Float64, len_new)

        for i in 1:len_new
            left_bin_left_edge = p.bins[center_idx - i - 1]
            left_bin_right_edge = p.bins[center_idx - i]

            right_bin_left_edge = p.bins[center_idx + i]
            right_bin_right_edge = p.bins[center_idx + i + 1]

            if (
                left_bin_left_edge + right_bin_right_edge > 1e-2 ||
                left_bin_right_edge + right_bin_left_edge > 1e-2
            )
                throw(ErrorException(
                    "Could not fold projection as bins are not symmetric around 0"
                ))
            end

            spectrum[i] = (
                p.spectrum[center_idx - i - 1] + p.spectrum[center_idx + i]
            )
            energies[i] = 0.5 * (right_bin_left_edge + right_bin_right_edge)
        end

        FoldedProjection(energies, spectrum)
    else
        throw(ErrorException("Projection is missing ecal and bins"))
    end
end


function to_dbspectrum(p::Projection)::DBSpectrum
    if isnothing(p.ecal)
        throw(ErrorException(
                "Could not convert projection to DBSpectrum due to missing "
                * "energy calibration"
        ))
    end

    DBSpectrum(
        UInt64.(p.spectrum),  # This will throw an InexactError
        EnergyDetector(
            p.parent_detector_name,
            p.ecal,
            nothing,
            nothing
        )
    )
end


mutable struct CDBSpectrum
    spectrum::Matrix{<:Unsigned}
    detpair::EnergyDetectorPair
    counts::UInt64
    peak::Union{Nothing, Matrix{Float64}}
    peak_bnds::Union{Nothing, Tuple{Tuple{Int, Int}, Tuple{Int, Int}}}
    peak_counts::Union{Nothing, Float64}
    peak_params::Dict{String, SimpleFitParam}
    lineshape_params::Dict{String, LineshapeParam}

    function CDBSpectrum(
        spectrum::AbstractMatrix{<:Integer},
        detpair::EnergyDetectorPair,
    )
        new(
            min_spectrum2d(spectrum),
            detpair,
            sum(spectrum),
            nothing,
            nothing,
            nothing,
            Dict{String, SimpleFitParam}(),
            Dict{String, LineshapeParam}(),
        )
    end
end


function project_axes(c::CDBSpectrum, onto_axis::Axis)::Projection
    if onto_axis == Axis_FirstDet
        Projection(
            vec(Float64.(sum(c.spectrum; dims=1))),
            c.detpair.first_det.name,
            c.detpair.first_det.ecal,
            nothing,
        )
    else
        Projection(
            vec(Float64.(sum(c.spectrum; dims=2))),
            c.detpair.second_det.name,
            c.detpair.second_det.ecal,
            nothing,
        )
    end
end


function _fit_2d_peak!(c::CDBSpectrum, bg_model::BackgroundModel)
    if isnothing(c.peak) || isnothing(c.peak_bnds)
        throw(ErrorException(
            "Attempting analysis on background subtracted peak, but no "
            * "subtraction has been performed yet"
        ))
    end

    ecal_1 = something(
        c.detpair.first_det.corrected_ecal, c.detpair.first_det.ecal
    )
    ecal_2 = something(
        c.detpair.second_det.corrected_ecal, c.detpair.second_det.ecal
    )

    x = c.peak_bnds[2][1]:c.peak_bnds[2][2] .|> i -> from_index(ecal_2, i)
    y = c.peak_bnds[1][1]:c.peak_bnds[1][2] .|> i -> from_index(ecal_1, i)

    if bg_model == BGModel_None
        i0, j0 = Tuple(argmax(c.peak))

        max_val = c.peak[i0, j0]

        y0 = from_index(ecal_1, i0 + c.peak_bnds[1][1] - 1)
        x0 = from_index(ecal_2, j0 + c.peak_bnds[2][1] - 1)

        c.peak_params = fit_gauss2d(
            x, y, c.peak, [max_val, x0, y0, 0.7, 1.8, -0.78],
        )
    end
end


function _subtract_bg!(c::CDBSpectrum, bg_model::BackgroundModel)
    if bg_model == BGModel_None
        return
    end
end


function extract_peak!(
    c::CDBSpectrum,
    peak_window_kev::Tuple{Float64, Float64},
    bg_model::BackgroundModel,
)
    ecal_1 = something(
        c.detpair.first_det.corrected_ecal, c.detpair.first_det.ecal
    )
    ecal_2 = something(
        c.detpair.second_det.corrected_ecal, c.detpair.second_det.ecal
    )

    first_row_e = to_index_float(ecal_1, M_E_KEV - peak_window_kev[1] / 2)
    last_row_e = to_index_float(ecal_1, M_E_KEV + peak_window_kev[1] / 2)

    first_col_e = to_index_float(ecal_2, M_E_KEV - peak_window_kev[2] / 2)
    last_col_e = to_index_float(ecal_2, M_E_KEV + peak_window_kev[2] / 2)

    first_row = min(max(Int(floor(first_row_e)), 1), size(c.spectrum)[1])
    last_row = min(max(Int(floor(last_row_e)), 1), size(c.spectrum)[1])

    first_col = min(max(Int(floor(first_col_e)), 1), size(c.spectrum)[2])
    last_col = min(max(Int(floor(last_col_e)), 1), size(c.spectrum)[2])

    c.peak = Float64.(c.spectrum[first_row:last_row, first_col:last_col])

    c.peak_bnds = ((first_row, last_row), (first_col, last_col))

    _subtract_bg!(c, bg_model)

    peak_counts = sum(c.peak)

    c.peak_counts = peak_counts
end


function correct_ecal!(c::CDBSpectrum, order::EcalCorrectionOrder)
    if order == EcalCorrOrder_None
        return
    end

    if isnothing(c.peak)
        extract_peak!(c, (4., 4.), BGModel_None)
    end

    if !haskey(c.peak_params, "x0") || !haskey(c.peak_params, "y0")
        _fit_2d_peak!(c, BGModel_None)
    end

    x0 = c.peak_params["x0"].val
    y0 = c.peak_params["y0"].val

    ecal_1 = copy(c.detpair.first_det.ecal)
    ecal_2 = copy(c.detpair.second_det.ecal)

    if order == EcalCorrOrder_Zeroth
        ecal_1.offset += M_E_KEV - y0
        ecal_2.offset += M_E_KEV - x0
    else
        ecal_1.scale *= (M_E_KEV - ecal_1.offset) / (y0 - ecal_1.offset)
        ecal_2.scale *= (M_E_KEV - ecal_2.offset) / (x0 - ecal_2.offset)
    end

    c.detpair.first_det.corrected_ecal = ecal_1
    c.detpair.second_det.corrected_ecal = ecal_2
end


function _get_max_projection_length(
    v::Float64,
    ecal::Tuple{LinearCalibration, LinearCalibration},
    peak_bnds::Tuple{Tuple{Int, Int}, Tuple{Int, Int}},
)::Float64
    i1, j1 = uv_to_ij(0., v, ecal)
    i2, j2 = uv_to_ij(0., -v, ecal)
    m = -ecal[2].scale / ecal[1].scale
    t1 = i1 - m * j1
    t2 = i2 - m * j2

    ((i_min, i_max), (j_min, j_max)) = peak_bnds

    i_min = i_min - 0.49
    i_max = i_max + 0.49
    j_min = j_min - 0.49
    j_max = j_max + 0.49

    x1 = ij_to_uv(i_min, (i_min - t1) / m, ecal)
    x2 = ij_to_uv(i_min, (i_min - t2) / m, ecal)
    x3 = ij_to_uv(m * j_min + t1, j_min, ecal)
    x4 = ij_to_uv(m * j_min + t2, j_min, ecal)
    x5 = ij_to_uv(i_max, (i_max - t1) / m, ecal)
    x6 = ij_to_uv(i_max, (i_max - t2) / m, ecal)
    x7 = ij_to_uv(m * j_max + t1, j_max, ecal)
    x8 = ij_to_uv(m * j_max + t2, j_max, ecal)

    minimum(abs(v[1]) for v in [x1, x2, x3, x4, x5, x6, x7, x8])
end


function _calculate_polygon_boundaries(
    c::CDBSpectrum,
    area::Area,
    ecal::Tuple{LinearCalibration, LinearCalibration},
    peak_bnds::Tuple{Tuple{Int, Int}, Tuple{Int, Int}},
)::ConvexToVerticesCounterClockwise
    if isa(area, DiagonalArea)
        width_cel = to_kev(area.width_cel, c.detpair.eres)
        width_cml = to_kev(area.width_cml, c.detpair.eres)

        if isinf(width_cel)
            width_cel = _get_max_projection_length(
                0.5 * width_cml, ecal, peak_bnds
            )
        end

        if any(!isfinite(x) for x in [width_cel, width_cml])
            throw(ErrorException("Invalid boundaries for area"))
        end

        convert_parallelogram(width_cel, width_cml, 0., 0., ecal, peak_bnds)
    elseif isa(area, DiagonalAreaWithOffset)
        width_cel = to_kev(area.width_cel, c.detpair.eres)
        width_cml = to_kev(area.width_cml, c.detpair.eres)
        offset_cel = to_kev(area.offset_cel, c.detpair.eres)
        offset_cml = to_kev(area.offset_cml, c.detpair.eres)

        if isinf(width_cel)
            width_cel = _get_max_projection_length(
                0.5 * width_cml, ecal, peak_bnds
            )
        end

        if any(!isfinite(x) for x in [
            width_cel, width_cml, offset_cel, offset_cml
        ])
            throw(ErrorException("Invalid boundaries for area"))
        end

        convert_parallelogram(
            width_cel, width_cml, offset_cel, offset_cml, ecal, peak_bnds
        )
    elseif isa(area, AxisAlignedArea)
        y_min, y_max = area.first_det_bnds
        x_min, x_max = area.second_det_bnds

        if isinf(y_min)
            y_min = from_index(ecal[1], peak_bnds[1][1] - 0.49)
        end
        if isinf(y_max)
            y_max = from_index(ecal[1], peak_bnds[1][2] + 0.49)
        end
        if isinf(x_min)
            x_min = from_index(ecal[2], peak_bnds[2][1] - 0.49)
        end
        if isinf(x_max)
            x_max = from_index(ecal[2], peak_bnds[2][2] + 0.49)
        end

        if any(
            !isfinite(x) for x in [y_min, y_max, x_min, x_max]
        )
            throw(ErrorException("Invalid boundaries for area"))
        end

        convert_rectangle((y_min, y_max), (x_min, x_max), ecal, peak_bnds)
    else
        radius_cel = to_kev(area.radius_cel, c.detpair.eres)
        radius_cml = to_kev(area.radius_cml, c.detpair.eres)

        if !isfinite(radius_cel) || !isfinite(radius_cml)
            throw(ErrorException("Invalid boundaries for area"))
        end

        convert_ellipse(radius_cel, radius_cml, ecal, peak_bnds)
    end
end


function _blend_and_sum(
    c::CDBSpectrum,
    weights::Matrix{UInt8},
    upper_row::Int,
    left_col::Int,
    nrows::Int,
    ncols::Int,
    in_corrected_peak::Bool,
)::Float64
    if in_corrected_peak
        if isnothing(c.peak)
            throw(ErrorException(
                "Attempting analysis on background subtracted peak, but no "
                * "subtraction has been performed yet"
            ))
        end

        view = c.peak[upper_row:upper_row+nrows-1, left_col:left_col+ncols-1]

        sum(view .* weights) / 255.
    else
        view = c.spectrum[upper_row:upper_row+nrows-1, left_col:left_col+ncols-1]

        sum(UInt64.(view) .* weights) / 255.
    end
end


function integrate(c::CDBSpectrum, area::Area, in_corrected_peak::Bool)::Float64
    ecal_1 = something(
        c.detpair.first_det.corrected_ecal, c.detpair.first_det.ecal
    )
    ecal_2 = something(
        c.detpair.second_det.corrected_ecal, c.detpair.second_det.ecal
    )
    ecal = (ecal_1, ecal_2)

    integral_bnds = if in_corrected_peak
        if isnothing(c.peak) || isnothing(c.peak_bnds)
            throw(ErrorException(
                "Attempting analysis on background subtracted peak, but no "
                * "subtraction has been performed yet"
            ))
        end
        c.peak_bnds
    else
        ((1, size(c.spectrum)[1]), (1, size(c.spectrum)[2]))
    end

    nrows = integral_bnds[1][2] - integral_bnds[1][1] + 1
    ncols = integral_bnds[2][2] - integral_bnds[2][1] + 1

    polygon = _calculate_polygon_boundaries(c, area, ecal, integral_bnds)

    vertices = to_vertices(polygon)

    left_col, right_col, lower_row, upper_row = get_bounding_box(vertices)

    if right_col > ncols - 1 || lower_row > nrows - 1
        throw(ErrorException("Polygon vertices are out of matrix bounds"))
    end

    vertices = [Vertex(v.x - left_col + 1, v.y - upper_row + 1) for v in vertices]

    nrows_view = lower_row - upper_row + 1
    ncols_view = right_col - left_col + 1

    weights = agg_aa(nrows_view, ncols_view, Polygon(vertices))

    _blend_and_sum(
        c,
        weights,
        upper_row,
        left_col,
        nrows_view,
        ncols_view,
        in_corrected_peak
    )
end


function calc_lineshape_param!(
    c::CDBSpectrum,
    definition::LineshapeParamDefinition,
)::LineshapeParam
    ecal_1 = something(
        c.detpair.first_det.corrected_ecal, c.detpair.first_det.ecal
    )
    ecal_2 = something(
        c.detpair.second_det.corrected_ecal, c.detpair.second_det.ecal
    )
    ecal = (ecal_1, ecal_2)

    integral_bnds = if definition.in_corrected_peak
        if isnothing(c.peak) || isnothing(c.peak_bnds)
            throw(ErrorException(
                "Attempting analysis on background subtracted peak, but no "
                * "subtraction has been performed yet"
            ))
        end
        c.peak_bnds
    else
        ((1, size(c.spectrum)[1]), (1, size(c.spectrum)[2]))
    end

    num_polygons = [
        _calculate_polygon_boundaries(c, area, ecal, integral_bnds)
        for area in definition.num
    ]
    denom_polygons = [
        _calculate_polygon_boundaries(c, area, ecal, integral_bnds)
        for area in definition.denom
    ]

    for field in [num_polygons, denom_polygons]
        for i in 1:length(field) - 1
            for j in i+1:length(field)
                intersection = intersect_convex_polygons(
                    field[i], field[j]
                )

                if length(intersection.vertices) > 0
                    throw(ErrorException(
                        "Lineshape numerator or denominator areas overlap"
                    ))
                end
            end
        end
    end

    num_areas = [
        integrate(c, area, definition.in_corrected_peak)
        for area in definition.num
    ]
    denom_areas = [
        integrate(c, area, definition.in_corrected_peak)
        for area in definition.denom
    ]

    common_areas = []

    for num_area in num_polygons
        for denom_area in denom_polygons
            intersection = intersect_convex_polygons(num_area, denom_area)

            if length(intersection.vertices) == 0
                continue
            end

            left_col, right_col, lower_row, upper_row = get_bounding_box(
                to_vertices(intersection)
            )

            intersection.vertices = [
                Vertex(v.x - left_col + 1, v.y - upper_row + 1)
                for v in intersection.vertices
            ]

            nrows_view = lower_row - upper_row + 1
            ncols_view = right_col - left_col + 1

            weights = agg_aa(nrows_view, ncols_view, intersection)

            integral = _blend_and_sum(
                c,
                weights,
                upper_row,
                left_col,
                nrows_view,
                ncols_view,
                definition.in_corrected_peak,
            )

            push!(common_areas, integral)
        end
    end

    num = sum(num_areas; init=zero(Float64))
    denom = sum(denom_areas; init=zero(Float64))
    common = sum(common_areas; init=zero(Float64))

    val = num / denom
    err = val * sqrt(
        (num - common) / num^2 +
        (denom - common) / denom^2 +
        common * ((denom - num) / (denom * num))^2
    )

    ls_param = LineshapeParam(
        val,
        err,
        deepcopy(definition.num),
        deepcopy(definition.denom),
        definition.in_corrected_peak,
    )

    c.lineshape_params[definition.name] = ls_param

    ls_param
end


function default_analyze!(c::CDBSpectrum)
    correct_ecal!(c, EcalCorrOrder_First)

    extract_peak!(c, (10., 10.), BGModel_None)

    for param in STD_LINESHAPE_PARAMS
        calc_lineshape_param!(c, param)
    end

    c.peak = nothing
    c.peak_bnds = nothing
end


function project_diagonal(
    c::CDBSpectrum,
    bins::ProjectionBins,
    width::Unit,
    max_length::Unit,
    in_corrected_peak::Bool,
)::Projection
    ecal_1 = something(
        c.detpair.first_det.corrected_ecal, c.detpair.first_det.ecal
    )
    ecal_2 = something(
        c.detpair.second_det.corrected_ecal, c.detpair.second_det.ecal
    )
    ecal = (ecal_1, ecal_2)

    integral_bnds = if in_corrected_peak
        if isnothing(c.peak) || isnothing(c.peak_bnds)
            throw(ErrorException(
                "Attempting analysis on background subtracted peak, but no "
                * "subtraction has been performed yet"
            ))
        end
        c.peak_bnds
    else
        ((1, size(c.spectrum)[1]), (1, size(c.spectrum)[2]))
    end

    width_kev = to_kev(width, c.detpair.eres)
    max_length_kev = to_kev(max_length, c.detpair.eres)

    if isa(bins, LinearProjectionBins)
        u_bin = to_kev(bins.bin_width, c.detpair.eres)
        u_max = _get_max_projection_length(
            width_kev / 2, ecal, integral_bnds
        )

        u_max = min(u_max, max_length_kev)

        num_bins = Int(floor(2 * u_max / u_bin))

        if num_bins < 2
            throw(ErrorException("Invalid bins for projection"))
        end

        max_offset = (num_bins / 2 - 0.5) * u_bin

        ecal = LinearCalibration(-max_offset, u_bin)

        spectrum = [
            integrate(
                c,
                DiagonalAreaWithOffset(
                    KeV(u_bin), width, KeV(from_index(ecal, i)), KeV(0)
                ),
                in_corrected_peak
            ) for i in 1:num_bins
        ]

        Projection(spectrum, c.detpair.name, ecal, nothing)
    else
        if length(bins.bins) < 2
            throw(ErrorException("Invalid bins for projection"))
        end

        bins_kev = [to_kev(x, c.detpair.eres) for x in bins.bins]

        spectrum = zeros(Float64, length(bins_kev) - 1)

        for i in 1:length(bins_kev) - 1
            left_edge = bins_kev[i]
            right_edge = bins_kev[i + 1]

            width_cel = right_edge - left_edge
            offset = 0.5 * (left_edge + right_edge)

            spectrum[i] = integrate(
                c,
                DiagonalAreaWithOffset(
                    KeV(width_cel), width, KeV(offset), KeV(0)
                ),
                in_corrected_peak
            )
        end

        Projection(spectrum, c.detpair.name, nothing, bins_kev)
    end
end
