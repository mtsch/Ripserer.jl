export ImageCohomology
struct ImageCohomology
    subset::Int
end

export ImageCohomology2
struct ImageCohomology2
    subset::Int
end

struct ImageCoboundaryMatrix{
    Flip,T<:Number,F,S<:AbstractCell,R<:ReducedMatrix,B<:Chain{T},C<:Chain{T},I
} <: AbstractCoboundaryMatrix{true}
    filtration::F
    reduced::R
    buffer::B     # stores the columns that were added to current chain
    chain::C      # current column as it's being reduced
    columns_to_reduce::Vector{S}
    columns_to_skip::Vector{S}
    subset::Int
    birth_intervals::Dict{I,PersistenceInterval}
end
function ImageCoboundaryMatrix{Flip,T,F,S,R,B,C,I}(
    ::Type{T}, filtration::F, columns_to_reduce, columns_to_skip, subset, intervals, flip
) where {Flip,T,F,S,R,B,C,I}
    return ImageCoboundaryMatrix(
        T, filtration, columns_to_reduce, columns_to_skip, subset, intervals, flip
    )
end
function ImageCoboundaryMatrix(
    ::Type{T}, filtration, columns_to_reduce, columns_to_skip, subset, intervals, flip
) where {T}
    S = eltype(columns_to_reduce)
    I = index_type(S)
    B = birth_type(S)
    if flip
        C = Simplex{dim(S) + 1, Tuple{Bool,B}, I}
    else
        C = simplex_type(filtration, dim(S) + 1)
    end

    reduced = ReducedMatrix{C,T,S}()
    buffer = Chain{T,S}()

    sizehint!(reduced, length(columns_to_reduce))
    chain = Chain{T,C}()

    birth_intervals = Dict(index(birth_simplex(int)) => int for int in intervals)

    return ImageCoboundaryMatrix{
        flip,T,typeof(filtration),S,typeof(reduced),typeof(buffer),typeof(chain),I
    }(
        filtration, reduced, buffer, chain, columns_to_reduce, columns_to_skip, subset, birth_intervals
    )
end

dim(::ImageCoboundaryMatrix{<:Any,<:Any,<:Any,S}) where {S} = dim(S)
field_type(::ImageCoboundaryMatrix{<:Any,F}) where {F} = F
is_flipped(::ImageCoboundaryMatrix{F}) where {F} = F

function sort_columns!(matrix::ImageCoboundaryMatrix{false})
    sort!(matrix.columns_to_reduce; rev=true, by=x -> _add_image_birth(x, matrix.subset))
end
function sort_columns!(matrix::ImageCoboundaryMatrix{true})
    sort!(matrix.columns_to_reduce; rev=true)
end

function interval(matrix::ImageCoboundaryMatrix{false}, column, pivot, cutoff, reps)
    if reps
        representative = collect_cocycle!(matrix, column, pivot)
    else
        representative = nothing
    end
    if !isnothing(pivot)
        return _image_interval(matrix, column, simplex(pivot), cutoff, representative)
    else
        return _image_interval(matrix, column, pivot, cutoff, representative)
    end
end
interval(::ImageCoboundaryMatrix{true}, _, _, _, _) = nothing

struct FlippedImageCoboundary{C<:Coboundary}
    coboundary::C
    subset::Int
end
function Base.iterate(ic::FlippedImageCoboundary, args...)
    next = iterate(ic.coboundary, args...)
    if !isnothing(next)
        sx, st = next
        return _add_image_birth(sx, ic.subset), st
    else
        return nothing
    end
end
function coboundary(matrix::ImageCoboundaryMatrix{false}, simplex::Simplex)
    return coboundary(matrix.filtration, simplex)
end
function coboundary(matrix::ImageCoboundaryMatrix{true}, simplex::Simplex)
    return FlippedImageCoboundary(coboundary(matrix.filtration, simplex), matrix.subset)
end

function next_matrix(matrix::ImageCoboundaryMatrix{F}, verbose, intervals, flip) where {F}
    C = simplex_type(matrix.filtration, dim(matrix) + 1)
    new_to_reduce = C[]
    new_to_skip = C[]
    sizehint!(new_to_skip, length(matrix.reduced))

    if verbose
        progbar = ProgressUnknown("Assembling columns:")
    end
    for simplex in columns_to_reduce(
        matrix.filtration,
        Iterators.flatten((matrix.columns_to_reduce, matrix.columns_to_skip)),
    )
        if F
            sx = _add_image_birth(simplex, matrix.subset)
        else
            sx = simplex
        end
        if haskey(matrix.reduced, sx)
            push!(new_to_skip, abs(simplex))
        else
            push!(new_to_reduce, abs(simplex))
        end
        verbose && next!(
            progbar;
            showvalues=(
                ("cleared", length(new_to_skip)), ("to reduce", length(new_to_reduce))
            ),
        )
    end
    verbose && finish!(progbar)

    if flip
        ImageCoboundaryMatrix(
            field_type(matrix), matrix.filtration, new_to_reduce, new_to_skip, matrix.subset, intervals, !F
        )
    else
        ImageCoboundaryMatrix(
            field_type(matrix), matrix.filtration, new_to_reduce, new_to_skip, matrix.subset, intervals, F
        )
    end
end

###########################################################################################
function _ripserer(
    ip::ImageCohomology, filtration, cutoff, verbose, field, dim_max, reps, implicit
)
    img_result = PersistenceDiagram[]

    @prog_println verbose "Subset:"
    subfiltration = subset(filtration, ip.subset)

    sub_edges = edges(subfiltration)
    sub_vertices = _sx_vertices(subfiltration, 0)
    sub_matrix = CoboundaryMatrix{true}(
        field, subfiltration, sub_vertices, empty(sub_vertices)
    )
    sub_intervals = compute_intervals!(sub_matrix, cutoff, verbose, _reps(reps, dim))

    # Image
    @prog_println verbose "Image:"
    img_vertices = _sx_vertices(filtration, 0)
    img_matrix = ImageCoboundaryMatrix(
        field, filtration, img_vertices, empty(img_vertices), ip.subset, sub_intervals, false
    )
    img_intervals = compute_intervals!(img_matrix, cutoff, verbose, _reps(reps, dim))

    push!(img_result, img_intervals)

    for dim in 1:dim_max
        @prog_println verbose "Subset:"
        sub_matrix = next_matrix(sub_matrix, verbose)
        sub_intervals = compute_intervals!(sub_matrix, cutoff, verbose, _reps(reps, dim))

        @prog_println verbose "Image:"
        img_matrix = next_matrix(img_matrix, verbose, ip.subset, sub_intervals, false)
        # No clearing
        append!(img_matrix.columns_to_reduce, img_matrix.columns_to_skip)
        img_intervals = compute_intervals!(img_matrix, cutoff, verbose, _reps(reps, dim))

        push!(img_result, img_intervals)
    end
    return img_result
end

function _ripserer(
    ip::ImageCohomology2, filtration, cutoff, verbose, field, dim_max, reps, implicit
)
    img_result = PersistenceDiagram[]

    @prog_println verbose "Subset:"
    subfiltration = subset(filtration, ip.subset)

    sub_edges = edges(subfiltration)
    sub_vertices = _sx_vertices(subfiltration, 0)
    sub_matrix = CoboundaryMatrix{true}(
        field, subfiltration, sub_vertices, empty(sub_vertices)
    )
    sub_intervals = compute_intervals!(
        sub_matrix, cutoff, verbose, _reps(reps, dim)
    )

    # Image
    @prog_println verbose "Image A:"
    to_reduce = _sx_vertices(filtration, 0)
    to_skip = empty(to_reduce)

    img_matrix_A = ImageCoboundaryMatrix(
        field, filtration, to_reduce, to_skip, ip.subset, sub_intervals, false
    )
    img_intervals_A = compute_intervals!(img_matrix_A, cutoff, verbose, _reps(reps, dim))

    @prog_println verbose "Image B:"
    img_matrix_B = ImageCoboundaryMatrix(
        field, filtration, to_reduce, to_skip, ip.subset, sub_intervals, true
    )
    img_intervals_B = compute_intervals!(img_matrix_B, cutoff, verbose, _reps(reps, dim))

    push!(img_result, img_intervals_A)

    for dim in 1:dim_max
        @prog_println verbose "Subset:"
        sub_matrix = next_matrix(sub_matrix, verbose)
        sub_intervals = compute_intervals!(
            sub_matrix, cutoff, verbose, _reps(reps, dim); sort_columns=true,
        )

        if dim != dim_max || iseven(dim)
            @prog_println verbose "Image A:"
            img_matrix_A = next_matrix(img_matrix_A, verbose, sub_intervals, true)
            img_intervals_A = compute_intervals!(
                img_matrix_A, cutoff, verbose, _reps(reps, dim),
            )
        end

        if dim ≠ dim_max || isodd(dim)
            @prog_println verbose "Image B:"
            img_matrix_B = next_matrix(img_matrix_B, verbose, sub_intervals, true)
            img_intervals_B = compute_intervals!(
                img_matrix_B, cutoff, verbose, _reps(reps, dim),
            )
        end

        if iseven(dim)
            push!(img_result, img_intervals_A)
        else
            push!(img_result, img_intervals_B)
        end
    end
    return img_result
end
