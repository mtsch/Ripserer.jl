"""

"""
struct ImageBoundaryMatrix{
    T<:Number,
    F,
    S1<:Simplex,
    S2<:Simplex,
    R<:ReducedMatrix,
    B<:Chain{T},
    C<:Chain{T},
    I<:Integer
}
    filtration::F
    reduced::R
    buffer::B
    chain::C
    birth_candidates::Vector{S1}
    columns_to_reduce::Vector{S2}
    zeroed::Set{S2}
    subset::Int
    birth_intervals::Dict{I,PersistenceInterval}
    infinite_intervals::Bool
end

function ImageBoundaryMatrix(
    ::Type{T}, filtration, birth_candidates, columns_to_reduce, subset, intervals;
    infinite_intervals=true,
) where {T}
    S2 = eltype(columns_to_reduce)
    I = index_type(S2)
    B = birth_type(S2)
    S1 = Simplex{dim(S2) - 1, Tuple{Bool,B}, I}

    reduced = ReducedMatrix{S1,T,S1}()
    sizehint!(reduced, length(columns_to_reduce))

    buffer = Chain{T,S1}()
    chain = Chain{T,S1}()

    birth_intervals = Dict(index(birth_simplex(int)) => int for int in intervals)

    return ImageBoundaryMatrix{
        T,typeof(filtration),S1,S2,typeof(reduced),typeof(buffer),typeof(chain),I
    }(
        filtration, reduced, buffer, chain,
        birth_candidates, columns_to_reduce, Set{S2}(), subset, birth_intervals,
        infinite_intervals,
    )
end

field_type(::ImageBoundaryMatrix{T}) where {T} = T
dim(::ImageBoundaryMatrix{<:Any,<:Any,S}) where {S} = dim(S)
ordering(::ImageBoundaryMatrix) = Base.Order.Reverse

is_implicit(::ImageBoundaryMatrix) = false
is_cohomology(::ImageBoundaryMatrix) = false

function attempt_early_stop!(matrix::ImageBoundaryMatrix, i, columns)
    if length(matrix.reduced.column_index) ≥ length(matrix.birth_candidates)
        # At this point, all potential births have been found. The rest of the columns
        # should be marked as zeroed.
        for j in i:length(columns)
            push!(matrix.zeroed, columns[j])
        end
        return true
    else
        return false
    end
end

function coboundary(matrix::ImageBoundaryMatrix, simplex::AbstractCell)
    return ImageBoundary(boundary(matrix.filtration, simplex), matrix.subset)
end

function append_infinite_intervals!(intervals, matrix::ImageBoundaryMatrix)
    if matrix.infinite_intervals && length(matrix.birth_candidates) ≠ length(intervals)
        for simplex in matrix.birth_candidates
            if !haskey(matrix.reduced, simplex)
                connection = get(matrix.birth_intervals, index(simplex), nothing)
                if isnothing(connection)
                    continue
                end
                push!(
                    intervals,
                    PersistenceInterval(
                        birth(simplex)[2],
                        Inf,
                        (; birth_simplex=simplex, death_simplex=nothing, mixup=connection),
                    ),
                )
            end
        end
    end
    return intervals
end

function mark_zero_column!(matrix::ImageBoundaryMatrix, column_index)
    if matrix.infinite_intervals
        push!(matrix.zeroed, column_index)
    end
    return nothing
end

function next_matrix(matrix::ImageBoundaryMatrix, intervals)
    birth_candidates = filter(matrix.columns_to_reduce) do sx
        sx in matrix.zeroed
    end
    birth_candidates = map(σ -> _add_image_birth(σ, matrix.subset), birth_candidates)
    columns = simplex_type(matrix.filtration, dim(matrix) + 2)[]
    for col in columns_to_reduce(matrix.filtration, matrix.columns_to_reduce)
        push!(columns, abs(col))
    end

    return ImageBoundaryMatrix(
        field_type(matrix), matrix.filtration, birth_candidates, columns,
        matrix.subset, intervals
    )
end

"""
    ImageBoundary

Wraps over `Boundary` and adds a flag to the births indicating whether the simplex belongs
to the subset or not.
"""
struct ImageBoundary{D,I,B<:Boundary{D,I}}
    boundary::B
    subset::Int
end

function Base.iterate(ib::ImageBoundary{D}, k=1) where {D}
    next = iterate(ib.boundary, k)
    if !isnothing(next)
        sx, k = next
        modified_sx = _add_image_birth(sx, ib.subset)
        return modified_sx, k
    else
        return nothing
    end
end

"""
    interval(matrix, column, pivot, cutoff, reps)

Construct a persistence interval.
"""
function interval(matrix::ImageBoundaryMatrix, column, pivot, cutoff, reps)
    if isnothing(pivot)
        # In homology, birth simplex is nothing when column is fully reduced.
        return nothing
    else
        birth_simplex, death_simplex = simplex(pivot), column
    end

    birth_time = Float64(birth(birth_simplex)[2])
    death_time = isnothing(death_simplex) ? Inf : Float64(birth(death_simplex))

    # Image interval only if birth simplex is also a birth simplex in the subfiltration
    # Record the intervals they are connected to.
    connection = get(matrix.birth_intervals, index(birth_simplex), nothing)
    if isnothing(connection)
        return nothing
    end

    if reps
        rep = (; representative=collect_cocycle!(matrix, column, pivot))
    else
        rep = NamedTuple()
    end
    meta = (
        ;
        birth_simplex=birth_simplex,
        death_simplex=death_simplex,
        mixup=connection,
        rep...,
    )
    return PersistenceInterval(birth_time, death_time, meta)
end


###
### Actual computation
###
export ImageHomology

struct ImageHomology
    subset::Int
end

function _sx_vertices(filtration, subset)
    vs = vertices(filtration)
    if subset > 0
        return collect(
            _add_image_birth(simplex(filtration, Val(0), (i,)), subset) for i in vs
        )
    else
        return collect(simplex(filtration, Val(0), (i,)) for i in vs)
    end
end

function _ripserer(
    ip::ImageHomology, filtration, cutoff, verbose, field, dim_max, reps, implicit
)
    sub_result = PersistenceDiagram[]
    img_result = PersistenceDiagram[]

    @prog_println verbose "Subset:"
    subfiltration = subset(filtration, ip.subset)

    sub_edges = edges(subfiltration)
    sub_matrix = BoundaryMatrix{false}(
        field, subfiltration, _sx_vertices(subfiltration, 0), edges(subfiltration)
    )
    sub_intervals = compute_intervals!(sub_matrix, 0, verbose, _reps(reps, dim))

    # Image
    @prog_println verbose "Image:"
    img_matrix = ImageBoundaryMatrix(
        field, filtration, _sx_vertices(filtration, ip.subset), edges(filtration),
        ip.subset, sub_intervals,
    )
    img_intervals = compute_intervals!(img_matrix, cutoff, verbose, _reps(reps, dim))

    push!(sub_result, sub_intervals)
    push!(img_result, img_intervals)

    for dim in 1:dim_max
        @prog_println verbose "Subset:"
        sub_matrix = next_matrix(sub_matrix)
        sub_intervals = compute_intervals!(sub_matrix, cutoff, verbose, _reps(reps, dim))

        @prog_println verbose "Image:"
        img_matrix = next_matrix(img_matrix, sub_intervals)
        img_intervals = compute_intervals!(img_matrix, cutoff, verbose, _reps(reps, dim))

        useful_cols = setdiff(img_matrix.columns_to_reduce, img_matrix.zeroed)
        @show useful_cols

        push!(sub_result, sub_intervals)
        push!(img_result, img_intervals)
    end
    return sub_result, img_result
end
