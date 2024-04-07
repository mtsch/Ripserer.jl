"""
Write a blurb here.
"""
struct Homology
    implicit::Bool

    Homology(; implicit=true) = new(implicit)
end

"""
    BoundaryMatrix{I}

This `struct` is used to compute homology. The `I` parameter sets whether the implicit
algoritm is used or not.
"""
struct BoundaryMatrix{
    I,T<:Number,F,S1<:AbstractCell,S2<:AbstractCell,R<:ReducedMatrix,B<:Chain{T},C<:Chain{T}
}
    filtration::F
    reduced::R
    buffer::B
    chain::C
    birth_candidates::Vector{S1} # used to find infinite intervals
    columns_to_reduce::Vector{S2}
    zeroed::Set{S2}
    infinite_intervals::Bool
end

function BoundaryMatrix{I}(
    ::Type{T}, filtration, birth_candidates, columns_to_reduce; infinite_intervals=true
) where {I,T}
    S2 = eltype(columns_to_reduce)
    S1 = simplex_type(filtration, dim(S2) - 1)

    if !I
        reduced = ReducedMatrix{S1,T,S1}()
        buffer = Chain{T,S1}()
    else
        reduced = ReducedMatrix{S1,T,S2}()
        buffer = Chain{T,S2}()
    end
    sizehint!(reduced, length(columns_to_reduce))
    chain = Chain{T,S1}()

    return BoundaryMatrix{
        I,T,typeof(filtration),S1,S2,typeof(reduced),typeof(buffer),typeof(chain)
    }(
        filtration,
        reduced,
        buffer,
        chain,
        birth_candidates,
        columns_to_reduce,
        Set{S2}(),
        infinite_intervals,
    )
end

field_type(::BoundaryMatrix{<:Any,T}) where {T} = T
dim(::BoundaryMatrix{<:Any,<:Any,<:Any,S}) where {S} = dim(S)
ordering(::BoundaryMatrix) = Base.Order.Reverse

is_implicit(::BoundaryMatrix{I}) where {I} = I
is_cohomology(::BoundaryMatrix) = false

function sort_columns!(matrix::BoundaryMatrix)
    sort!(matrix.columns_to_reduce)
end

function attempt_early_stop!(matrix::BoundaryMatrix, i, columns)
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

# The naming here is not ideal...
function coboundary(matrix::BoundaryMatrix, simplex::AbstractCell)
    return boundary(matrix.filtration, simplex)
end

function mark_zero_column!(matrix::BoundaryMatrix, column_index)
    if matrix.infinite_intervals
        push!(matrix.zeroed, column_index)
    end
    return nothing
end

function append_infinite_intervals!(intervals, matrix::BoundaryMatrix)
    if matrix.infinite_intervals && length(matrix.birth_candidates) ≠ length(intervals)
        for simplex in matrix.birth_candidates
            if !haskey(matrix.reduced, simplex)
                push!(
                    intervals,
                    PersistenceInterval(
                        birth(simplex),
                        Inf,
                        (; birth_simplex=simplex, death_simplex=nothing),
                    ),
                )
            end
        end
    end
    return intervals
end

function next_matrix(matrix::BoundaryMatrix{I}) where {I}
    birth_candidates = filter(matrix.columns_to_reduce) do sx
        sx in matrix.zeroed
    end
    columns = simplex_type(matrix.filtration, dim(matrix) + 2)[]
    for col in columns_to_reduce(matrix.filtration, matrix.columns_to_reduce)
        push!(columns, abs(col))
    end

    return BoundaryMatrix{I}(
        field_type(matrix), matrix.filtration, birth_candidates, columns
    )
end

function _ripserer(
    ::Val{:homology}, filtration, cutoff, verbose, field, dim_max, reps, implicit
)
    result = PersistenceDiagram[]
    zeroth, to_reduce, to_skip = zeroth_intervals(
        filtration, cutoff, verbose, field, _reps(reps, 0)
    )
    push!(result, zeroth)
    if dim_max > 0
        birth_candidates = to_reduce
        columns = simplex_type(filtration, 2)[]
        for col in columns_to_reduce(filtration, Iterators.flatten((to_reduce, to_skip)))
            push!(columns, abs(col))
        end
        matrix = BoundaryMatrix{implicit}(field, filtration, birth_candidates, columns)
        for dim in 1:dim_max
            push!(
                result,
                compute_intervals!(
                    matrix, cutoff, verbose, _reps(reps, dim) ;sort_columns=true
                ),
            )
            if dim < dim_max
                matrix = next_matrix(matrix)
            end
        end
    end
    return result
end
