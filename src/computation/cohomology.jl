abstract type AbstractCoboundaryMatrix{I} end

is_implicit(::AbstractCoboundaryMatrix{I}) where {I} = I
is_cohomology(::AbstractCoboundaryMatrix) = true
ordering(::AbstractCoboundaryMatrix) = Base.Order.Forward

attempt_early_stop!(::AbstractCoboundaryMatrix, _, _) = false

function coboundary(matrix::AbstractCoboundaryMatrix, simplex::AbstractCell)
    return coboundary(matrix.filtration, simplex)
end

function next_matrix(matrix::AbstractCoboundaryMatrix{I}, verbose, args...) where {I}
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
        if haskey(matrix.reduced, simplex)
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

    return typeof(matrix)(
        field_type(matrix), matrix.filtration, new_to_reduce, new_to_skip, args...
    )
end

"""
    CoboundaryMatrix{I}

This `struct` is used to compute cohomology. The `I` parameter sets whether the implicit
algoritm is used or not.
"""
struct CoboundaryMatrix{
    I,T<:Number,F,S<:AbstractCell,R<:ReducedMatrix,B<:Chain{T},C<:Chain{T}
} <: AbstractCoboundaryMatrix{I}
    filtration::F
    reduced::R
    buffer::B     # stores the columns that were added to current chain
    chain::C      # current column as it's being reduced
    columns_to_reduce::Vector{S}
    columns_to_skip::Vector{S}
end

function CoboundaryMatrix{I,T,F,S,R,B,C}(
    ::Type{T}, filtration::F, columns_to_reduce, columns_to_skip
) where {I,T,F,S,R,B,C}
    return CoboundaryMatrix{I}(T, filtration, columns_to_reduce, columns_to_skip)
end
function CoboundaryMatrix{I}(
    ::Type{T}, filtration, columns_to_reduce, columns_to_skip
) where {I,T}
    S = eltype(columns_to_reduce)
    C = simplex_type(filtration, dim(S) + 1)
    if I
        reduced = ReducedMatrix{C,T,S}()
        buffer = Chain{T,S}()
    else
        reduced = ReducedMatrix{C,T,C}()
        buffer = Chain{T,C}()
    end
    sizehint!(reduced, length(columns_to_reduce))
    chain = Chain{T,C}()

    return CoboundaryMatrix{
        I,T,typeof(filtration),S,typeof(reduced),typeof(buffer),typeof(chain)
    }(
        filtration, reduced, buffer, chain, columns_to_reduce, columns_to_skip
    )
end

field_type(::CoboundaryMatrix{<:Any,T}) where {T} = T
dim(cm::CoboundaryMatrix{<:Any,<:Any,<:Any,S}) where {S} = dim(S)

function sort_columns!(matrix::CoboundaryMatrix)
    sort!(matrix.columns_to_reduce; rev=true)
end

function _ripserer(
    ::Val{:cohomology}, filtration, cutoff, verbose, field, dim_max, reps, implicit
)
    result = PersistenceDiagram[]
    zeroth, to_reduce, to_skip = zeroth_intervals(
        filtration, cutoff, verbose, field, _reps(reps, 0)
    )
    push!(result, zeroth)
    if dim_max > 0
        matrix = CoboundaryMatrix{implicit}(field, filtration, to_reduce, to_skip)
        for dim in 1:dim_max
            push!(
                result,
                compute_intervals!(
                    matrix, cutoff, verbose, _reps(reps, dim); sort_columns=dim > 1
                ),
            )
            if dim < dim_max
                matrix = next_matrix(matrix, verbose)
            end
        end
    end
    return result
end
