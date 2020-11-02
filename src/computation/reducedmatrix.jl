"""
    ReducedMatrix{I<:AbstractSimplex, C<:Chain}

Representation of the reduced part of the matrix. Is indexed by `I`, indexing returns views
into a `Chain` of type `C`.

!!! warning
    This type does no input validation. Doing things like adding columns that already exist
    is undefined.
"""
struct ReducedMatrix{I<:AbstractSimplex,F,S,C<:Chain{F,S}}
    column_index::Dict{I,Int}
    indices::Vector{Int}
    values::C

    function ReducedMatrix{I,F,S}() where {I,F,S}
        column_index = Dict{S,Int}()
        indices = Int[1]
        values = Chain{F,S}()
        return new{I,F,S,typeof(values)}(column_index, indices, values)
    end
end

function Base.sizehint!(matrix::ReducedMatrix, size)
    sizehint!(matrix.column_index, size)
    sizehint!(matrix.indices, size + 1)
    sizehint!(matrix.values, size)
    return matrix
end

function Base.setindex!(matrix::ReducedMatrix, column::Chain, simplex)
    @assert sign(simplex) == 1
    n = length(column)
    if n > 0
        append!(matrix.values, column)
        push!(matrix.indices, matrix.indices[end] + n)
        matrix.column_index[simplex] = length(matrix.indices) - 1
    end
    return column
end

Base.length(matrix::ReducedMatrix) = length(matrix.indices) - 1
function Base.haskey(matrix::ReducedMatrix{I}, simplex::I) where {I}
    return haskey(matrix.column_index, simplex)
end
function Base.haskey(matrix::ReducedMatrix{I}, element::AbstractChainElement{I}) where {I}
    return haskey(matrix.column_index, simplex(element))
end

function Base.getindex(matrix::ReducedMatrix{S}, element::AbstractChainElement{S}) where {S}
    return matrix[simplex(element)]
end
function Base.getindex(matrix::ReducedMatrix{S}, simplex::S) where {S}
    index = get(matrix.column_index, simplex, 0)
    if index == 0
        from = 0
        to = -1
    else
        from = matrix.indices[index]
        to = matrix.indices[index + 1] - 1
    end
    # Constructing directly is equivalent to Base.unsafe_view
    return SubArray(matrix.values, (from:to,))
end

"""
    bulk_insert!(matrix::ReducedMatrix, pairs)

Insert apparent pairs into reduced matrix. Columns are inserted so `matrix[τ] == [σ]`, where
`pairs` is a collection of tuples `(σ, τ)`.
"""
function bulk_insert!(matrix::ReducedMatrix{<:Any,F,S}, pairs) where {F,S}
    n = length(pairs)
    n_values = length(matrix.values)
    resize!(matrix.values, n_values + n)
    n_indices = length(matrix.indices)
    last_index = last(matrix.indices)
    resize!(matrix.indices, n_indices + n)
    @inbounds for (i, (σ, τ)) in enumerate(pairs)
        matrix.indices[n_indices + i] = last_index + i
        matrix.values[n_values + i] = ChainElement{S,F}(σ, sign(τ))
        matrix.column_index[abs(τ)] = n_indices + i - 1
    end
    return matrix
end
function bulk_insert!(matrix::ReducedMatrix, ::Tuple{})
    return matrix
end
