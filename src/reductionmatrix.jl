"""
    ReductionMatrix{S<:AbstractSimplex, E<:AbstractChainElement}

Representation of the reduced part of the matrix. Is indexed by `S`, indexing iterates
values of type `E`. `is_homology` controls the ordering of values.

Changes to the matrix are recorded in a buffer with `record!`. Once a column is ready to be
added, use `commit!` to move the changes from the buffer to the matrix and create a
column. Changes can also be discarded with `undo!`.
"""
struct ReductionMatrix{S<:AbstractSimplex, E<:AbstractChainElement}
    column_index::Dict{S, Int}
    indices::Vector{Int}
    values::Vector{E}
    buffer::Vector{E}
    is_homology::Bool

    function ReductionMatrix{S, E}(is_homology) where {S, E}
        return new{S, E}(Dict{S, Int}(), Int[1], E[], E[], is_homology)
    end
end

"""
    record!(matrix::ReductionMatrix, element)

Record the operation that was performed in the buffer.
"""
function record!(matrix::ReductionMatrix, element)
    push!(matrix.buffer, element)
    return element
end

"""
    commit!(matrix::ReductionMatrix, simplex, factor)

Commit the changes that were `record!`ed, creating a new column indexed by `simplex`. This
sorts the buffer and adds duplicates together. All entries are multiplied by `factor`.
"""
function commit!(matrix::ReductionMatrix, simplex, factor)
    isempty(matrix.buffer) && return matrix

    sort!(matrix.buffer, alg=QuickSort, rev=matrix.is_homology)
    i = 0
    @inbounds prev = matrix.buffer[1]
    @inbounds for j in 2:length(matrix.buffer)
        current = matrix.buffer[j]
        if current == prev
            prev += current
        else
            if !iszero(prev)
                i += 1
                matrix.buffer[i] = prev * factor
            end
            prev = current
        end
    end
    if !iszero(prev)
        i += 1
        @inbounds matrix.buffer[i] = prev * factor
    end
    if i > 0
        resize!(matrix.buffer, i)
        append!(matrix.values, matrix.buffer)
        push!(matrix.indices, i + 1)
        matrix.column_index[abs(simplex)] = length(matrix.indices) - 1
    end

    empty!(matrix.buffer)
    return matrix
end

"""
    undo!(matrix::ReductionMatrix)

Undo recorded changes.
"""
function undo!(matrix::ReductionMatrix)
    empty!(matrix.buffer)
    return matrix
end

Base.eltype(::Type{<:ReductionMatrix{<:Any, E}}) where E = E
Base.length(matrix::ReductionMatrix) = length(matrix.indices) - 1

function Base.getindex(matrix::ReductionMatrix{C}, element::AbstractChainElement{C}) where C
    index = get(matrix.column_index, simplex(element), 0)
    return RMColumnIterator{typeof(matrix)}(matrix, index)
end

function Base.getindex(matrix::ReductionMatrix{C}, simplex::C) where C
    index = get(matrix.column_index, simplex, 0)
    return RMColumnIterator{typeof(matrix)}(matrix, index)
end

"""
    RMColumnIterator{R}

An iterator over a column of a `ReductionMatrix`.
"""
struct RMColumnIterator{R<:ReductionMatrix}
    matrix::R
    index::Int
end

Base.IteratorSize(::Type{RMColumnIterator}) = Base.HasLength()
Base.IteratorEltype(::Type{RMColumnIterator}) = Base.HasEltype()
Base.eltype(::Type{<:RMColumnIterator{R}}) where R = eltype(R)
function Base.length(iter::RMColumnIterator)
    index = iter.index
    if index == 0
        return 0
    else
        return iter.matrix.indices[index + 1] - iter.matrix.indices[index]
    end
end

function Base.iterate(iter::RMColumnIterator)
    if iter.index == 0
        return nothing
    else
        return iterate(iter, 1)
    end
end

function Base.iterate(iter::RMColumnIterator, i)
    indices = iter.matrix.indices
    index = i + indices[iter.index] - 1
    if index ≥ indices[iter.index + 1]
        return nothing
    else
        return (iter.matrix.values[index], i + 1)
    end
end
