"""
    ReducedMatrix{S<:AbstractSimplex, E<:AbstractChainElement}

Representation of the reduced part of the matrix. Is indexed by `S`, indexing iterates
values of type `E`. `is_homology` controls the ordering of values.

Changes to the matrix are recorded in a buffer with `record!`. Once a column is ready to be
added, use `commit!` to move the changes from the buffer to the matrix and create a
column. Changes can also be discarded with `discard!`.
"""
struct ReducedMatrix{S<:AbstractSimplex,E<:AbstractChainElement,O<:Base.Ordering}
    column_index::Dict{S,Int}
    indices::Vector{Int}
    values::Vector{E}
    buffer::Vector{E}
    ordering::O

    function ReducedMatrix{S,E}(ordering::O) where {S,E,O}
        return new{S,E,O}(Dict{S,Int}(), Int[1], E[], E[], ordering)
    end
end

function Base.sizehint!(matrix::ReducedMatrix, size)
    sizehint!(matrix.column_index, size)
    sizehint!(matrix.indices, size + 1)
    sizehint!(matrix.values, size)
    return matrix
end

"""
    record!(matrix::ReducedMatrix, element)
    record!(matrix::ReducedMatrix, elements, factor)
    record!(matrix::ReducedMatrix, chain, factor)

Record the operation that was performed in the buffer. If the second argument is a
`WorkingChain`, it is emptied into the buffer.
"""
function record!(matrix::ReducedMatrix{<:Any,E}, simplex::AbstractSimplex) where {E}
    return push!(matrix.buffer, E(simplex))
end

function record!(matrix::ReducedMatrix, elements, factor)
    i = length(matrix.buffer)
    resize!(matrix.buffer, length(matrix.buffer) + length(elements))
    @inbounds for element in elements
        i += 1
        matrix.buffer[i] = element * factor
    end
    return elements
end

function record!(matrix::ReducedMatrix, chain::WorkingChain)
    while (element = pop!(chain)) ≢ nothing
        push!(matrix.buffer, element)
    end
    return chain
end

"""
    collect_buffer!(matrix::ReducedMatrix{<:Any, E}, factor=one(E))

Sort buffer and add duplicates together. Return buffer.
"""
function collect_buffer!(matrix::ReducedMatrix{<:Any,E}, factor=one(E)) where {E}
    if !isempty(matrix.buffer)
        sort!(matrix.buffer; alg=QuickSort, order=matrix.ordering)
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
        resize!(matrix.buffer, i)
    end
    return matrix.buffer
end

"""
    commit!(matrix::ReducedMatrix, simplex, factor)

Commit the changes that were `record!`ed, creating a new column indexed by `simplex`. This
sorts the buffer and adds duplicates together. All entries are multiplied by `factor`.
"""
function commit!(matrix::ReducedMatrix, simplex, factor)
    @assert sign(simplex) == 1
    collect_buffer!(matrix, factor)

    n = length(matrix.buffer)
    if n > 0
        append!(matrix.values, matrix.buffer)
        push!(matrix.indices, matrix.indices[end] + n)
        matrix.column_index[simplex] = length(matrix.indices) - 1
    end

    return matrix
end

"""
    clear_buffer!(matrix::ReducedMatrix)

Empty buffer.
"""
function clear_buffer!(matrix::ReducedMatrix)
    empty!(matrix.buffer)
    return matrix
end

"""
    bulk_add!(matrix::ReducedMatrix{<:Any, E}, pairs) where E

Insert apparent pairs into reduced matrix. Columns are inserted so `matrix[τ] == [σ]`, where
`pairs` is a collection of tuples `(σ, τ)`.
"""
function bulk_add!(matrix::ReducedMatrix{<:Any,E}, pairs) where {E}
    n = length(pairs)
    n_values = length(matrix.values)
    resize!(matrix.values, n_values + n)
    n_indices = length(matrix.indices)
    last_index = last(matrix.indices)
    resize!(matrix.indices, n_indices + n)
    @inbounds for (i, (σ, τ)) in enumerate(pairs)
        matrix.indices[n_indices + i] = last_index + i
        matrix.values[n_values + i] = E(σ, sign(τ))
        matrix.column_index[abs(τ)] = n_indices + i - 1
    end
    return matrix
end
function bulk_add!(matrix::ReducedMatrix, ::Tuple{})
    return matrix
end

Base.length(matrix::ReducedMatrix) = length(matrix.indices) - 1

Base.haskey(matrix::ReducedMatrix, simplex) = haskey(matrix.column_index, simplex)

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
