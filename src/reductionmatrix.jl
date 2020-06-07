"""
    ReducedMatrix{S<:AbstractSimplex, E<:AbstractChainElement}

Representation of the reduced part of the matrix. Is indexed by `S`, indexing iterates
values of type `E`. `is_homology` controls the ordering of values.

Changes to the matrix are recorded in a buffer with `record!`. Once a column is ready to be
added, use `commit!` to move the changes from the buffer to the matrix and create a
column. Changes can also be discarded with `discard!`.
"""
struct ReducedMatrix{S<:AbstractSimplex, E<:AbstractChainElement, O<:Base.Ordering}
    column_index::Dict{S, Int}
    indices::Vector{Int}
    values::Vector{E}
    buffer::Vector{E}
    ordering::O

    function ReducedMatrix{S, E}(is_cohomology) where {S, E}
        ordering = is_cohomology ? Base.Order.Forward : Base.Order.Reverse
        return new{S, E, typeof(ordering)}(Dict{S, Int}(), Int[1], E[], E[], ordering)
    end
end

"""
    record!(matrix::ReducedMatrix, element)

Record the operation that was performed in the buffer.
"""
function record!(matrix::ReducedMatrix, element)
    push!(matrix.buffer, element)
    return element
end

"""
    commit!(matrix::ReducedMatrix, simplex, factor)

Commit the changes that were `record!`ed, creating a new column indexed by `simplex`. This
sorts the buffer and adds duplicates together. All entries are multiplied by `factor`.
"""
function commit!(matrix::ReducedMatrix, simplex, factor)
    isempty(matrix.buffer) && return matrix

    sort!(matrix.buffer, alg=QuickSort, order=matrix.ordering)
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
        push!(matrix.indices, matrix.indices[end] + i)
        matrix.column_index[abs(simplex)] = length(matrix.indices) - 1
    end

    empty!(matrix.buffer)
    return matrix
end

"""
    discard!(matrix::ReducedMatrix)

Undo recorded changes.
"""
function discard!(matrix::ReducedMatrix)
    empty!(matrix.buffer)
    return matrix
end

Base.eltype(::Type{<:ReducedMatrix{<:Any, E}}) where E = E
Base.length(matrix::ReducedMatrix) = length(matrix.indices) - 1

function Base.getindex(matrix::ReducedMatrix{C}, element::AbstractChainElement{C}) where C
    matrix[simplex(element)]
end

function Base.getindex(matrix::ReducedMatrix{C}, simplex::C) where C
    index = get(matrix.column_index, abs(simplex), 0)
    return RMColumnIterator{typeof(matrix)}(matrix, index)
end

"""
    RMColumnIterator{R}

An iterator over a column of a `ReducedMatrix`.
"""
struct RMColumnIterator{R<:ReducedMatrix}
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

# ======================================================================================== #
struct WorkingBoundary{E<:AbstractChainElement, O<:Base.Ordering}
    heap::Vector{E}
    ordering::O

    function WorkingBoundary{E}(is_cohomology) where E
        ordering = is_cohomology ? Base.Order.Forward : Base.Order.Reverse
        new{E, typeof(ordering)}(E[], ordering)
    end
end


Base.empty!(col::WorkingBoundary) = empty!(col.heap)
Base.isempty(col::WorkingBoundary) = isempty(col.heap)

function Base.sizehint!(col::WorkingBoundary, size)
    sizehint!(col.heap, size)
    return col
end

function _pop_pivot!(column::WorkingBoundary)
    isempty(column) && return nothing
    heap = column.heap

    pivot = heappop!(heap, column.ordering)
    while !isempty(heap)
        if iszero(pivot)
            pivot = heappop!(heap, column.ordering)
        elseif first(heap) == pivot
            pivot += heappop!(heap, column.ordering)
        else
            break
        end
    end
    return iszero(pivot) ? nothing : pivot
end

"""
    get_pivot!(column)

Return the pivot of the column - the element with the lowest diameter.
"""
function get_pivot!(column::WorkingBoundary)
    pivot = _pop_pivot!(column)
    if !isnothing(pivot)
        heappush!(column.heap, pivot, column.ordering)
    end
    return pivot
end

function Base.push!(column::WorkingBoundary{E}, simplex::AbstractSimplex) where E
    push!(column, E(simplex))
end
function Base.push!(column::WorkingBoundary{E}, element::E) where E
    heap = column.heap
    @inbounds if !isempty(heap) && heap[1] == element
        heap[1] += element
        if iszero(heap[1])
            heappop!(heap, column.ordering)
        end
    else
        heappush!(heap, element, column.ordering)
    end
    return column
end

function nonheap_push!(column::WorkingBoundary{E}, simplex::AbstractSimplex) where E
    push!(column.heap, E(simplex))
end
function nonheap_push!(column::WorkingBoundary{E}, element::E) where E
    push!(column.heap, element)
end

repair!(column::WorkingBoundary) = heapify!(column.heap, column.ordering)
Base.first(column::WorkingBoundary) = isempty(column) ? nothing : first(column.heap)

function move!(column::WorkingBoundary{E}) where E
    dst = E[]
    while (pivot = _pop_pivot!(column)) ≠ nothing
        push!(dst, pivot)
    end
    return dst
end

# ======================================================================================== #
struct ReductionMatrix{
    Co, Field, Filtration, Simplex, SimplexElem, Face, FaceElem, O<:Base.Ordering
}
    filtration::Filtration
    reduced::ReducedMatrix{Face, SimplexElem, O}
    working_boundary::WorkingBoundary{FaceElem, O}
    columns_to_reduce::Vector{Simplex}
    columns_to_skip::Vector{Simplex}
end

function ReductionMatrix{Co, Field}(
    filtration::Filtration, columns_to_reduce, columns_to_skip
) where {Co, Field, Filtration}

    Simplex = eltype(columns_to_reduce)
    Face = Co ? coface_type(Simplex) : face_type(Simplex)
    O = Co ? typeof(Base.Order.Forward) : typeof(Base.Order.Reverse)
    SimplexElem = chain_element_type(Simplex, Field)
    FaceElem = chain_element_type(Face, Field)

    reduced = ReducedMatrix{Face, SimplexElem}(Co)
    working_boundary = WorkingBoundary{FaceElem}(Co)

    return ReductionMatrix{Co, Field, Filtration, Simplex, SimplexElem, Face, FaceElem, O}(
        filtration,
        reduced,
        working_boundary,
        columns_to_reduce,
        columns_to_skip
    )
end

"""
    co_boundary(matrix, simplex)

Iterate over the (co)boundary of the `simplex`. Chooses between the boundary and the
coboundary based on the matrix type.
"""
function co_boundary(matrix::ReductionMatrix{true}, simplex::AbstractSimplex)
    return coboundary(matrix.filtration, simplex)
end
function co_boundary(matrix::ReductionMatrix{false}, simplex::AbstractSimplex)
    return boundary(matrix.filtration, simplex)
end

simplex_type(::ReductionMatrix{<:Any, <:Any, <:Any, S}) where S = S
simplex_element(::ReductionMatrix{<:Any, <:Any, <:Any, <:Any, E}) where E = E
face_element(::ReductionMatrix{<:Any, <:Any, <:Any, <:Any, <:Any, <:Any, E}) where E = E
dim(::ReductionMatrix{true, <:Any, <:Any, S}) where S = dim(S)
dim(::ReductionMatrix{false, <:Any, <:Any, S}) where S = dim(S) - 1

function initialize_boundary!(matrix::ReductionMatrix{Co}, column_simplex) where Co
    emergent_check = matrix.filtration isa AbstractFlagFiltration
    empty!(matrix.working_boundary)
    for face in co_boundary(matrix, column_simplex)
        if emergent_check && diam(face) == diam(column_simplex)
            if isempty(matrix.reduced[face])
                empty!(matrix.working_boundary)
                return face_element(matrix)(face)
            end
            emergent_check = false
        end
        nonheap_push!(matrix.working_boundary, face)
    end
    if isempty(matrix.working_boundary)
        return nothing
    else
        repair!(matrix.working_boundary)
        return first(matrix.working_boundary)
    end
end

function add!(matrix::ReductionMatrix, column, factor)
    for element in column
        record!(matrix.reduced, factor * element)
        for face in co_boundary(matrix, simplex(element))
            push!(
                matrix.working_boundary,
                face_element(matrix)(face, factor * coefficient(element))
            )
        end
    end
    return matrix
end

function reduce_column!(matrix::ReductionMatrix, column_simplex, cutoff, reps)
    pivot = initialize_boundary!(matrix, column_simplex)

    while !isnothing(pivot)
        column = matrix.reduced[pivot]
        isempty(column) && break

        add!(matrix, column, -coefficient(pivot))
        pivot = get_pivot!(matrix.working_boundary)
    end
    if isnothing(pivot)
        discard!(matrix.reduced)
    else
        record!(matrix.reduced, simplex_element(matrix)(column_simplex))
        isempty(matrix.reduced.buffer) && error()
        commit!(matrix.reduced, simplex(pivot), inv(coefficient(pivot)))
    end

    return interval(matrix, column_simplex, pivot, cutoff, reps)
end

function interval(matrix::ReductionMatrix{true}, column_simplex, pivot, cutoff, reps)
    birth = diam(column_simplex)
    if isnothing(pivot)
        death = Inf
    else
        death = diam(pivot)
    end

    if reps && !isfinite(death)
        return PersistenceInterval(birth, death, simplex_element(matrix)[])
    elseif reps && death - birth > cutoff
        representative = collect(matrix.reduced[pivot])
        return PersistenceInterval(birth, death, representative)
    elseif !isfinite(death) || death - birth > cutoff
        return PersistenceInterval(birth, death)
    else
        return nothing
    end
end

function interval(matrix::ReductionMatrix{false}, column_simplex, pivot, cutoff, reps)
    if isnothing(pivot)
        return nothing
    else
        birth = diam(pivot)
        death = diam(column_simplex)
    end

    if reps && !isfinite(death)
        return PersistenceInterval(birth, death, face_element(matrix)[])
    elseif reps && death - birth > cutoff
        representative = move!(matrix.working_boundary)
        return PersistenceInterval(birth, death, representative)
    elseif !isfinite(death) || death - birth > cutoff
        return PersistenceInterval(birth, death)
    else
        return nothing
    end
end

function compute_intervals!(matrix::ReductionMatrix{Co}, cutoff, reps, progress) where Co
    if reps
        if Co
            intervals = PersistenceInterval{Vector{simplex_element(matrix)}}[]
        else
            intervals = PersistenceInterval{Vector{face_element(matrix)}}[]
        end
    else
        intervals = PersistenceInterval{Nothing}[]
    end
    if progress
        progbar = Progress(
            length(matrix.columns_to_reduce), desc="Computing $(dim(matrix))d intervals... "
        )
    end
    for column in matrix.columns_to_reduce
        interval = reduce_column!(matrix, column, cutoff, reps)
        if !isnothing(interval)
            push!(intervals, interval)
        end
        progress && next!(progbar)
    end
    return sort!(
        PersistenceDiagram(dim(matrix), intervals, threshold(matrix.filtration))
    )
end

simplex_name(::Type{<:Simplex{2}}) = "triangles"
simplex_name(::Type{<:Simplex{3}}) = "tetrahedra"
simplex_name(::Type{<:AbstractSimplex{D}}) where D = "$D-simplices"

function next_matrix(matrix::ReductionMatrix{Co, Field}, progress) where {Co, Field}
    C = coface_type(simplex_type(matrix))
    new_to_reduce = C[]
    new_to_skip = C[]

    if progress
        progbar = Progress(
            length(matrix.columns_to_reduce) + length(matrix.columns_to_skip),
            desc="Assembling...             "
        )
    end
    for simplex in Iterators.flatten((matrix.columns_to_reduce, matrix.columns_to_skip))
        for coface in coboundary(matrix.filtration, simplex, Val(false))
            # Clearing optimization only enabled for cohomology.
            if Co && !isempty(matrix.reduced[coface])
                push!(new_to_skip, abs(coface))
            else
                push!(new_to_reduce, abs(coface))
            end
        end
        progress && next!(progbar)
    end
    sort!(new_to_reduce, rev=Co)
    if progress
        printstyled("Assembled $(length(new_to_reduce)) $(simplex_name(C)).\n", color=:green)
    end

    return ReductionMatrix{Co, Field}(matrix.filtration, new_to_reduce, new_to_skip)
end
