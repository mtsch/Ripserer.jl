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

    function ReducedMatrix{S, E}(ordering::O) where {S, E, O}
        return new{S, E, O}(Dict{S, Int}(), Int[1], E[], E[], ordering)
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

Record the operation that was performed in the buffer.
"""
function record!(matrix::ReducedMatrix, element::AbstractChainElement)
    push!(matrix.buffer, element)
    return element
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

"""
    commit!(matrix::ReducedMatrix, simplex, factor)

Commit the changes that were `record!`ed, creating a new column indexed by `simplex`. This
sorts the buffer and adds duplicates together. All entries are multiplied by `factor`.
"""
function commit!(matrix::ReducedMatrix, simplex, factor)
    @assert sign(simplex) == 1
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
        matrix.column_index[simplex] = length(matrix.indices) - 1
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

Base.haskey(matrix::ReducedMatrix, simplex) = haskey(matrix.column_index, abs(simplex))

function Base.getindex(matrix::ReducedMatrix{S}, element::AbstractChainElement{S}) where S
    return matrix[simplex(element)]
end
function Base.getindex(matrix::ReducedMatrix{S}, simplex::S) where S
    index = get(matrix.column_index, abs(simplex), 0)
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

# ======================================================================================== #
struct WorkingBoundary{E<:AbstractChainElement, O<:Base.Ordering}
    heap::Vector{E}
    ordering::O

    function WorkingBoundary{E}(ordering::O) where {E, O}
        new{E, O}(E[], ordering)
    end
end

Base.empty!(col::WorkingBoundary) = empty!(col.heap)
Base.isempty(col::WorkingBoundary) = isempty(col.heap)

function Base.sizehint!(col::WorkingBoundary, size)
    sizehint!(col.heap, size)
    return col
end

function Base.pop!(column::WorkingBoundary)
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

function Base.push!(column::WorkingBoundary{E}, simplex::AbstractSimplex) where E
    push!(column, E(simplex))
end
function Base.push!(column::WorkingBoundary{E}, element::E) where E
    heappush!(column.heap, element, column.ordering)
end

function nonheap_push!(column::WorkingBoundary{E}, simplex::AbstractSimplex) where E
    push!(column.heap, E(simplex))
end
function nonheap_push!(column::WorkingBoundary{E}, element::E) where E
    push!(column.heap, element)
end

repair!(column::WorkingBoundary) = heapify!(column.heap, column.ordering)
Base.first(column::WorkingBoundary) = first(column.heap)

function move!(column::WorkingBoundary{E}) where E
    dst = E[]
    while (pivot = pop!(column)) ≠ nothing
        push!(dst, pivot)
    end
    return dst
end

# ======================================================================================== #
# TODO: the following code is pretty messy. There must be a way to handle all this with a
# bit less type magic.
struct ReductionMatrix{
    Cohomology, T, Filtration, Simplex, SimplexElem, Facet, FacetElem, O<:Base.Ordering
}
    filtration::Filtration
    reduced::ReducedMatrix{Facet, SimplexElem, O}
    working_boundary::WorkingBoundary{FacetElem, O}
    columns_to_reduce::Vector{Simplex}
    columns_to_skip::Vector{Simplex}
end

function ReductionMatrix{C, T}(
    filtration::Filtration, columns_to_reduce, columns_to_skip
) where {C, T, Filtration}

    Simplex = eltype(columns_to_reduce)
    facet_dim = dim(Simplex) + (C ? 1 : -1)
    Facet = simplex_type(filtration, facet_dim)
    ordering = C ? Base.Order.Forward : Base.Order.Reverse
    O = typeof(ordering)
    SimplexElem = chain_element_type(Simplex, T)
    FacetElem = chain_element_type(Facet, T)

    reduced = ReducedMatrix{Facet, SimplexElem}(ordering)
    sizehint!(reduced, length(columns_to_reduce))
    working_boundary = WorkingBoundary{FacetElem}(ordering)

    return ReductionMatrix{C, T, Filtration, Simplex, SimplexElem, Facet, FacetElem, O}(
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

is_cohomology(::ReductionMatrix{C}) where C = C
field_type(::ReductionMatrix{<:Any, F}) where F = F
simplex_type(::ReductionMatrix{<:Any, <:Any, <:Any, S}) where S = S
simplex_element(::ReductionMatrix{<:Any, <:Any, <:Any, <:Any, E}) where E = E
facet_element(::ReductionMatrix{<:Any, <:Any, <:Any, <:Any, <:Any, <:Any, E}) where E = E
dim(::ReductionMatrix{true, <:Any, <:Any, S}) where S = dim(S)
dim(::ReductionMatrix{false, <:Any, <:Any, S}) where S = dim(S) - 1

function initialize_boundary!(matrix::ReductionMatrix, column_simplex)
    # TODO: can emergent pairs be safely be enabled for all kinds of complexes?
    # An argument that disables it could be added.
    empty!(matrix.working_boundary)
    for facet in co_boundary(matrix, column_simplex)
        # Checking this on every facet helps if more than one facet has the same diameter.
        if (is_cohomology(matrix) &&
            diam(facet) == diam(column_simplex) &&
            !haskey(matrix.reduced, facet)
            )
            empty!(matrix.working_boundary)
            return facet_element(matrix)(facet)
        end
        nonheap_push!(matrix.working_boundary, facet)
    end
    if isempty(matrix.working_boundary)
        return nothing
    else
        repair!(matrix.working_boundary)
        return pop!(matrix.working_boundary)
    end
end

function add!(matrix::ReductionMatrix, column, pivot)
    factor = -coefficient(pivot)
    record!(matrix.reduced, column, factor)
    for element in column
        for facet in co_boundary(matrix, simplex(element))
            abs(simplex(pivot)) == abs(facet) && continue
            push!(
                matrix.working_boundary,
                facet_element(matrix)(facet, coefficient(element) * factor)
            )
        end
    end
    return matrix
end

function reduce_column!(matrix::ReductionMatrix, column_simplex)
    pivot = initialize_boundary!(matrix, column_simplex)

    while !isnothing(pivot)
        column = matrix.reduced[pivot]
        isempty(column) && break

        add!(matrix, column, pivot)
        pivot = pop!(matrix.working_boundary)
    end
    if isnothing(pivot)
        discard!(matrix.reduced)
    else
        record!(matrix.reduced, simplex_element(matrix)(column_simplex))
        commit!(matrix.reduced, simplex(pivot), inv(coefficient(pivot)))
    end

    return pivot
end

function birth_death(::ReductionMatrix{true}, column, pivot)
    return diam(column), isnothing(pivot) ? Inf : diam(pivot)
end
function birth_death(::ReductionMatrix{false}, column, pivot)
    return isnothing(pivot) ? Inf : diam(pivot), diam(column)
end

function interval(
    ::Type{<:PersistenceInterval}, matrix::ReductionMatrix, column, pivot, cutoff
)
    birth, death = birth_death(matrix, column, pivot)
    return death - birth > cutoff ? PersistenceInterval(birth, death) : nothing
end
function interval(
    ::Type{R}, matrix::ReductionMatrix{true}, column, pivot, cutoff
) where R<:RepresentativeInterval
    birth, death = birth_death(matrix, column, pivot)

    if !isfinite(death)
        return R(
            PersistenceInterval(birth, death), column, nothing, eltype(matrix.reduced)[]
        )
    elseif death - birth > cutoff
        return R(
            PersistenceInterval(birth, death),
            column, simplex(pivot), collect(matrix.reduced[pivot])
        )
    else
        return nothing
    end
end
function interval(
    ::Type{R}, matrix::ReductionMatrix{false}, column, pivot, cutoff
) where R<:RepresentativeInterval
    birth, death = birth_death(matrix, column, pivot)

    if death - birth > cutoff
        push!(matrix.working_boundary, pivot)
        return R(
            PersistenceInterval(birth, death),
            simplex(pivot), column, move!(matrix.working_boundary)
        )
    else
        return nothing
    end
end

function interval_eltype(matrix::ReductionMatrix, ::Val{true})
    critical_birth_type = simplex_type(matrix.filtration, dim(matrix))
    critical_death_type = simplex_type(matrix.filtration, dim(matrix) + 1)
    representative_type = chain_element_type(critical_birth_type, field_type(matrix))
    return RepresentativeInterval{
        PersistenceInterval,
        critical_birth_type,
        Union{critical_death_type, Nothing},
        Vector{representative_type},
    }
end
function interval_eltype(::ReductionMatrix, ::Val{false})
    return PersistenceInterval
end

function compute_intervals!(
    matrix::ReductionMatrix, cutoff, progress, ::Val{reps}
) where {reps}
    if progress
        progbar = Progress(
            length(matrix.columns_to_reduce);
            desc="Computing $(dim(matrix))d intervals... ",
        )
    end
    intervals = interval_eltype(matrix, Val(reps))[]

    for column in matrix.columns_to_reduce
        pivot = reduce_column!(matrix, column)
        int = interval(eltype(intervals), matrix, column, pivot, cutoff)
        !isnothing(int) && push!(intervals, int)
        progress && next!(progbar; showvalues=((:n_intervals, length(intervals)),))
    end

    return sort!(PersistenceDiagram(
        dim(matrix),
        map(int -> postprocess_interval(matrix.filtration, int), intervals),
        threshold(matrix.filtration),
    ))
end

simplex_name(::Type{<:Simplex{2}}) = "triangles"
simplex_name(::Type{<:Simplex{3}}) = "tetrahedra"
simplex_name(::Type{<:AbstractSimplex{D}}) where D = "$D-simplices"

function next_matrix(matrix::ReductionMatrix, progress)
    C = simplex_type(matrix.filtration, dim(simplex_type(matrix)) + 1)
    new_to_reduce = C[]
    new_to_skip = C[]
    is_cohomology(matrix) && sizehint!(new_to_skip, length(matrix.reduced))

    if progress
        progbar = Progress(
            length(matrix.columns_to_reduce) + length(matrix.columns_to_skip),
            desc="Assembling...             "
        )
    end
    for simplex in Iterators.flatten((matrix.columns_to_reduce, matrix.columns_to_skip))
        for cofacet in coboundary(matrix.filtration, simplex, Val(false))
            # Clearing optimization only enabled for cohomology.
            if is_cohomology(matrix) && haskey(matrix.reduced, cofacet)
                push!(new_to_skip, abs(cofacet))
            else
                push!(new_to_reduce, abs(cofacet))
            end
        end
        progress && next!(progbar)
    end
    progress && printstyled(
        stderr, "Assembled $(length(new_to_reduce)) $(simplex_name(C)). Sorting... ",
        color=:green
    )
    sort!(new_to_reduce, rev=is_cohomology(matrix))
    progress && printstyled(stderr, "done.\n", color=:green)

    return ReductionMatrix{is_cohomology(matrix), field_type(matrix)}(
        matrix.filtration, new_to_reduce, new_to_skip
    )
end
