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
function record!(matrix::ReducedMatrix{<:Any, E}, simplex::AbstractSimplex) where E
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
    bulk_add!(matrix::ReducedMatrix{<:Any, E}, pairs) where E

Insert apparent pairs into reduced matrix. Columns are inserted so `matrix[τ] == [σ]`, where
`pairs` is a collection of tuples `(σ, τ)`.
"""
function bulk_add!(matrix::ReducedMatrix{<:Any, E}, pairs) where E
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
    matrix
end
function bulk_add!(matrix::ReducedMatrix, ::Tuple{})
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

Base.haskey(matrix::ReducedMatrix, simplex) = haskey(matrix.column_index, simplex)

function Base.getindex(matrix::ReducedMatrix{S}, element::AbstractChainElement{S}) where S
    return matrix[simplex(element)]
end
function Base.getindex(matrix::ReducedMatrix{S}, simplex::S) where S
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

# ======================================================================================== #
struct WorkingCoboundary{E<:AbstractChainElement, O<:Base.Ordering}
    heap::Vector{E}
    ordering::O

    function WorkingCoboundary{E}(ordering::O) where {E, O}
        new{E, O}(E[], ordering)
    end
end

Base.empty!(col::WorkingCoboundary) = empty!(col.heap)
Base.isempty(col::WorkingCoboundary) = isempty(col.heap)

function Base.sizehint!(col::WorkingCoboundary, size)
    sizehint!(col.heap, size)
    return col
end

function Base.pop!(column::WorkingCoboundary)
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

function Base.push!(column::WorkingCoboundary{E}, element::E) where E
    heappush!(column.heap, element, column.ordering)
end

function nonheap_push!(column::WorkingCoboundary{E}, simplex::AbstractSimplex) where E
    push!(column.heap, E(simplex))
end
function nonheap_push!(column::WorkingCoboundary{E}, element::E) where E
    push!(column.heap, element)
end

repair!(column::WorkingCoboundary) = heapify!(column.heap, column.ordering)
Base.first(column::WorkingCoboundary) = first(column.heap)

function move!(column::WorkingCoboundary{E}) where E
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
    Cohomology, T, Filtration, Simplex, SimplexElem, Cofacet, CofacetElem, O<:Base.Ordering
}
    filtration::Filtration
    reduced::ReducedMatrix{Cofacet, SimplexElem, O}
    working_coboundary::WorkingCoboundary{CofacetElem, O}
    columns_to_reduce::Vector{Simplex}
    columns_to_skip::Vector{Simplex}
end

function ReductionMatrix{C, T}(
    filtration::Filtration, columns_to_reduce, columns_to_skip
) where {C, T, Filtration}

    Simplex = eltype(columns_to_reduce)
    cofacet_dim = dim(Simplex) + (C ? 1 : -1)
    Cofacet = simplex_type(filtration, cofacet_dim)
    ordering = C ? Base.Order.Forward : Base.Order.Reverse
    O = typeof(ordering)
    SimplexElem = chain_element_type(Simplex, T)
    CofacetElem = chain_element_type(Cofacet, T)

    reduced = ReducedMatrix{Cofacet, SimplexElem}(ordering)
    sizehint!(reduced, length(columns_to_reduce))
    working_coboundary = WorkingCoboundary{CofacetElem}(ordering)

    return ReductionMatrix{C, T, Filtration, Simplex, SimplexElem, Cofacet, CofacetElem, O}(
        filtration,
        reduced,
        working_coboundary,
        columns_to_reduce,
        columns_to_skip
    )
end

"""
    coboundary(matrix, simplex)

Iterate over the (co)boundary of the `simplex`. Chooses between the boundary and the
coboundary based on the matrix type.
"""
function coboundary(matrix::ReductionMatrix{true}, simplex::AbstractSimplex)
    return coboundary(matrix.filtration, simplex)
end
function coboundary(matrix::ReductionMatrix{false}, simplex::AbstractSimplex)
    return boundary(matrix.filtration, simplex)
end

is_cohomology(::ReductionMatrix{C}) where C = C
field_type(::ReductionMatrix{<:Any, F}) where F = F
simplex_type(::ReductionMatrix{<:Any, <:Any, <:Any, S}) where S = S
simplex_element(::ReductionMatrix{<:Any, <:Any, <:Any, <:Any, E}) where E = E
cofacet_element(::ReductionMatrix{<:Any, <:Any, <:Any, <:Any, <:Any, <:Any, E}) where E = E
dim(::ReductionMatrix{true, <:Any, <:Any, S}) where S = dim(S)
dim(::ReductionMatrix{false, <:Any, <:Any, S}) where S = dim(S) - 1

function initialize_coboundary!(matrix::ReductionMatrix, column_simplex)
    empty!(matrix.working_coboundary)
    # Emergent pairs: we are looking for pairs of simplices (σ, τ) where σ is the youngest
    # facet of τ and τ is the oldest cofacet of σ. These pairs give birth to persistence
    # intervals with zero length and can be skipped.

    # This implementation of this optimization only works if (co)boundary simplices are
    # returned in the correct order and if the birth times of σ and τ are the same.
    emergent_check = emergent_pairs(matrix.filtration)
    for cofacet in coboundary(matrix, column_simplex)
        if emergent_check && birth(cofacet) == birth(column_simplex)
            emergent_check = false
            if !haskey(matrix.reduced, cofacet)
                return cofacet_element(matrix)(cofacet)
            end
        end
        nonheap_push!(matrix.working_coboundary, cofacet)
    end
    if isempty(matrix.working_coboundary)
        return nothing
    else
        repair!(matrix.working_coboundary)
        return pop!(matrix.working_coboundary)
    end
end

function add!(matrix::ReductionMatrix, column, pivot)
    factor = -coefficient(pivot)
    for element in column
        for cofacet in coboundary(matrix, simplex(element))
            simplex(pivot) == cofacet && continue
            push!(
                matrix.working_coboundary,
                cofacet_element(matrix)(cofacet, coefficient(element) * factor)
            )
        end
    end
    return matrix
end

function reduce_column!(matrix::ReductionMatrix, column_simplex)
    pivot = initialize_coboundary!(matrix, column_simplex)

    while !isnothing(pivot)
        column = matrix.reduced[pivot]
        isempty(column) && break

        add!(matrix, column, pivot)
        record!(matrix.reduced, column, -coefficient(pivot))
        pivot = pop!(matrix.working_coboundary)
    end
    if isnothing(pivot)
        discard!(matrix.reduced)
    else
        record!(matrix.reduced, column_simplex)
        commit!(matrix.reduced, simplex(pivot), inv(coefficient(pivot)))
    end

    return pivot
end

function birth_death(::ReductionMatrix{true}, column, pivot)
    return (
        Float64(birth(column)),
        column,
        isnothing(pivot) ? Inf : Float64(birth(pivot)),
        isnothing(pivot) ? nothing : simplex(pivot),
    )
end
function birth_death(::ReductionMatrix{false}, column, pivot)
    return (
        isnothing(pivot) ? Inf : Float64(birth(pivot)),
        isnothing(pivot) ? nothing : simplex(pivot),
        Float64(birth(column)),
        column,
    )
end

function add_interval!(
    intervals, matrix::ReductionMatrix, column, pivot, cutoff, ::Val{reps}
) where reps
    birth_time, birth_sx, death_time, death_sx = birth_death(matrix, column, pivot)
    if death_time - birth_time > cutoff
        if reps && is_cohomology(matrix)
            if isfinite(death_time)
                rep = (;representative=collect(matrix.reduced[pivot]))
            else
                rep = (;representative=eltype(matrix.reduced)[])
            end
        elseif reps
            push!(matrix.working_coboundary, pivot)
            rep = (;representative=move!(matrix.working_coboundary))
        else
            rep = NamedTuple()
        end
        meta = (;
            birth_simplex=birth_sx,
            death_simplex=death_sx,
            rep...,
        )
        push!(intervals, PersistenceInterval(
            birth_time, death_time, meta
        ))
    end
end

function compute_intervals!(
    matrix::ReductionMatrix, cutoff, progress, ::Val{reps}
) where {reps}
    # Set up result.
    intervals = interval_type(
        matrix.filtration, Val(dim(matrix)), Val(reps), field_type(matrix)
    )[]

    # Apparent pair stuff.
    if is_cohomology(matrix)
        columns, apparent = find_apparent_pairs(
            matrix.filtration, matrix.columns_to_reduce, progress
        )
        bulk_add!(matrix.reduced, apparent)
        foreach(apparent) do (σ, τ)
            add_interval!(intervals, matrix, σ, cofacet_element(matrix)(τ), cutoff, Val(reps))
        end
    else
        columns = matrix.columns_to_reduce
    end

    # Interval computation.
    progress && printstyled(
        stderr, "$(length(columns)) $(simplex_name(eltype(columns))) to reduce. Sorting... ",
        color=:green
    )
    # One-dimensional columns are already sorted.
    sort_t = time_ns()
    if !is_cohomology(matrix) || dim(matrix) > 1
        sort!(columns, rev=is_cohomology(matrix))
    end
    elapsed = round((time_ns() - sort_t) / 1e9, digits=3)
    progress && printstyled(stderr, "done. ($elapsed seconds)\n", color=:green)

    if progress
        progbar = Progress(
            length(columns);
            desc="Computing $(dim(matrix))d intervals... ",
        )
    end
    for column in columns
        pivot = reduce_column!(matrix, column)
        add_interval!(intervals, matrix, column, pivot, cutoff, Val(reps))
        progress && next!(progbar; showvalues=((:intervals, length(intervals)),))
    end

    return postprocess_diagram(
        matrix.filtration, PersistenceDiagram(
            intervals;
            threshold=Float64(threshold(matrix.filtration)),
            dim=dim(matrix),
            field_type=field_type(matrix),
            filtration=matrix.filtration,
        )
    )
end

simplex_name(::Type{<:Simplex{2}}) = "triangles"
simplex_name(::Type{<:Simplex{3}}) = "tetrahedra"
simplex_name(::Type{<:AbstractSimplex{D}}) where D = "$D-simplices"

function next_matrix(matrix::ReductionMatrix, progress)
    new_dim = dim(simplex_type(matrix)) + 1
    C = simplex_type(matrix.filtration, dim(simplex_type(matrix)) + 1)
    new_to_reduce = C[]
    new_to_skip = C[]
    is_cohomology(matrix) && sizehint!(new_to_skip, length(matrix.reduced))

    if progress
        progbar = ProgressUnknown("Assembling columns:")
    end
    for simplex in columns_to_reduce(
        matrix.filtration,
        Iterators.flatten((matrix.columns_to_reduce, matrix.columns_to_skip)),
    )
        if is_cohomology(matrix) && haskey(matrix.reduced, simplex)
            push!(new_to_skip, abs(simplex))
        else
            push!(new_to_reduce, abs(simplex))
        end
        progress && next!(progbar; showvalues=(
            ("cleared", length(new_to_skip)),
            ("to reduce", length(new_to_reduce)),
        ))
    end
    progress && print(stderr, "\r")

    return ReductionMatrix{is_cohomology(matrix), field_type(matrix)}(
        matrix.filtration, new_to_reduce, new_to_skip
    )
end
