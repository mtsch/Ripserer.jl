"""
    WorkingChain{E<:AbstractChainElement, O<:Base.Ordering}

This structure hold the chain that is currently being reduced. Elements are added to the
chain with `push!`. The pivot element is the minimal element with respect to `O` and can be
extracted with `pop!`. `nonheap_push!` followed by `repair!` can also be used to add
elemetns to the chain.

`move!` is used to empty the chain and transfer its elements to a new array.
"""
struct WorkingChain{E<:AbstractChainElement, O<:Base.Ordering}
    heap::Vector{E}
    ordering::O

    function WorkingChain{E}(ordering::O) where {E, O}
        new{E, O}(E[], ordering)
    end
end

Base.empty!(col::WorkingChain) = empty!(col.heap)
Base.isempty(col::WorkingChain) = isempty(col.heap)

function Base.sizehint!(col::WorkingChain, size)
    sizehint!(col.heap, size)
    return col
end

function Base.pop!(column::WorkingChain)
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

function Base.push!(column::WorkingChain{E}, element::E) where E
    heappush!(column.heap, element, column.ordering)
end

function nonheap_push!(column::WorkingChain{E}, simplex::AbstractSimplex) where E
    push!(column.heap, E(simplex))
end
function nonheap_push!(column::WorkingChain{E}, element::E) where E
    push!(column.heap, element)
end

repair!(column::WorkingChain) = heapify!(column.heap, column.ordering)
Base.first(column::WorkingChain) = first(column.heap)

function move!(column::WorkingChain{E}) where E
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
    working_coboundary::WorkingChain{CofacetElem, O}
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
    working_coboundary = WorkingChain{CofacetElem}(ordering)

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
