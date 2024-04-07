"""
    initialize_column!(matrix, column_index)

Initialize the column indexed by `column_index` by adding its (co)boundary to
`matrix.chain`. This is where the emergent pairs optimization gets triggered for implicit
versions of the algorithm.

Return the pivot as a [`ChainElement`](@ref).
"""
function initialize_column!(matrix, column_index)
    empty!(matrix.chain)
    # Emergent pairs: we are looking for pairs of simplices (σ, τ) where σ is the youngest
    # facet of τ and τ is the oldest cofacet of σ. These pairs give birth to persistence
    # intervals with zero length and can be skipped.

    # This implementation of this optimization only works if (co)boundary simplices are
    # returned in the correct order and if the birth times of σ and τ are the same.
    emergent_check = emergent_pairs(matrix.filtration) && is_implicit(matrix)
    for cofacet in coboundary(matrix, column_index)
        if emergent_check && birth(cofacet) == birth(column_index)
            emergent_check = false
            if !haskey(matrix.reduced, cofacet)
                return ChainElement{typeof(cofacet),field_type(matrix)}(cofacet)
            end
        end
        push!(matrix.chain, cofacet)
    end
    if isempty(matrix.chain)
        return nothing
    else
        heapify!(matrix.chain, ordering(matrix))
        return heappop!(matrix.chain, ordering(matrix))
    end
end

"""
    add!(matrix, column, pivot)

Add already `column` multiplied by `-coefficient(pivot)` to `matrix.chain`.

"""
function add!(matrix, column, pivot)
    add!(Val(is_implicit(matrix)), matrix, column, pivot)
    return nothing
end
# Implicit version
function add!(::Val{true}, matrix, column, pivot)
    factor = -coefficient(pivot)
    for element in column
        for cofacet in coboundary(matrix, simplex(element))
            simplex(pivot) == cofacet && continue
            heappush!(
                matrix.chain, (cofacet, coefficient(element) * factor), ordering(matrix)
            )
        end
    end
    append!(matrix.buffer, (s, c * factor) for (s, c) in column)
    return nothing
end
# Explicit version
function add!(::Val{false}, matrix, column, pivot)
    factor = -coefficient(pivot)
    for element in column
        # The pivot is not stored in the column, so there is no need to check for it.
        heappush!(matrix.chain, element * factor, ordering(matrix))
    end
    return nothing
end

"""
    finalize!(matrix, column_index, pivot)

After reduction is done, finalize the current column being reduced by adding it to the
reduced matrix.
"""
function finalize!(matrix, column_index, pivot)
    finalize!(Val(is_implicit(matrix)), matrix, column_index, pivot)
    return nothing
end
# Implicit version
function finalize!(::Val{true}, matrix, column_index, pivot)
    push!(matrix.buffer, column_index)
    matrix.reduced[simplex(pivot)] = clean!(
        matrix.buffer, ordering(matrix), inv(coefficient(pivot))
    )
    return nothing
end
# Explicit version
function finalize!(::Val{false}, matrix, _, pivot)
    matrix.reduced[simplex(pivot)] = clean!(
        matrix.chain, ordering(matrix), inv(coefficient(pivot))
    )
    return nothing
end

"""
    mark_zero_column!(matrix, column_index)

This function is called on a (co)boundary matrix when a column is completely reduced.
"""
function mark_zero_column!(_, _)
    return nothing
end

"""
    reduce_column!(matrix, column_index)

Reduce `column_index` by repeatedly adding other columns to it. Once nothing more can be
added, `finalize!` the column.
"""
function reduce_column!(matrix, column_index)
    empty!(matrix.buffer)
    pivot = initialize_column!(matrix, column_index)
    while !isnothing(pivot)
        column = matrix.reduced[pivot]
        isempty(column) && break

        add!(matrix, column, pivot)
        pivot = heappop!(matrix.chain, ordering(matrix))
    end
    if !isnothing(pivot)
        finalize!(matrix, column_index, pivot)
    else
        mark_zero_column!(matrix, column_index)
    end

    return pivot
end

"""
    collect_cocycle!(matrix, column, pivot)

Collect the representative (co)cycle.
"""
function collect_cocycle!(matrix, column, pivot)
    if is_cohomology(matrix)
        if !is_implicit(matrix)
            error("representative cocycles for explicit cohomology not supported")
        elseif isnothing(pivot)
            push!(matrix.buffer, column)
            return copy(clean!(matrix.buffer, ordering(matrix)))
        else
            return Chain(matrix.reduced[pivot])
        end
    else
        if is_implicit(matrix)
            heappush!(matrix.chain, pivot, ordering(matrix))
            return heapmove!(matrix.chain, ordering(matrix))
        else
            rep = Chain(matrix.reduced[pivot])
            push!(rep, pivot)
            return clean!(rep, ordering(matrix))
        end
    end
end

"""
    interval(matrix, column, pivot, cutoff, reps)

Construct a persistence interval.
"""
function interval(matrix, column, pivot, cutoff, reps)
    if is_cohomology(matrix)
        birth_simplex = column
        death_simplex = isnothing(pivot) ? nothing : simplex(pivot)
    elseif isnothing(pivot)
        # In homology, birth simplex is nothing when column is fully reduced.
        return nothing
    else
        birth_simplex, death_simplex = simplex(pivot), column
    end
    birth_time = Float64(birth(birth_simplex))
    death_time = isnothing(death_simplex) ? Inf : Float64(birth(death_simplex))
    if death_time - birth_time > cutoff
        if reps
            rep = (; representative=collect_cocycle!(matrix, column, pivot))
        else
            rep = NamedTuple()
        end
        meta = (; birth_simplex=birth_simplex, death_simplex=death_simplex, rep...)
        return PersistenceInterval(birth_time, death_time, meta)
    else
        return nothing
    end
end

"""
    handle_apparent_pairs!(matrix, intervals, cutoff, verbose, reps)

Handle apparent pair precomputation, if defined for `matrix.filtration`. Only does anything
for implicit cohomology. Resulting intervals (if any) are `push!`ed to `intervals`.
"""
function handle_apparent_pairs!(matrix, intervals, cutoff, verbose, reps)
    coho = Val(is_cohomology(matrix))
    impl = Val(is_implicit(matrix))
    return handle_apparent_pairs!(coho, impl, matrix, intervals, cutoff, verbose, reps)
end
# Implicit cohomology version
function handle_apparent_pairs!(
    ::Val{true}, ::Val{true}, matrix, intervals, cutoff, verbose, reps
)
    columns, apparent = find_apparent_pairs(
        matrix.filtration, matrix.columns_to_reduce, verbose
    )
    bulk_insert!(matrix.reduced, apparent)
    for (σ, τ) in apparent
        τ_elem = ChainElement{typeof(τ),field_type(matrix)}(τ)
        int = interval(matrix, σ, τ_elem, cutoff, reps)
        !isnothing(int) && push!(intervals, int)
    end
    return columns
end
# Other versions don't support this
function handle_apparent_pairs!(::Val, ::Val, matrix, _, _, _, _)
    return matrix.columns_to_reduce
end

"""
    compute_intervals!(matrix, cutoff, verbose, reps)

Compute all intervals by fully reducing `matrix`.
"""
function compute_intervals!(matrix, cutoff, verbose, reps; sort_columns=true)
    intervals = PersistenceInterval[]

    #TODO: remove this feature
    columns = handle_apparent_pairs!(matrix, intervals, cutoff, verbose, reps)

    @prog_print(
        verbose,
        fmt_number(length(columns)),
        " ",
        simplex_name(eltype(columns)),
        " to reduce.",
    )
    if sort_columns
        sort_t = time_ns()
        sort_columns!(matrix)
        elapsed = round((time_ns() - sort_t) / 1e9; digits=3)
        @prog_println verbose " Sorted in " elapsed "s."
    else
        @prog_println verbose
    end

    if verbose
        progbar = Progress(length(columns); desc="Computing $(dim(matrix))d intervals... ")
    end
    for (i, column) in enumerate(columns)
        if attempt_early_stop!(matrix, i, columns)
            break
        end
        pivot = reduce_column!(matrix, column)
        int = interval(matrix, column, pivot, cutoff, reps)
        !isnothing(int) && push!(intervals, int)

        verbose && next!(progbar; showvalues=((:intervals, length(intervals)),))
    end
    verbose && finish!(progbar)
    if !is_cohomology(matrix)
        append_infinite_intervals!(intervals, matrix)
    end

    diagram = PersistenceDiagram(
        sort!(intervals; by=i -> (persistence(i), birth(i)));
        threshold=Float64(threshold(matrix.filtration)),
        dim=dim(matrix),
        field=field_type(matrix),
        filtration=matrix.filtration,
    )
    return postprocess_diagram(matrix.filtration, diagram)
end
