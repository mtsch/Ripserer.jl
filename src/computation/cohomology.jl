struct CoboundaryMatrix{T, F, S, R, C}
    filtration::F
    #reduced::ReducedMatrix{Cofacet, SimplexElem, O}
    reduced::R
    #chain::WorkingChain{CofacetElem, O}
    chain::C
    columns_to_reduce::Vector{S}
    columns_to_skip::Vector{S}
end

function CoboundaryMatrix(::Type{T}, filtration, columns_to_reduce, columns_to_skip) where T
    Simplex = eltype(columns_to_reduce)
    Cofacet = simplex_type(filtration, dim(Simplex) + 1)
    ordering = Base.Order.Forward
    SimplexElem = chain_element_type(Simplex, T)
    CofacetElem = chain_element_type(Cofacet, T)

    reduced = ReducedMatrix{Cofacet, SimplexElem}(ordering)
    sizehint!(reduced, length(columns_to_reduce))
    chain = WorkingChain{CofacetElem}(ordering)

    return CoboundaryMatrix{T, typeof(filtration), Simplex, typeof(reduced), typeof(chain)}(
        filtration, reduced, chain, columns_to_reduce, columns_to_skip
    )
end

function coboundary(matrix::CoboundaryMatrix, simplex::AbstractSimplex)
    return coboundary(matrix.filtration, simplex)
end

field_type(::CoboundaryMatrix{T}) where T = T
simplex_type(::CoboundaryMatrix{<:Any, <:Any, S}) where S = S
simplex_element(::CoboundaryMatrix{T, <:Any, S}) where {T, S} = chain_element_type(S, T)
dim(cm::CoboundaryMatrix) = dim(simplex_type(cm))
cofacet_type(cm::CoboundaryMatrix{<:Any, F}) where F = simplex_type(F, dim(cm) + 1)
cofacet_element(cm::CoboundaryMatrix{T}) where {T} = chain_element_type(cofacet_type(cm), T)

function initialize_coboundary!(matrix::CoboundaryMatrix, column_simplex)
    empty!(matrix.chain)
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
        nonheap_push!(matrix.chain, cofacet)
    end
    if isempty(matrix.chain)
        return nothing
    else
        repair!(matrix.chain)
        return pop!(matrix.chain)
    end
end

function add!(matrix::CoboundaryMatrix, column, pivot)
    factor = -coefficient(pivot)
    for element in column
        for cofacet in coboundary(matrix, simplex(element))
            simplex(pivot) == cofacet && continue
            push!(
                matrix.chain,
                cofacet_element(matrix)(cofacet, coefficient(element) * factor)
            )
        end
    end
    return matrix
end

function reduce_column!(matrix::CoboundaryMatrix, column_simplex)
    pivot = initialize_coboundary!(matrix, column_simplex)

    while !isnothing(pivot)
        column = matrix.reduced[pivot]
        isempty(column) && break

        add!(matrix, column, pivot)
        record!(matrix.reduced, column, -coefficient(pivot))
        pivot = pop!(matrix.chain)
    end
    if isnothing(pivot)
        discard!(matrix.reduced)
    else
        record!(matrix.reduced, column_simplex)
        commit!(matrix.reduced, simplex(pivot), inv(coefficient(pivot)))
    end

    return pivot
end

function add_interval!(
    intervals, matrix::CoboundaryMatrix, birth_sx, death_sx, cutoff, ::Val{reps}
) where reps
    birth_time = Float64(birth(birth_sx))
    death_time = !isnothing(death_sx) ? Float64(birth(death_sx)) : Inf
    if death_time - birth_time > cutoff
        if reps && isfinite(death_time)
            rep = (;representative=collect(matrix.reduced[death_sx]))
        elseif reps
            rep = (;representative=eltype(matrix.reduced)[])
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
    return nothing
end

function compute_intervals!(
    matrix::CoboundaryMatrix, cutoff, progress, ::Val{reps}
) where {reps}
    ###
    ### Set up output.
    ###
    intervals = interval_type(
        matrix.filtration, Val(dim(matrix)), Val(reps), field_type(matrix)
    )[]

    ###
    ### Apparent pair stuff.
    ###
    columns, apparent = find_apparent_pairs(
        matrix.filtration, matrix.columns_to_reduce, progress
    )
    bulk_add!(matrix.reduced, apparent)
    foreach(apparent) do (σ, τ)
        add_interval!(intervals, matrix, σ, cofacet_element(matrix)(τ), cutoff, Val(reps))
    end

    ###
    ### Interval computation.
    ###
    prog_print(
        progress, length(columns), " ", (simplex_name(eltype(columns))), " to reduce."
    )
    # One-dimensional columns are already sorted.
    if dim(matrix) > 1
        prog_print(progress, " Sorting...")
        sort_t = time_ns()
        sort!(columns, rev=true)
        elapsed = round((time_ns() - sort_t) / 1e9, digits=3)
        prog_println(progress, "done. (", elapsed, "seconds)")
    else
        prog_println(progress)
    end

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

function next_matrix(matrix::CoboundaryMatrix, progress)
    new_dim = dim(matrix) + 1
    C = cofacet_type(matrix)
    new_to_reduce = C[]
    new_to_skip = C[]
    sizehint!(new_to_skip, length(matrix.reduced))

    if progress
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
        progress && next!(progbar; showvalues=(
            ("cleared", length(new_to_skip)),
            ("to reduce", length(new_to_reduce)),
        ))
    end
    prog_print(progress, '\r')

    return CoboundaryMatrix(
        field_type(matrix), matrix.filtration, new_to_reduce, new_to_skip
    )
end

function cohomology(
    filtration, cutoff, progress, ::Type{F}, ::Val{dim_max}, ::Val{reps}
) where {F, dim_max, reps}
    result = PersistenceDiagram[]
    zeroth, to_reduce, to_skip = zeroth_intervals(
        filtration, cutoff, progress, F, Val(reps)
    )
    push!(result, zeroth)
    if dim_max > 0
        matrix = CoboundaryMatrix(F, filtration, to_reduce, to_skip)
        for dim in 1:dim_max
            push!(result, compute_intervals!(matrix, cutoff, progress, Val(reps)))
            if dim < dim_max
                matrix = next_matrix(matrix, progress)
            end
        end
    end
    return result
end
