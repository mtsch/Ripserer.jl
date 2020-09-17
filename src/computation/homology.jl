struct BoundaryMatrix{T, F, S, R, C}
    filtration::F
    reduced::R
    chain::C
    columns_to_reduce::Vector{S}
end

function BoundaryMatrix(::Type{T}, filtration, columns_to_reduce, precomputed=()) where T
    #TODO: move precomputed to reduced
    Simplex = typeof(first(columns_to_reduce))
    Facet = simplex_type(filtration, dim(Simplex) - 1)
    ordering = Base.Order.Reverse
    SimplexElem = chain_element_type(Simplex, T)
    FacetElem = chain_element_type(Facet, T)

    columns = Simplex[]
    foreach(columns_to_reduce) do c
        push!(columns, c)
    end

    reduced = ReducedMatrix{Facet, FacetElem}(ordering)
    sizehint!(reduced, length(columns) + length(precomputed))
    chain = WorkingChain{FacetElem}(ordering)

    return BoundaryMatrix{T, typeof(filtration), Simplex, typeof(reduced), typeof(chain)}(
        filtration, reduced, chain, columns
    )
end

function boundary(matrix::BoundaryMatrix, simplex::AbstractSimplex)
    return boundary(matrix.filtration, simplex)
end

field_type(::BoundaryMatrix{T}) where T = T
simplex_type(::BoundaryMatrix{<:Any, <:Any, S}) where S = S
simplex_element(::BoundaryMatrix{T, <:Any, S}) where {T, S} = chain_element_type(S, T)
dim(bm::BoundaryMatrix) = dim(simplex_type(bm)) - 1
facet_type(bm::BoundaryMatrix{<:Any, F}) where F = simplex_type(F, dim(bm) - 1)
facet_element(bm::BoundaryMatrix{T}) where {T} = chain_element_type(facet_type(bm), T)

function initialize_boundary!(matrix::BoundaryMatrix, column_simplex)
    empty!(matrix.chain)
    for facet in boundary(matrix, column_simplex)
        nonheap_push!(matrix.chain, facet)
    end
    if isempty(matrix.chain)
        return nothing
    else
        sort!(matrix.chain, alg=InsertionSort)
        return pop!(matrix.chain)
    end
end

function add!(matrix::BoundaryMatrix, column, pivot)
    factor = -coefficient(pivot)
    for element in column
        simplex(pivot) == simplex(element) && continue
        push!(matrix.chain, element * factor)
    end
    return matrix
end

function reduce_column!(matrix::BoundaryMatrix, column_simplex)
    pivot = initialize_boundary!(matrix, column_simplex)

    while !isnothing(pivot)
        column = matrix.reduced[pivot]
        isempty(column) && break

        add!(matrix, column, pivot)
        #record!(matrix.reduced, column, -coefficient(pivot))
        pivot = pop!(matrix.chain)
    end

    if isnothing(pivot)
        discard!(matrix.reduced)
    else
        push!(matrix.chain, pivot)
        record!(matrix.reduced, matrix.chain, -coefficient(pivot))
        commit!(matrix.reduced, simplex(pivot), inv(coefficient(pivot)))
    end

    return pivot
end

function add_interval!(intervals, matrix::BoundaryMatrix, birth_sx, death_sx, cutoff)
    if !isnothing(birth_sx)
        birth_time = Float64(birth(birth_sx))
        death_time = Float64(birth(death_sx))
        if death_time - birth_time > cutoff
            meta = (;
                    birth_simplex=birth_sx,
                    death_simplex=death_sx,
                    representative=matrix.reduced[birth_sx],
                    )
            push!(intervals, PersistenceInterval(
                birth_time, death_time, meta
            ))
        end
    end
    return nothing
end

function compute_intervals!(matrix::BoundaryMatrix, cutoff, progress)
    ###
    ### Set up result.
    ###
    intervals = interval_type(
        matrix.filtration, Val(dim(matrix)), Val(true), field_type(matrix)
    )[]

    columns = matrix.columns_to_reduce

    ###
    ### Interval computation.
    ###
    prog_print(
        progress, length(columns), " ", (simplex_name(eltype(columns))), " to reduce."
    )
    sort_t = time_ns()
    prog_print(progress, " Sorting...")
    sort!(columns)
    elapsed = round((time_ns() - sort_t) / 1e9, digits=3)
    prog_println(progress, "done. (", elapsed, "seconds)")

    if progress
        progbar = Progress(
            length(columns);
            desc="Computing $(dim(matrix))d intervals... ",
        )
    end
    for column in columns
        pivot = reduce_column!(matrix, column)
        add_interval!(intervals, matrix, pivot, column, cutoff)
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

function homology(
    filtration, cutoff, progress, ::Type{F}, ::Val{dim_max}, ::Val{reps}
) where {F, dim_max, reps}
    result = PersistenceDiagram[]
    zeroth, to_reduce, to_skip = zeroth_intervals(
        filtration, cutoff, progress, F, Val(reps)
    )
    push!(result, zeroth)
    if dim_max > 0
        simplices = columns_to_reduce(filtration, Iterators.flatten((to_reduce, to_skip)))
        for dim in 1:dim_max
            if isempty(simplices)
                return result
            end
            matrix = BoundaryMatrix(F, filtration, simplices)
            push!(result, compute_intervals!(matrix, cutoff, progress))
            if dim < dim_max
                simplices = columns_to_reduce(filtration, simplices)
            end
        end
    end
    return result
end
