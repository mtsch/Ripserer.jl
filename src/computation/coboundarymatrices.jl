function initialize_coboundary!(matrix, column)
    initialize_coboundary!(Val(is_cohomology(matrix)), matrix, column)
end
# Cohomology version
function initialize_coboundary!(::Val{true}, matrix, column)
    empty!(matrix.chain)
    # Emergent pairs: we are looking for pairs of simplices (σ, τ) where σ is the youngest
    # facet of τ and τ is the oldest cofacet of σ. These pairs give birth to persistence
    # intervals with zero length and can be skipped.

    # This implementation of this optimization only works if (co)boundary simplices are
    # returned in the correct order and if the birth times of σ and τ are the same.
    emergent_check = emergent_pairs(matrix.filtration) && is_implicit(matrix)
    for cofacet in coboundary(matrix, column)
        if emergent_check && birth(cofacet) == birth(column)
            emergent_check = false
            if !haskey(matrix.reduced, cofacet)
                return chain_element_type(matrix)(cofacet)
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
# Homology version
function initialize_coboundary!(::Val{false}, matrix, column)
    empty!(matrix.chain)
    for facet in coboundary(matrix, column)
        nonheap_push!(matrix.chain, facet)
    end
    if isempty(matrix.chain)
        return nothing
    else
        repair!(matrix.chain)
        return pop!(matrix.chain)
    end
end

function add!(matrix, column, pivot)
    add!(Val(is_implicit(matrix)), matrix, column, pivot)
end
# Implicit version
function add!(::Val{true}, matrix, column, pivot)
    factor = -coefficient(pivot)
    for element in column
        for cofacet in coboundary(matrix, simplex(element))
            simplex(pivot) == cofacet && continue
            push!(
                matrix.chain,
                chain_element_type(matrix)(cofacet, coefficient(element) * factor),
            )
        end
    end
    record!(matrix.reduced, column, factor)
end
# Explicit version
function add!(::Val{false}, matrix, column, pivot)
    factor = -coefficient(pivot)
    for element in column
        # The pivot is not stored in the column, so there is no need to check for it.
        push!(matrix.chain, element * factor)
    end
end

function finalize!(matrix, column, pivot)
    finalize!(Val(is_implicit(matrix)), matrix, column, pivot)
end
# Implicit version
function finalize!(::Val{true}, matrix, column, pivot)
    record!(matrix.reduced, column)
    commit!(matrix.reduced, simplex(pivot), inv(coefficient(pivot)))
end
# Explicit version
function finalize!(::Val{false}, matrix, column, pivot)
    record!(matrix.reduced, matrix.chain)
    commit!(matrix.reduced, simplex(pivot), inv(coefficient(pivot)))
end

function reduce_column!(matrix, column_to_reduce)
    pivot = initialize_coboundary!(matrix, column_to_reduce)

    while !isnothing(pivot)
        column = matrix.reduced[pivot]
        isempty(column) && break

        add!(matrix, column, pivot)
        pivot = pop!(matrix.chain)
    end
    if isnothing(pivot)
        discard!(matrix.reduced)
    else
        finalize!(matrix, column_to_reduce, pivot)
    end

    return pivot
end

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
    if !reps
        representative = nothing
    elseif is_cohomology(matrix)
        if isnothing(pivot)
            representative = simplex_type(matrix.filtration, dim(matrix))[]
        elseif is_implicit(matrix)
            representative = collect(matrix.reduced[pivot])
        else
            elem_t = chain_element_type(simplex_type(matrix.filtration, dim(matrix)), field_type(matrix))
            tmp_chain = WorkingChain{elem_t}(Base.Order.Forward)
            for elem in matrix.reduced[pivot]
                for facet in boundary(matrix.filtration, simplex(elem))
                    push!(tmp_chain, elem_t(facet, coefficient(elem)))
                end
            end
            representative = move!(tmp_chain)
        end
    else
        if is_implicit(matrix)
            representative = pushfirst!(move!(matrix.chain), birth_simplex)
        else
            representative = pushfirst!(collect(matrix.reduced[pivot]), birth_simplex)
        end
    end
    return interval(
        Val(is_implicit(matrix)), birth_simplex, death_simplex, cutoff, representative
    )
end
# Implicit version
function interval(::Val{true}, birth_simplex, death_simplex, cutoff, representative)
    birth_time = Float64(birth(birth_simplex))
    death_time = isnothing(death_simplex) ? Inf : Float64(birth(death_simplex))
    if death_time - birth_time > cutoff
        if !isnothing(representative)
            rep = (;representative=representative)
        else
            rep = NamedTuple()
        end
        meta = (;birth_simplex=birth_simplex, death_simplex=death_simplex, rep...)
        return PersistenceInterval(birth_time, death_time, meta)
    else
        return nothing
    end
end
# Explicit version
function interval(::Val{false}, birth_simplex, death_simplex, cutoff, representative)
    birth_time = Float64(birth(birth_simplex))
    death_time = isnothing(death_simplex) ? Inf : Float64(birth(death_simplex))
    if death_time - birth_time > cutoff
        if !isnothing(representative)
            rep = (;representative=representative)
        else
            rep = NamedTuple()
        end
        meta = (;birth_simplex=birth_simplex, death_simplex=death_simplex, rep...)
        return PersistenceInterval(birth_time, death_time, meta)
    else
        return nothing
    end
end

function handle_apparent_pairs!(matrix, intervals, cutoff, progress, reps)
    coho = Val(is_cohomology(matrix))
    impl = Val(is_implicit(matrix))
    return handle_apparent_pairs!(coho, impl, matrix, intervals, cutoff, progress, reps)
end
# Implicit cohomology version
function handle_apparent_pairs!(
    ::Val{true}, ::Val{true}, matrix, intervals, cutoff, progress, reps
)
    columns, apparent = find_apparent_pairs(
        matrix.filtration, matrix.columns_to_reduce, progress
    )
    bulk_add!(matrix.reduced, apparent)
    for (σ, τ) in apparent
        int = interval(matrix, σ, chain_element_type(matrix)(τ), cutoff, reps)
        !isnothing(int) && push!(intervals, int)
    end
    return columns
end
# Other versions
function handle_apparent_pairs!(::Val, ::Val, matrix, _, _, _, _)
    return matrix.columns_to_reduce
end

function compute_intervals!(matrix, cutoff, progress, reps)
    ###
    ### Set up output.
    ###
    intervals = interval_type(
        matrix.filtration, Val(dim(matrix)), Val(reps), field_type(matrix)
    )[]

    ###
    ### Apparent pair stuff.
    ###
    columns = handle_apparent_pairs!(matrix, intervals, cutoff, progress, reps)

    ###
    ### Interval computation.
    ###
    prog_print(
        progress, length(columns), " ", (simplex_name(eltype(columns))), " to reduce."
    )
    # One-dimensional columns in cohomology are already sorted.
    if !is_cohomology(matrix) || dim(matrix) > 1
        prog_print(progress, " Sorting... ")
        sort_t = time_ns()
        sort!(columns, rev=is_cohomology(matrix))
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
        int = interval(matrix, column, pivot, cutoff, reps)
        !isnothing(int) && push!(intervals, int)

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

###
###
###
struct CoboundaryMatrix{I, T, F, S, R, C}
    filtration::F
    reduced::R
    chain::C
    columns_to_reduce::Vector{S}
    columns_to_skip::Vector{S}
end

function CoboundaryMatrix{I}(
    ::Type{T}, filtration, columns_to_reduce, columns_to_skip
) where {I, T}
    Simplex = eltype(columns_to_reduce)
    Cofacet = simplex_type(filtration, dim(Simplex) + 1)
    ordering = Base.Order.Forward
    SimplexElem = chain_element_type(Simplex, T)
    CofacetElem = chain_element_type(Cofacet, T)

    if I
        reduced = ReducedMatrix{Cofacet, SimplexElem}(ordering)
    else
        reduced = ReducedMatrix{Cofacet, CofacetElem}(ordering)
    end
    sizehint!(reduced, length(columns_to_reduce))
    chain = WorkingChain{CofacetElem}(ordering)

    return CoboundaryMatrix{
        I, T, typeof(filtration), Simplex, typeof(reduced), typeof(chain)
    }(
        filtration, reduced, chain, columns_to_reduce, columns_to_skip
    )
end


field_type(::CoboundaryMatrix{<:Any, T}) where T = T
dim(cm::CoboundaryMatrix{<:Any, <:Any, <:Any, S}) where S = dim(S)
function chain_element_type(bm::CoboundaryMatrix{<:Any, T, F}) where {T, F}
    chain_element_type(simplex_type(F, dim(bm) + 1), T)
end

is_implicit(::CoboundaryMatrix{I}) where I = I
is_cohomology(::CoboundaryMatrix) = true

coboundary(matrix::CoboundaryMatrix, simplex::AbstractSimplex) = coboundary(matrix.filtration, simplex)

function next_matrix(matrix::CoboundaryMatrix{I}, progress) where I
    C = simplex_type(matrix.filtration, dim(matrix) + 1)
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

    return CoboundaryMatrix{I}(
        field_type(matrix), matrix.filtration, new_to_reduce, new_to_skip
    )
end

###
###
###
struct BoundaryMatrix{I, T, F, S, R, C}
    filtration::F
    reduced::R
    chain::C
    columns_to_reduce::Vector{S}
end

function BoundaryMatrix{I}(::Type{T}, filtration, columns_to_reduce) where {I, T}
    Simplex = typeof(first(columns_to_reduce))
    Facet = simplex_type(filtration, dim(Simplex) - 1)
    ordering = Base.Order.Reverse
    SimplexElem = chain_element_type(Simplex, T)
    FacetElem = chain_element_type(Facet, T)

    columns = Simplex[]
    foreach(columns_to_reduce) do c
        push!(columns, c)
    end

    if !I
        reduced = ReducedMatrix{Facet, FacetElem}(ordering)
    else
        reduced = ReducedMatrix{Facet, SimplexElem}(ordering)
    end
    sizehint!(reduced, length(columns))
    chain = WorkingChain{FacetElem}(ordering)

    return BoundaryMatrix{
        I, T, typeof(filtration), Simplex, typeof(reduced), typeof(chain)
    }(
        filtration, reduced, chain, columns
    )
end

field_type(::BoundaryMatrix{<:Any, T}) where T = T
dim(bm::BoundaryMatrix{<:Any, <:Any, <:Any, S}) where S = dim(S) - 1
function chain_element_type(bm::BoundaryMatrix{<:Any, T, F}) where {T, F}
    chain_element_type(simplex_type(F, dim(bm)), T)
end

is_implicit(::BoundaryMatrix{I}) where I = I
is_cohomology(::BoundaryMatrix) = false

function coboundary(matrix::BoundaryMatrix, simplex::AbstractSimplex)
    return boundary(matrix.filtration, simplex)
end
