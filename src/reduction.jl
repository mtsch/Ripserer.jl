"""
    ReductionState

This structure represents the reduction matrix in the current dimension. A new one is
created for every dimension.

# Fields:

* `filtration`: the filtration we are analyzing.
* `reduction_matrix`: the reduction matrix. Each column of the matrix records the operations
  that were performed when reducing the column.
* `working_column`: the current working column, the column we are currently reducing.
* `reduction_entries`: this is where we record which simplices we added to the working
  column.
"""
struct ReductionState{
    Field,
    S<:AbstractSimplex, SE<:AbstractChainElement{S, Field},
    C<:AbstractSimplex, CE<:AbstractChainElement{C, Field},
    F,
}
    filtration        ::F
    reduction_matrix  ::ReductionMatrix{C, SE}
    working_column    ::Column{CE}
    reduction_entries ::Column{SE}

    function ReductionState{Field, S}(filtration::F) where {Field, S<:AbstractSimplex, F}
        SE = chain_element_type(S, Field)
        C = coface_type(S)
        CE = chain_element_type(C, Field)

        new{Field, S, SE, C, CE, F}(
            filtration, ReductionMatrix{C, SE}(), Column{CE}(), Column{SE}(),
        )
    end
end

simplex_type(rs::ReductionState{<:Any, S}) where S =
    S
coface_element(rs::ReductionState{<:Any, <:Any, <:Any, <:Any, CE}) where CE =
    CE

"""
    add!(rs::ReductionState, current_pivot)

Add column with column in `rs.reduction_matrix` indexed by `current_pivot` and multiplied by
the correct factor to `rs.working_column`. Also record the addition in
`rs.reduction_entries`.
"""
function add!(rs::ReductionState, current_pivot)
    λ = -coefficient(current_pivot)
    for element in rs.reduction_matrix[current_pivot]
        push!(rs.reduction_entries, λ * element)
        for coface in coboundary(rs.filtration, simplex(element))
            push!(rs.working_column, coface_element(rs)(coface, λ * coefficient(element)))
        end
    end
    pivot(rs.working_column)
end

"""
    initialize!(rs::ReductionState, column_simplex)

Initialize the columns in `rs`. Empty both `rs.working_column` and `rs.reduction_entries`,
then push the coboundary of `column_simplex` to `rs.working_column`.
"""
function initialize!(rs::ReductionState, column_simplex::AbstractSimplex)
    empty!(rs.working_column)
    empty!(rs.reduction_entries)

    for coface in coboundary(rs.filtration, column_simplex)
        if diam(coface) == diam(column_simplex) && !has_column(rs.reduction_matrix, coface)
            empty!(rs.working_column)
            return coface_element(rs)(coface)
        end
        push!(rs.working_column, coface)
    end
    pivot(rs.working_column)
end

"""
    reduce_working_column!(rs::ReductionState, column_simplex, Val(reps))

Reduce the working column by adding other columns to it until it has the lowest pivot or is
reduced. Record it in the reduction matrix and return the persistence interval. If `reps` is
`true`, add representative cocycles to the interval.
"""
function reduce_working_column!(
    rs::ReductionState{F, S}, column_simplex, cutoff, ::Val{reps}
) where {F, S, reps}
    current_pivot = initialize!(rs, column_simplex)

    while !isnothing(current_pivot) && has_column(rs.reduction_matrix, current_pivot)
        current_pivot = add!(rs, current_pivot)
    end
    if isnothing(current_pivot)
        death = ∞
    else
        insert_column!(rs.reduction_matrix, current_pivot)
        push!(rs.reduction_entries, column_simplex)
        move_mul!(rs.reduction_matrix, rs.reduction_entries, inv(coefficient(current_pivot)))
        death = diam(simplex(current_pivot))
    end
    birth = diam(column_simplex)

    if reps && !isfinite(death)
        PersistenceInterval(birth, death, chain_element_type(S, F)[])
    elseif reps && death - birth > cutoff
        representative = collect(rs.reduction_matrix[current_pivot])
        PersistenceInterval(birth, death, representative)
    elseif !isfinite(death) || death - birth > cutoff
        PersistenceInterval(birth, death)
    else
        nothing
    end
end

"""
    compute_pairs!(rs::ReductionState, columns)

Compute persistence intervals by reducing `columns`, a collection of simplices. Return
`PersistenceDiagram`. Only keep intervals with desired birth/death `cutoff` and return
representative cocycles if `reps` is `true`.
"""
function compute_intervals!(
    rs::ReductionState{F, S}, columns, cutoff, ::Val{reps}
) where {S, F, reps}
    T = dist_type(rs.filtration)
    if reps
        intervals = PersistenceInterval{T, Vector{chain_element_type(S, F)}}[]
    else
        intervals = PersistenceInterval{T, Nothing}[]
    end
    for column in columns
        interval = reduce_working_column!(rs, column, cutoff, Val(reps))
        if !isnothing(interval)
            push!(intervals, interval)
        end
    end
    sort!(PersistenceDiagram(dim(eltype(columns)), intervals))
end

"""
    assemble_columns!(rs::ReductionState, columns, simplices)

Assemble columns that need to be reduced in the next dimension. Apply clearing optimization.
"""
function assemble_columns!(rs::ReductionState{<:Any, S}, simplices) where S
    C = coface_type(S)
    columns = C[]
    new_simplices = C[]

    for simplex in simplices
        for coface in coboundary(rs.filtration, simplex, Val(false))
            push!(new_simplices, abs(coface))
            if !has_column(rs.reduction_matrix, coface)
                push!(columns, abs(coface))
            end
        end
    end
    sort!(columns, rev=true)
    columns, new_simplices
end

"""
    zeroth_intervals(filtration, cutoff, field_type, ::Val{reps})

Compute 0-dimensional persistent homology using Kruskal's Algorithm.

Only keep intervals with desired birth/death `cutoff`. Compute homology with coefficients in
`field_type`. If `reps` is `true`, compute representative cocycles.
"""
function zeroth_intervals(filtration, cutoff, field_type, ::Val{reps}) where reps
    T = dist_type(filtration)
    V = vertex_type(filtration)
    CE = chain_element_type(V, field_type)
    dset = DisjointSetsWithBirth([birth(filtration, v) for v in 1:n_vertices(filtration)])
    if reps
        intervals = PersistenceInterval{T, Vector{CE}}[]
    else
        intervals = PersistenceInterval{T, Nothing}[]
    end
    simplices = edge_type(filtration)[]
    columns = edge_type(filtration)[]

    for sx in edges(filtration)
        u, v = vertices(sx)
        i = find_root!(dset, u)
        j = find_root!(dset, v)
        push!(simplices, sx)
        if i ≠ j

            # According to the elder rule, the vertex with the lower birth will fall
            # into a later interval.
            dead = birth(dset, i) > birth(dset, j) ? i : j
            if diam(sx) - birth(dset, dead) > cutoff
                if reps
                    representative = map(find_leaves!(dset, dead)) do w
                        CE(V((w,), birth(filtration, w)))
                    end
                else
                    representative = nothing
                end
                interval = PersistenceInterval(birth(dset, dead), diam(sx), representative)
                push!(intervals, interval)
            end
            union!(dset, i, j)
        else
            push!(columns, sx)
        end
    end
    for v in 1:n_vertices(filtration)
        if find_root!(dset, v) == v
            if reps
                representative = map(find_leaves!(dset, v)) do w
                    CE(V((w,), birth(filtration, w)))
                end
            else
                representative = nothing
            end
            push!(intervals, PersistenceInterval(birth(dset, v), ∞, representative))
        end
    end
    reverse!(columns)
    sort!(PersistenceDiagram(0, intervals)), columns, simplices
end

"""
    nth_intervals(filtration, columns, simplices, cutoff, field_type, ::Val{reps}; next=true)

Compute the ``n``-th intervals of persistent cohomology. The ``n`` is determined from the
`eltype` of `columns`. If `next` is `true`, assemble columns for the next dimension.

Only keep intervals with desired birth/death `cutoff`. Compute homology with coefficients in
`field_type`. If `reps` is `true`, compute representative cocycles.
"""
function nth_intervals(
    filtration, columns::Vector{S}, simplices, cutoff, field_type, ::Val{reps}; next=true,
) where {S<:AbstractSimplex, reps}

    rs = ReductionState{field_type, S}(filtration)
    sizehint!(rs.reduction_matrix, length(columns))
    intervals = compute_intervals!(rs, columns, cutoff, Val(reps))
    if next
        (intervals, assemble_columns!(rs, simplices)...)
    else
        (intervals, nothing, nothing)
    end
end

"""
    ripserer(dists::AbstractMatrix{T}; kwargs...)
    ripserer(points; metric=Euclidean(), births, kwargs...)

Compute the persistent homology of metric space represented by `dists` or `points` and
`metric`. `points` must be an array of bitstypes, such as `NTuple`s or `SVectors`.

# Keyoword Arguments

* `dim_max`: compute persistent homology up to this dimension. Defaults to `1`.
* `modulus`: compute persistent homology with coefficients in the prime field of integers
  mod `modulus`. Defaults to `2`.
* `field_type`: use this type of field of coefficients. Defaults to `Mod{modulus}`.
* `threshold`: compute persistent homology up to diameter smaller than threshold.
  For non-sparse Rips filtrations, it defaults to radius of input space.
* `sparse`: if `true`, use `SparseRipsFiltration`. Defaults to `false`. If the `dists`
  argument is a sparse matrix, it overrides this option.
* `cutoff`: only keep intervals with `persistence(interval) > cutoff`. Defaults to `0`.
* `representatives`: if `true`, return representative cocycles along with persistence
  intervals. Defaults to `false`.
* `metric`: when calculating persistent homology from points, any metric from
  [`Distances.jl`](https://github.com/JuliaStats/Distances.jl) can be used. Defaults to
  `Euclidean()`.
* `births`: when calculating persistent homology from points, births can be used to add
  birth times to vertices. Defaults to all births equal to `0`.
"""
function ripserer(
    dists::AbstractMatrix;
    dim_max=1,
    sparse=false,
    cutoff=0,
    representatives=false,
    modulus=2,
    field_type=Mod{modulus},
    kwargs..., # kwargs for filtration
)
    if sparse || issparse(dists)
        filtration = SparseRipsFiltration(dists; kwargs...)
    else
        filtration = RipsFiltration(dists; kwargs...)
    end
    ripserer(
        filtration;
        dim_max=dim_max,
        representatives=representatives,
        cutoff=cutoff,
        field_type=field_type,
    )
end

function ripserer(points; metric=Euclidean(), births=nothing, kwargs...)
    dists = distances(metric, points, births)
    ripserer(dists; kwargs...)
end

"""
    ripserer(filtration::AbstractFiltration; dim_max=1)

Compute persistent homology from `filtration` object.
"""
function ripserer(
    filtration::AbstractFiltration;
    dim_max=1, representatives=false, cutoff=0, field_type=Mod{2}
)
    ripserer(filtration, cutoff, field_type, Val(dim_max), Val(representatives))
end

function ripserer(filtration, cutoff, field_type, ::Val{0}, ::Val{reps}) where reps
    diagram, _, _ = zeroth_intervals(filtration, cutoff, field_type, Val(reps))
    [diagram]
end

@generated function ripserer(
    filtration, cutoff, field_type, ::Val{dim_max}, ::Val{reps}
) where {dim_max, reps}
    # We unroll the loop over 1:dim_max to ensure type stability.
    # Generated code looks something like:
    #
    # ints_0, cols_1, sxs_1 = zeroth_intervals(filtration, ...)
    # ints_1, cols_2, sxs_2 = nth_itervals(filtration, cols_1, sxs_1, ...)
    # ...
    # ints_dim_max, _, _ = nth_itervals(filtration, cols_dim_max-1, sxs_dim_max-1, ..., next=false)
    #
    # [ints_0, ints_1, ..., ints_dim_max]
    ints = [Symbol("ints_", i) for i in 1:dim_max]
    cols = [Symbol("cols_", i) for i in 1:dim_max]
    sxs = [Symbol("sxs_", i) for i in 1:dim_max]

    expr = quote
        ints_0, $(cols[1]), $(sxs[1]) = zeroth_intervals(
            filtration, cutoff, field_type, Val(reps)
        )
    end
    for i in 1:dim_max-1
        expr = quote
            $expr
            $(ints[i]), $(cols[i+1]), $(sxs[i+1]) = nth_intervals(
                filtration, $(cols[i]), $(sxs[i]), cutoff, field_type, Val(reps)
            )
        end
    end
    quote
        $expr
        $(ints[dim_max]), _, _ =
            nth_intervals(
                filtration, $(cols[dim_max]), $(sxs[dim_max]), cutoff, field_type, Val(reps),
                next=false
            )

        [ints_0, $(ints...)]
    end
end
