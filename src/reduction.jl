# reduction matrix ======================================================================= #
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
simplex_element(rs::ReductionState{<:Any, <:Any, SE}) where SE =
    SE
coface_element(rs::ReductionState{<:Any, <:Any, <:Any, <:Any, CE}) where CE =
    CE

"""
    add!(rs::ReductionState, index)

Add column with column `index` multiplied by the correct factor to `rs.working_column`.
Also record the addition in `rs.reduction_entries`.
"""
function add!(rs::ReductionState, current_pivot)
    λ = -coef(current_pivot)
    for element in rs.reduction_matrix[current_pivot]
        push!(rs.reduction_entries, λ * element)
        for coface in coboundary(rs.filtration, simplex(element))
            push!(rs.working_column, coface_element(rs)(coface, λ * coef(element)))
        end
    end
    pivot(rs.working_column)
end

"""
    initialize!(rs::ReductionState, column_simplex)

Initialize `rs.working_column` by emptying it and `reduction_entries` and pushing the
coboundary of `column_simplex` to `rs.working_column`.
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
    reduce_working_column!(rs::ReductionState, res, column_simplex)

Reduce the working column by adding other columns to it until it has the lowest pivot or is
reduced. Record it in the reduction matrix and return the persistence interval.
"""
function reduce_working_column!(
    rs::ReductionState,
    column_simplex::AbstractSimplex,
    ::Val{representatives},
) where representatives
    current_pivot = initialize!(rs, column_simplex)

    while !isnothing(current_pivot) && has_column(rs.reduction_matrix, current_pivot)
        current_pivot = add!(rs, current_pivot)
    end
    if isnothing(current_pivot)
        death = ∞
    else
        insert_column!(rs.reduction_matrix, current_pivot)
        push!(rs.reduction_entries, column_simplex)
        move!(rs.reduction_matrix, rs.reduction_entries, times=inv(coef(current_pivot)))
        death = diam(simplex(current_pivot))
    end
    birth = diam(column_simplex)
    if representatives
        if !isnothing(current_pivot)
            representative = representatives(rs.reduction_matrix, current_pivot, death)
            PersistenceInterval(birth, death, representative)
        else
            PersistenceInterval(birth, death, eltype(rs.reduction_matrix)[])
        end
    else
        PersistenceInterval(birth, death)
    end
end

"""
    compute_pairs!(rs::ReductionState, columns)

Compute persistence intervals by reducing `columns`, a collection of simplices.
"""
function compute_intervals!(
    rs::ReductionState, columns, ratio, ::Val{representatives},
) where representatives
    T = dist_type(rs.filtration)
    if representatives
        SE = simplex_element(rs)
        intervals = PersistenceInterval{T, Vector{SE}}[]
    else
        intervals = PersistenceInterval{T, Nothing}[]
    end
    for column in columns
        interval = reduce_working_column!(rs, column, Val(representatives))
        if death(interval) > birth(interval) * ratio
            push!(intervals, interval)
        end
    end
    sort!(PersistenceDiagram(dim(eltype(columns)), intervals))
end

"""
    assemble_columns!(rs::ReductionState, columns, simplices)

Assemble columns that need to be reduced in the next dimension. Apply clearing optimization.
"""
function assemble_columns!(rs::ReductionState, simplices)
    S = coface_type(simplex_type(rs))
    columns = S[]
    new_simplices = S[]

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
    zeroth_intervals(filtration)

Compute 0-dimensional persistent homology using Kruskal's Algorithm.
If `filtration` is sparse, also return a vector of all 1-simplices with diameter below
threshold.
"""
function zeroth_intervals(filtration, ratio, field_type, ::Val{representatives}) where representatives
    T = dist_type(filtration)
    V = vertex_type(filtration)
    dset = DisjointSetsWithBirth([birth(filtration, v) for v in 1:n_vertices(filtration)])
    if representatives
        intervals = PersistenceInterval{T, Vector{V}}[]
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
            if representatives
                representative = map(
                    x -> V((x,), birth(filtration, x)) => one(field_type), find_leaves!(dset, i),
                )
            else
                representative = nothing
            end
            interval = PersistenceInterval(
                max(birth(dset, i), birth(dset, j)), diam(sx), representative,
            )
            if death(interval) > birth(interval) * ratio
                push!(intervals, interval)
            end
            union!(dset, i, j)
        else
            push!(columns, sx)
        end
    end
    for v in 1:n_vertices(filtration)
        if find_root!(dset, v) == v
            if representatives
                representative = map(x -> V((x,), birth(dset, x)), find_leaves!(dset, v))
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
    nth_intervals(filtration, columns, simplices; next=true)

Compute the ``n``-th intervals of persistent cohomology. The ``n`` is determined from the
`eltype` of `columns`. If `next` is `true`, assemble columns for the next dimension.
"""
function nth_intervals(
    filtration, columns::Vector{S}, simplices, ratio, field_type, ::Val{representatives}; next=true,
) where {S<:AbstractSimplex, representatives}

    rs = ReductionState{field_type, S}(filtration)
    sizehint!(rs.reduction_matrix, length(columns))
    intervals = compute_intervals!(rs, columns, ratio, Val(representatives))
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
* `field_type`: use this type of field of coefficients. Defaults to `PrimeField{modulus}`.
* `threshold`: compute persistent homology up to diameter smaller than threshold.
  For non-sparse Rips filtrations, it defaults to radius of input space.
* `sparse`: if `true`, use `SparseRipsFiltration`. Defaults to `false`. If the `dists`
  argument is a sparse matrix, it overrides this option.
* `ratio`: only keep intervals with `death(interval) > birth(interval) * ratio`.
  Defaults to `1`.
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
    ratio=1,
    representatives=false,
    modulus=2,
    field_type=PrimeField{modulus},
    kwargs...,
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
        ratio=ratio,
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
ripserer(filtration::AbstractFiltration;
         dim_max=1, representatives=false, ratio=1, field_type=PrimeField{2}) =
    ripserer(filtration, ratio, field_type, Val(dim_max), Val(representatives))

function ripserer(filtration::AbstractFiltration, ratio, field_type, ::Val{0}, ::Val{C}) where C
    diagram, _, _ = zeroth_intervals(filtration, ratio, field_type, Val(C))
    [diagram]
end

@generated function ripserer(
    filtration::AbstractFiltration, ratio, field_type, ::Val{D}, ::Val{C},
) where {D, C}
    # We unroll the loop over 1:D to ensure type stability.
    # Generated code looks something like:
    #
    # ints_0, cols_1, sxs_1 = zeroth_intervals(filtration, ...)
    # ints_1, cols_2, sxs_2 = nth_itervals(filtration, cols_1, sxs_1, ...)
    # ...
    # ints_D, _, _ = nth_itervals(filtration, cols_D-1, sxs_D-1, ..., next=false)
    #
    # [ints_0, ints_1, ..., ints_D]
    ints = [Symbol("ints_", i) for i in 1:D]
    cols = [Symbol("cols_", i) for i in 1:D]
    sxs = [Symbol("sxs_", i) for i in 1:D]

    expr = quote
        ints_0, $(cols[1]), $(sxs[1]) = zeroth_intervals(filtration, ratio, field_type, Val(C))
    end
    for i in 1:D-1
        expr = quote
            $expr
            $(ints[i]), $(cols[i+1]), $(sxs[i+1]) =
                nth_intervals(filtration, $(cols[i]), $(sxs[i]), ratio, field_type, Val(C))
        end
    end
    quote
        $expr
        $(ints[D]), _, _ =
            nth_intervals(filtration, $(cols[D]), $(sxs[D]), ratio, field_type, Val(C), next=false)

        [ints_0, $(ints...)]
    end
end
