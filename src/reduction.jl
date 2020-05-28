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
        working_column = Column{CE}()
        reduction_entries = Column{SE}()

        new{Field, S, SE, C, CE, F}(
            filtration, ReductionMatrix{C, SE}(), working_column, reduction_entries
        )
    end
end

simplex_type(rs::ReductionState{<:Any, S}) where S = S
coface_element(rs::ReductionState{<:Any, <:Any, <:Any, <:Any, CE}) where CE = CE

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
    return pivot(rs.working_column)
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
    return pivot(rs.working_column)
end

"""
    reduce_working_column!(rs::ReductionState, column_simplex, cutoff, reps)

Reduce the working column by adding other columns to it until it has the lowest pivot or is
reduced. Record it in the reduction matrix and return the persistence intervals with
persistence larger than `cutoff`. If `reps` is `true`, add representative cocycles to the
interval.
"""
function reduce_working_column!(
    rs::ReductionState{F, S}, column_simplex, cutoff, reps
) where {F, S}
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
        return PersistenceInterval(birth, death, chain_element_type(S, F)[])
    elseif reps && death - birth > cutoff
        representative = collect(rs.reduction_matrix[current_pivot])
        return PersistenceInterval(birth, death, representative)
    elseif !isfinite(death) || death - birth > cutoff
        return PersistenceInterval(birth, death)
    else
        return nothing
    end
end

"""
    compute_pairs!(rs::ReductionState, columns, cutoff, reps, progress)

Compute persistence intervals by reducing `columns`, a collection of simplices. Return
`PersistenceDiagram`. Only keep intervals with persistence larger than `cutoff` and return
representative cocycles if `reps` is `true`. If `progress` is set, show a progress bar.
"""
function compute_intervals!(
    rs::ReductionState{F, S}, columns, cutoff, reps, progress
) where {S, F}
    T = dist_type(rs.filtration)
    if reps
        intervals = PersistenceInterval{T, Vector{chain_element_type(S, F)}}[]
    else
        intervals = PersistenceInterval{T, Nothing}[]
    end
    if progress
        progbar = Progress(length(columns), desc="Computing $(dim(S))d intervals... ")
    end
    for column in columns
        interval = reduce_working_column!(rs, column, cutoff, reps)
        if !isnothing(interval)
            push!(intervals, interval)
        end
        progress && next!(progbar)
    end
    return sort!(
        PersistenceDiagram(dim(eltype(columns)), threshold(rs.filtration), intervals)
    )
end

"""
    assemble_columns!(rs::ReductionState, unreduced, reduced, progress)

Assemble columns that need to be reduced in the next dimension. Apply clearing optimization.
"""
function assemble_columns!(
    rs::ReductionState{<:Any, S}, unreduced, reduced, progress
) where S
    C = coface_type(S)
    new_unreduced = C[]
    new_reduced = C[]

    if progress
        progbar = Progress(
            length(unreduced) + length(reduced), desc="Assembling columns...     "
        )
    end
    for cols in (unreduced, reduced)
        for simplex in cols
            for coface in coboundary(rs.filtration, simplex, Val(false))
                if !has_column(rs.reduction_matrix, coface)
                    push!(new_unreduced, abs(coface))
                else
                    push!(new_reduced, abs(coface))
                end
            end
            progress && next!(progbar)
        end
    end
    sort!(new_unreduced, rev=true)
    progress && printstyled("Assembled $(length(new_unreduced)) columns.\n", color=:green)
    return new_unreduced, new_reduced
end

"""
    nth_intervals(filtration, columns, simplices, cutoff, field_type, ::Val{reps}; next=true)

Compute the ``n``-th intervals of persistent cohomology. The ``n`` is determined from the
`eltype` of `columns`. If `assemble` is `true`, assemble columns for the next dimension.

Only keep intervals with persistence larger than `cutoff`.  Compute homology with
coefficients in `field_type`. If `reps` is `true`, compute representative cocycles. Show a
progress bar if `progress` is set.
"""
function nth_intervals(
    filtration, unreduced, reduced, cutoff, field_type, reps, assemble, progress
)
    rs = ReductionState{field_type, eltype(unreduced)}(filtration)
    sizehint!(rs.reduction_matrix, length(unreduced))
    intervals = compute_intervals!(rs, unreduced, cutoff, reps, progress)
    if assemble
        return (intervals, assemble_columns!(rs, unreduced, reduced, progress)...)
    else
        return (intervals, nothing, nothing)
    end
end

"""
    zeroth_representative(dset, vertex, reps, CE, V)

Collect the zero dimensional representative of `vertex`.
"""
function zeroth_representative(filtration, dset, vertex, reps, CE, V)
    if reps
        return map(find_leaves!(dset, vertex)) do u
            CE(V((u,), birth(filtration, u)))
        end
    else
        return nothing
    end
end

"""
    zeroth_intervals(filtration, cutoff, field_type, reps, progress)

Compute 0-dimensional persistent homology using Kruskal's Algorithm.

Only keep intervals with desired birth/death `cutoff`. Compute homology with coefficients in
`field_type`. If `reps` is `true`, compute representative cocycles. Show a progress bar if
`progress` is set.
"""
function zeroth_intervals(filtration, cutoff, field_type, reps, progress)
    T = dist_type(filtration)
    V = vertex_type(filtration)
    CE = chain_element_type(V, field_type)
    dset = DisjointSetsWithBirth([birth(filtration, v) for v in 1:n_vertices(filtration)])
    if reps
        intervals = PersistenceInterval{T, Vector{CE}}[]
    else
        intervals = PersistenceInterval{T, Nothing}[]
    end
    reduced = edge_type(filtration)[]
    unreduced = edge_type(filtration)[]
    simplices = edges(filtration)
    if progress
        progbar = Progress(length(simplices), desc="Computing 0d intervals... ")
    end
    for sx in simplices
        u, v = vertices(sx)
        i = find_root!(dset, u)
        j = find_root!(dset, v)
        if i ≠ j
            # According to the elder rule, the vertex with the lower birth will fall
            # into a later interval.
            dead = birth(dset, i) > birth(dset, j) ? i : j
            if diam(sx) - birth(dset, dead) > cutoff
                representative = zeroth_representative(filtration, dset, dead, reps, CE, V)
                interval = PersistenceInterval(birth(dset, dead), diam(sx), representative)
                push!(intervals, interval)
            end
            union!(dset, i, j)
            push!(reduced, sx)
        else
            push!(unreduced, sx)
        end
        progress && next!(progbar)
    end
    for v in 1:n_vertices(filtration)
        if find_root!(dset, v) == v
            representative = zeroth_representative(filtration, dset, v, reps, CE, V)
            push!(intervals, PersistenceInterval(birth(dset, v), ∞, representative))
        end
    end
    reverse!(unreduced)
    progress && printstyled("Assembled $(length(unreduced)) columns.\n", color=:green)
    return (
        sort!(PersistenceDiagram(0, threshold(filtration), intervals)), unreduced, reduced
    )
end

"""
    ripserer(dists::AbstractMatrix; kwargs...)
    ripserer(points; metric=Euclidean(), births, kwargs...)
    ripserer(filtration::AbstractFiltration; kwargs...)

Compute the persistent homology of metric space represented by `dists`, `points` and
`metric` or an `::AbstractFiltration`.

If using points, `points` must be an array of bitstypes, such as `NTuple`s or `SVectors`.

# Keyoword Arguments

* `dim_max`: compute persistent homology up to this dimension. Defaults to `1`.
* `modulus`: compute persistent homology with coefficients in the prime field of integers
  mod `modulus`. Defaults to `2`.
* `field_type`: use this type of field of coefficients. Defaults to `Mod{modulus}`.
* `threshold`: compute persistent homology up to diameter smaller than threshold.
  For non-sparse Rips filtrations, it defaults to radius of input space.
* `cutoff`: only keep intervals with `persistence(interval) > cutoff`. Defaults to `0`.
* `representatives`: if `true`, return representative cocycles along with persistence
  intervals. Defaults to `false`.
* `progress`: If `true`, show a progress bar. Defaults to `false`.
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
    progress=false,
    kwargs..., # kwargs for filtration
)
    if issparse(dists)
        filtration = SparseRips(dists; kwargs...)
    else
        filtration = Rips(dists; kwargs...)
    end
    return ripserer(
        filtration;
        dim_max=dim_max,
        representatives=representatives,
        cutoff=cutoff,
        field_type=field_type,
        progress=progress,
    )
end

function ripserer(points; metric=Euclidean(), births=nothing, kwargs...)
    dists = distances(metric, points, births)
    return ripserer(dists; kwargs...)
end

function ripserer(
    filtration::AbstractFiltration;
    dim_max=1, representatives=false, cutoff=0, field_type=Mod{2}, progress=false
)
    return ripserer(filtration, cutoff, field_type, Val(dim_max), representatives, progress)
end

function ripserer(
    filtration, cutoff, field_type, ::Val{dim_max}, reps, progress
) where dim_max
    res = PersistenceDiagram[]
    res_0, cols, sxs = zeroth_intervals(filtration, cutoff, field_type, reps, progress)
    push!(res, res_0)
    for dim in 1:dim_max
        res_n, cols, sxs = nth_intervals(
            filtration, cols, sxs, cutoff, field_type, reps, dim ≠ dim_max, progress)
        push!(res, res_n)
        # TODO: A lot of things get allocated and then destroyed in nth_intervals, so it
        # might make sense to do a GC.gc() here.
    end
    return res
end
