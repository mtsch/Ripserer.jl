"""
    DisjointSetsWithBirth{T}

Almost identical to `DataStructures.IntDisjointSets`, but keeps track of vertex birth times.
Has no `num_groups` method.
"""
struct DisjointSetsWithBirth{T}
    parents ::Vector{Int}
    ranks   ::Vector{Int}
    births  ::Vector{T}

    function DisjointSetsWithBirth(births::AbstractVector{T}) where T
        n = length(births)
        return new{T}(collect(1:n), fill(0, n), copy(births))
    end
end

Base.length(s::DisjointSetsWithBirth) =
    length(s.parents)

function DataStructures.find_root!(s::DisjointSetsWithBirth, x)
    parents = s.parents
    p = parents[x]
    @inbounds if parents[p] != p
        parents[x] = p = find_root!(s, p)
    end
    return p
end

"""
    find_leaves!(s::DisjointSetsWithBirth, x)

Find all leaves below `x`, i.e. vertices that have `x` as root.
"""
function find_leaves!(s::DisjointSetsWithBirth, x)
    leaves = Int[]
    for i in 1:length(s)
        find_root!(s, i) == x && push!(leaves, i)
    end
    return leaves
end

function Base.union!(s::DisjointSetsWithBirth, x, y)
    parents = s.parents
    xroot = find_root!(s, x)
    yroot = find_root!(s, y)
    xroot != yroot ? root_union!(s, xroot, yroot) : xroot
end

function DataStructures.root_union!(s::DisjointSetsWithBirth, x, y)
    parents = s.parents
    rks = s.ranks
    births = s.births
    @inbounds xrank = rks[x]
    @inbounds yrank = rks[y]

    if xrank < yrank
        x, y = y, x
    elseif xrank == yrank
        rks[x] += 1
    end
    @inbounds parents[y] = x
    @inbounds births[x] = min(births[x], births[y])
    return x
end

birth(dset::DisjointSetsWithBirth, i) = dset.births[i]

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
        intervals = PersistenceInterval{Vector{CE}}[]
    else
        intervals = PersistenceInterval{Nothing}[]
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
            push!(intervals, PersistenceInterval(birth(dset, v), Inf, representative))
        end
    end
    reverse!(unreduced)
    progress && printstyled("Assembled $(length(unreduced)) columns.\n", color=:green)
    return (
        sort!(PersistenceDiagram(0, intervals, threshold(filtration))), unreduced, reduced,
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
    result = PersistenceDiagram[]
    zeroth, to_reduce, to_skip = zeroth_intervals(filtration, cutoff, field_type, reps, progress)
    push!(result, zeroth)

    higer_intervals!(result, to_reduce, to_skip, filtration, cutoff, reps, progress, dim_max, field_type)

    return result
end
