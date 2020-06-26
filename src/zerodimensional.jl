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

Base.length(s::DisjointSetsWithBirth) = length(s.parents)

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

function birth_death(dset::DisjointSetsWithBirth, vertex, edge)
    return birth(dset, vertex), isnothing(edge) ? Inf : diam(edge)
end

function interval(
    ::Type{R}, dset::DisjointSetsWithBirth, filtration, vertex, edge, cutoff
) where R<:RepresentativeInterval
    birth, death = birth_death(dset, vertex, edge)
    if death - birth > cutoff
        rep = simplex.(Ref(filtration), Val(0), tuple.(find_leaves!(dset, vertex)))
        birth_simplex = simplex(filtration, Val(0), (findfirst(r -> diam(r) == birth, rep),))
        return R(PersistenceInterval(birth, death), birth_simplex, edge, rep)
    else
        return nothing
    end
end

function interval(
    ::Type{PersistenceInterval}, dset::DisjointSetsWithBirth, _, vertex, edge, cutoff
)
    birth, death = birth_death(dset, vertex, edge)
    if death - birth > cutoff
        return PersistenceInterval(birth, death)
    else
        return nothing
    end
end

"""
    zeroth_intervals(filtration, cutoff, progress, field_type, Val(reps))

Compute 0-dimensional persistent homology using Kruskal's Algorithm.

Only keep intervals with desired birth/death `cutoff`. Compute homology with coefficients in
`field_type`. If `reps` is `true`, compute representative cocycles. Show a progress bar if
`progress` is set.
"""
function zeroth_intervals(
    filtration, cutoff, progress, ::Type{F}, ::Val{reps}
) where {F, reps}
    V = vertex_type(filtration)
    CE = chain_element_type(V, F)
    dset = DisjointSetsWithBirth([birth(filtration, v) for v in 1:n_vertices(filtration)])
    if reps
        intervals = RepresentativeInterval{
            PersistenceInterval,
            vertex_type(filtration),
            Union{edge_type(filtration), Nothing},
            Vector{CE},
        }[]
    else
        intervals = PersistenceInterval[]
    end
    to_skip = edge_type(filtration)[]
    to_reduce = edge_type(filtration)[]
    simplices = sort!(edges(filtration))
    if progress
        progbar = Progress(length(simplices), desc="Computing 0d intervals... ")
    end
    for edge in simplices
        u, v = index.(vertices(edge))
        i = find_root!(dset, u)
        j = find_root!(dset, v)
        if i ≠ j
            # According to the elder rule, the vertex with the lower birth will fall
            # into a later interval.
            v = birth(dset, i) > birth(dset, j) ? i : j
            int = interval(eltype(intervals), dset, filtration, v, edge, cutoff)
            if !isnothing(int)
                push!(intervals, int)
            end
            union!(dset, i, j)
            push!(to_skip, edge)
        else
            push!(to_reduce, edge)
        end
        progress && next!(progbar)
    end
    for v in 1:n_vertices(filtration)
        if find_root!(dset, v) == v
            push!(intervals, interval(eltype(intervals), dset, filtration, v, nothing, 0))
        end
    end
    reverse!(to_reduce)
    progress && printstyled(stderr, "Assembled $(length(to_reduce)) edges.\n", color=:green)
    return (
        sort!(PersistenceDiagram(
            0,
            map(int -> postprocess_interval(filtration, int), intervals),
            threshold(filtration)
        )),
        to_reduce,
        to_skip,
    )
end
