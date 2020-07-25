"""
    DisjointSetsWithBirth{T}

Almost identical to `DataStructures.IntDisjointSets`, but keeps track of vertex birth times.
Has no `num_groups` method.
"""
struct DisjointSetsWithBirth{
    I<:AbstractArray,
    A<:AbstractArray,
    B<:AbstractArray,
    C<:AbstractArray,
}
    vertices::I
    parents::A
    ranks::B
    births::C

    function DisjointSetsWithBirth(vertices, births)
        parents = collect(vertices)
        ranks = similar(parents, Int)
        ranks .= 0
        births = collect(births)
        return new{typeof(vertices), typeof(parents), typeof(ranks), typeof(births)}(
            vertices, parents, ranks, births
        )
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
    leaves = eltype(s.parents)[]
    for i in s.vertices
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
    ranks = s.ranks
    births = s.births
    @inbounds xrank = ranks[x]
    @inbounds yrank = ranks[y]

    if xrank < yrank
        x, y = y, x
    elseif xrank == yrank
        ranks[x] += 1
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
    dset = DisjointSetsWithBirth(vertices(filtration), birth(filtration))
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
        progbar = Progress(
            length(simplices) + n_vertices(filtration);
            desc="Computing 0d intervals... ",
        )
    end
    for edge in simplices
        u, v = vertices(edge)
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
        progress && next!(progbar; showvalues=((:n_intervals, length(intervals)),))
    end
    for v in vertices(filtration)
        if find_root!(dset, v) == v
            push!(intervals, interval(eltype(intervals), dset, filtration, v, nothing, 0))
        end
        progress && next!(progbar; showvalues=((:n_intervals, length(intervals)),))
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
