"""
    DisjointSetsWithBirth{T}

Almost identical to `DataStructures.IntDisjointSets`, but keeps track of vertex birth times
and birth vertices.
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
        births = collect(zip(births, vertices))
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

function add_interval!(
    intervals, dset::DisjointSetsWithBirth, filtration, vertex, edge, cutoff, ::Val{reps}
) where reps
    birth_time, birth_vertex = birth(dset, vertex)
    death_time = isnothing(edge) ? Inf : birth(edge)
    if death_time - birth_time > cutoff
        birth_simplex = simplex(filtration, Val(0), (birth_vertex,))
        if reps
            rep = (;representative=sort!(
                [simplex(filtration, Val(0), (v,)) for v in find_leaves!(dset, vertex)]
            ))
        else
            rep = NamedTuple()
        end
        push!(intervals, PersistenceInterval(
            birth_time, death_time;
            birth_simplex=birth_simplex,
            death_simplex=edge,
            rep...,
        ))
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
    V = simplex_type(filtration, 0)
    CE = chain_element_type(V, F)
    dset = DisjointSetsWithBirth(vertices(filtration), births(filtration))

    intervals = interval_type(filtration, Val(0), Val(reps), F)[]

    to_skip = simplex_type(filtration, 1)[]
    to_reduce = simplex_type(filtration, 1)[]
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
            # According to the elder rule, the vertex with the higer birth will die first.
            last_vertex = birth(dset, i) > birth(dset, j) ? i : j
            add_interval!(intervals, dset, filtration, last_vertex, edge, cutoff, Val(reps))

            union!(dset, i, j)
            push!(to_skip, edge)
        else
            push!(to_reduce, edge)
        end
        progress && next!(progbar; showvalues=((:n_intervals, length(intervals)),))
    end
    for v in vertices(filtration)
        if find_root!(dset, v) == v && !isnothing(simplex(filtration, Val(0), (v,), 1))
            add_interval!(intervals, dset, filtration, v, nothing, cutoff, Val(reps))
        end
        progress && next!(progbar; showvalues=((:n_intervals, length(intervals)),))
    end
    reverse!(to_reduce)

    thresh = Float64(threshold(filtration))
    return (
        postprocess_diagram(
            filtration,
            PersistenceDiagram(
                sort!(intervals, by=persistence);
                threshold=thresh,
                dim=0,
                field_type=F,
                filtration=filtration,
            ),
        ),
        to_reduce,
        to_skip,
    )
end
