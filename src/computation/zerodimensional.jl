"""
    DisjointSetsWithBirth{T}

Almost identical to `DataStructures.IntDisjointSets`, but keeps track of vertex birth times
and birth vertices.
Has no `num_groups` method.
"""
struct DisjointSetsWithBirth{I<:AbstractArray,A<:AbstractArray,B<:AbstractArray}
    vertices::I
    parents::A
    births::B

    function DisjointSetsWithBirth(vertices, births)
        parents = collect(vertices)
        births = collect(zip(births, vertices))
        return new{typeof(vertices),typeof(parents),typeof(births)}(
            vertices, parents, births
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
    return xroot != yroot ? root_union!(s, xroot, yroot) : xroot
end

function DataStructures.root_union!(s::DisjointSetsWithBirth, x, y)
    parents = s.parents
    births = s.births

    if births[x] > births[y]
        x, y = y, x
    end
    @inbounds parents[y] = x
    @inbounds births[y] = births[x]

    return x
end

birth(dset::DisjointSetsWithBirth, i) = dset.births[i]

function interval(
    dset::DisjointSetsWithBirth, filtration, vertex, parent, edge, cutoff, reps, merge_tree
)
    birth_time, birth_vertex = birth(dset, vertex)
    death_time = isnothing(edge) ? Inf : birth(edge)
    persistence = death_time - birth_time
    # If computing merge tree, we need all intervals, however we don't want to store
    # the representatives of the intervals we'll eventually remove.
    if merge_tree || persistence > cutoff && !isnan(persistence)
        birth_simplex = simplex(filtration, Val(0), (birth_vertex,))
        if reps && (!merge_tree || persistence > cutoff && !isnan(persistence))
            rep = (;
                representative=sort!([
                    simplex(filtration, Val(0), (v,)) for v in find_leaves!(dset, vertex)
                ])
            )
        else
            rep = NamedTuple()
        end
        if merge_tree
            if isnothing(parent)
                parent_simplex = nothing
            else
                parent_simplex = simplex(filtration, Val(0), (parent,))
            end
            children = (; parent_simplex, children=PersistenceInterval[])
        else
            children = NamedTuple()
        end
        meta = (; birth_simplex, death_simplex=edge, rep..., children...)
        return PersistenceInterval(birth_time, death_time, meta)
    else
        return nothing
    end
end

function _build_merge_tree!(diagram, cutoff)
    sort!(diagram; by=x -> (persistence(x), birth_simplex(x)), rev=true)
    filtration = diagram.filtration
    interval_map = Dict(birth_simplex(int) => int for int in diagram)

    # walk top-down adding parents to all intervals and removing ones below the cutoff
    curr_i = 0
    for i in eachindex(diagram)
        int = diagram[i]
        if persistence(int) ≤ cutoff || isnan(persistence(int))
            continue
        end
        if !isnothing(int.parent_simplex)
            parent = interval_map[int.parent_simplex]
            while persistence(parent) ≤ cutoff || isnan(persistence(int))
                parent = interval_map[parent.parent_simplex]
            end
        else
            parent = nothing
        end
        int = PersistenceInterval(int.birth, int.death, (; int.meta..., parent))
        interval_map[int.birth_simplex] = int
        curr_i += 1
        diagram[curr_i] = int
    end
    resize!(diagram.intervals, curr_i)
    reverse!(diagram.intervals)

    # add children to their parent's children list
    for int in diagram
        !isnothing(int.parent) && push!(int.parent.children, int)
    end
    return diagram
end

"""
    zeroth_intervals(filtration, cutoff, verbose, field_type, reps)

Compute 0-dimensional persistent homology using Kruskal's Algorithm.

Only keep intervals with desired birth/death `cutoff`. Compute homology with coefficients in
`field_type`. If `reps` is `true`, compute representative cocycles. Show a progress bar if
`verbose` is set.
"""
function zeroth_intervals(
    filtration, cutoff, verbose, ::Type{F}, reps, merge_tree
) where {F}
    V = simplex_type(filtration, 0)
    CE = chain_element_type(V, F)
    dset = DisjointSetsWithBirth(vertices(filtration), births(filtration))

    intervals = PersistenceInterval[]

    to_skip = simplex_type(filtration, 1)[]
    to_reduce = simplex_type(filtration, 1)[]
    simplices = sort!(edges(filtration))
    if verbose
        progbar = Progress(
            length(simplices) + nv(filtration); desc="Computing 0d intervals... "
        )
    end
    for edge in simplices
        v, u = vertices(edge)
        i = find_root!(dset, u)
        j = find_root!(dset, v)
        if i ≠ j
            # According to the elder rule, the vertex with the higer birth will die first.
            if birth(dset, i) > birth(dset, j)
                parent_vertex = j
                last_vertex = i
            else
                parent_vertex = i
                last_vertex = j
            end
            int = interval(
                dset, filtration, last_vertex, parent_vertex, edge, cutoff, reps, merge_tree
            )
            !isnothing(int) && push!(intervals, int)

            union!(dset, i, j)
            push!(to_skip, edge)
        else
            push!(to_reduce, edge)
        end
        verbose && next!(progbar; showvalues=((:n_intervals, length(intervals)),))
    end
    for v in vertices(filtration)
        if find_root!(dset, v) == v && !isnothing(simplex(filtration, Val(0), (v,)))
            int = interval(dset, filtration, v, nothing, nothing, cutoff, reps, merge_tree)
            push!(intervals, int)
        end
        verbose && next!(progbar; showvalues=((:n_intervals, length(intervals)),))
    end
    reverse!(to_reduce)

    thresh = Float64(threshold(filtration))
    if merge_tree
        diagram = PersistenceDiagram(
            intervals; threshold=thresh, dim=0, field=F, filtration=filtration
        )
        _build_merge_tree!(diagram, cutoff)
    else
        diagram = PersistenceDiagram(
            sort!(intervals; by=persistence);
            threshold=thresh,
            dim=0,
            field=F,
            filtration=filtration,
        )
    end
    return (postprocess_diagram(filtration, diagram), to_reduce, to_skip)
end
