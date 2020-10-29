struct LazyFalses end
Base.getindex(::LazyFalses, i, j) = false

function _get_edge_neighbors!(buffer, graph::AbstractGraph, edge, removed=LazyFalses())
    s, d = src(edge), dst(edge)
    nu = neighbors(graph, s)
    nv = neighbors(graph, d)
    empty!(buffer)
    # Sorted version of intersect
    i, j = 1, 1
    while i ≤ length(nu) && j ≤ length(nv)
        u = nu[i]
        v = nv[j]
        if u == v
            if !removed[u, s] && !removed[u, d]
                push!(buffer, u)
            end
            i += 1; j += 1
        elseif u < v
            i += 1
        elseif v < u
            j += 1
        end
    end
    return buffer
end

function _get_adjacent_edges!(adjacent, graph::AbstractGraph, edge, removed=LazyFalses())
    for u in (src(edge), dst(edge))
        for v in neighbors(graph, u)
            if !removed[v, u]
                adjacent[u,v] = adjacent[v,u] = true
            end
        end
    end
    return adjacent
end

function _not_dominated(buffer, graph::AbstractGraph, edge, removed=LazyFalses())
    # N[e] ⊆ N[v] for some v ∉ edge
    s, d = src(edge), dst(edge)
    ne = _get_edge_neighbors!(buffer, graph, edge, removed)
    for w in ne
        (w == s || w == d) && continue

        nw = neighbors(graph, w)
        # Sorted version of issubset
        i, j = 1, 1
        while i ≤ length(ne) && j ≤ length(nw)
            u = ne[i]
            v = nw[j]
            if removed[v, w]
                #
            elseif u < v
                break
            elseif u == v
                i += 1
            end
            j += 1
        end
        i > length(ne) && return false
    end
    return true
end

function _add_core_edge!(core_edges, edge, val)
    s, d = src(edge), dst(edge)
    core_edges[s, d] = core_edges[d, s] = val
end

function _remove_edge!(removed, edge)
    s, d = src(edge), dst(edge)
    removed[s, d] = removed[d, s] = true
end

function _in(edgeset, edge::Edge)
    return edgeset[src(edge), dst(edge)]
end

function _core_graph(filtration::Rips{<:Any, T}; progress=false) where T
    filtration_edges = Tuple{Edge{Int}, T}[]
    for e in sort!(edges(filtration))
        v, u = vertices(e)
        push!(filtration_edges, (Edge{Int}(u, v), birth(e)))
    end

    buffer = Int[]
    removed = fill(false, nv(filtration), nv(filtration))
    core = zeros(T, nv(filtration), nv(filtration))
    adjacent = fill(false, nv(filtration), nv(filtration))
    graph = Graph(nv(filtration))
    # Add self loops so i ∈ neighbors(graph, i)
    for i in 1:nv(filtration)
        add_edge!(graph, i, i)
    end

    prev_time = filtration_edges[1][2]
    if progress
        progbar = Progress(length(filtration_edges), desc="Collapsing edges...")
    end
    for (i, (curr_edge, curr_time)) in enumerate(filtration_edges)
        add_edge!(graph, curr_edge)
        if _not_dominated(buffer, graph, curr_edge)
            _add_core_edge!(core, curr_edge, curr_time)

            adjacent .= false
            removed .= false
            _get_adjacent_edges!(adjacent, graph, curr_edge)
            # Some edges that were not dominated before may be dominated now.
            for j in (i - 1):-1:1
                edge, _ = filtration_edges[j]
                if iszero(core[src(edge), dst(edge)])
                    if _in(adjacent, edge) && _not_dominated(buffer, graph, edge, removed)
                        _add_core_edge!(core, edge, curr_time)
                        _get_adjacent_edges!(adjacent, graph, edge, removed)
                    else
                        _remove_edge!(removed, edge)
                    end
                end
            end
        end
        progress && next!(progbar)
    end
    return sparse(core)
end

"""
    EdgeCollapsedRips{I, T} <: AbstractRipsFiltration{I, T}

Perform a sequence of edge collapses on a filtration. This may significantly reduce
computation time and does not change the result. The speedup is especially apparent with
datasets that have a boundary, and with high-dimensional persistent homology computation.

The drawback is that doing the collapses themselves can be time-consuming. The construction
does not require a lot of memory (``\\mathcal{O}(n^2)`` for ``n`` vertices). This might
still be a good choice for large inputs if you are willing to wait but don't have enough
memory to compute persistent homology with `Rips`.

See the reference below for a description of the algorithm.

# Constructors

* `EdgeCollapsedRips(::AbstractRipsFiltration; progress=false, threshold=nothing )`:
  Collapse a given filtration. Setting `progress` shows a progress bar.

* `EdgeCollapsedRips(::EdgeCollapsedRips; progress=false, threshold=nothing)`: Allows
  changing `I` or `threshold` without recomputing.

* `EdgeCollapsedRips(arg; kwargs...)`: Use `arg` and `kwargs` to construct a `Rips`
  filtration and collapse it

* `EdgeCollapsedRips{I}(arg; kwargs...)`: Change the index type used to represent
  simplices. May be necessary for large inputs and high-dimensional computation.

# Examples

```jldoctest
julia> using Random; Random.seed!(1337);

julia> data = [tuple(rand(6)...) for _ in 1:100];

julia> rips = Rips(data)
Rips{Int64, Float64}(nv=100, sparse=false)

julia> length(Ripserer.edges(rips))
3906

julia> collapsed = EdgeCollapsedRips(data) # or EdgeCollapsedRips(rips)
EdgeCollapsedRips{Int64, Float64}(nv=100)

julia> length(Ripserer.edges(collapsed))
1419

julia> ripserer(rips) == ripserer(collapsed)
true

julia> ripserer(collapsed; dim_max=4)
5-element Vector{PersistenceDiagramsBase.PersistenceDiagram}:
 100-element 0-dimensional PersistenceDiagram
 75-element 1-dimensional PersistenceDiagram
 37-element 2-dimensional PersistenceDiagram
 14-element 3-dimensional PersistenceDiagram
 1-element 4-dimensional PersistenceDiagram

```

# Reference

Boissonnat, J. D., & Pritam, S. (2019). [Edge Collapse and Persistence of Flag
Complexes](https://hal.inria.fr/hal-02395227/).
"""
struct EdgeCollapsedRips{I,T} <: AbstractRipsFiltration{I,T}
    adj::SparseMatrixCSC{T,Int}
    threshold::T
end

function EdgeCollapsedRips(
    rips::AbstractRipsFiltration{I,T}; progress=false, threshold=nothing
) where {I,T}
    adj = _core_graph(rips; progress=progress)
    if isnothing(threshold)
        threshold = maximum(adj)
    end
    return EdgeCollapsedRips{I,T}(adj, T(threshold))
end
function EdgeCollapsedRips{I}(dist_or_points; progress=false, kwargs...) where {I}
    return EdgeCollapsedRips(Rips{I}(dist_or_points; kwargs...); progress=progress)
end
function EdgeCollapsedRips{I}(
    rips::EdgeCollapsedRips{<:Any,T}; threshold=rips.threshold
) where {I,T}
    return EdgeCollapsedRips{I,T}(rips.adj, T(threshold))
end
function EdgeCollapsedRips(arg; kwargs...)
    return EdgeCollapsedRips{Int}(arg; kwargs...)
end

threshold(rips::EdgeCollapsedRips) = rips.threshold
adjacency_matrix(rips::EdgeCollapsedRips) = rips.adj
