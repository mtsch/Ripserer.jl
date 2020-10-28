function _get_edge_neighbors!(buffer, graph::AbstractGraph, edge)
    nu = neighbors(graph, src(edge))
    nv = neighbors(graph, dst(edge))
    empty!(buffer)
    # Sorted version of intersect
    i, j = 1, 1
    while i ≤ length(nu) && j ≤ length(nv)
        u = nu[i]
        v = nv[j]
        if u == v
            push!(buffer, u)
            i += 1; j += 1
        elseif u < v
            i += 1
        elseif v < u
            j += 1
        end
    end
    return buffer
end

function _get_adjacent_edges!(edgeset, graph::AbstractGraph, edge)
    for u in (src(edge), dst(edge))
        for v in neighbors(graph, u)
            push!(edgeset, Edge(TupleTools.sort((u, v))...))
        end
    end
    return edgeset
end

function _is_dominated(buffer, graph::AbstractGraph, edge)
    # N[e] ⊆ N[v] for some v ∉ edge
    s, d = src(edge), dst(edge)
    ne = _get_edge_neighbors!(buffer, graph, edge)
    for w in ne
        (w == s || w == d) && continue

        nw = neighbors(graph, w)
        # Sorted version of issubset
        i, j = 1, 1
        while i ≤ length(ne) && j ≤ length(nw)
            u = ne[i]
            v = nw[j]
            if u < v
                break
            elseif u == v
                i += 1
            end
            j += 1
        end
        i > length(ne) && return true
    end
    return false
end

function _add_core_edge!(core_edges, edge, val)
    s, d = src(edge), dst(edge)
    core_edges[s, d] = core_edges[d, s] = val
end

function _back_pass!(buffer, core_edges, neighbor_edges, graph, parent, edges, i, time)
    empty!(neighbor_edges)
    _get_adjacent_edges!(neighbor_edges, graph, parent)
    graph′ = copy(graph)
    # Some edges that were not dominated before may be dominated now.
    for j in (i - 1):-1:1
        edge, _ = edges[j]
        if iszero(core_edges[src(edge), dst(edge)])
            if edge ∈ neighbor_edges && !_is_dominated(buffer, graph′, edge)
                _add_core_edge!(core_edges, edge, time)
                _get_adjacent_edges!(neighbor_edges, graph′, edge)
            else
                rem_edge!(graph′, edge)
            end
        end
    end
end

function _core_graph(filtration::Rips{<:Any, T}; eps=0, progress=false) where T
    #eps ≠ 0 && error("eps not implemented")
    filtration_edges = Tuple{Edge{Int}, T}[]
    for e in sort!(edges(filtration))
        v, u = vertices(e)
        push!(filtration_edges, (Edge{Int}(u, v), birth(e)))
    end

    buffer = Int[]
    core_edges = zeros(T, nv(filtration), nv(filtration))
    neighbor_edges = Set{Edge{Int}}()
    graph = Graph(nv(filtration))
    # Add self loops so i ∈ neighbors(graph, i)
    for i in 1:nv(filtration)
        add_edge!(graph, i, i)
    end

    prev_time = filtration_edges[1][2]
    if progress
        progbar = Progress(length(filtration_edges), desc="Collapsing edges...")
    end
    for (i, (edge, curr_time)) in enumerate(filtration_edges)
        add_edge!(graph, edge)
        was_dominated = false
        if !_is_dominated(buffer, graph, edge)
            _add_core_edge!(core_edges, edge, curr_time)
            was_dominated = true
        end
        # This part implements the approximate computation mentioned in reference.
        # We need to look at the next time for step, because we want the error to be smaller
        # than eps.
        step = filtration_edges[min(i + 1, end)][2] - prev_time
        if (iszero(eps) && was_dominated) || step > eps
            _back_pass!(
                buffer,
                core_edges,
                neighbor_edges,
                graph,
                edge,
                filtration_edges,
                i,
                curr_time,
            )
            prev_time = curr_time
        end
        progress && next!(progbar)
    end
    return sparse(core_edges)
end

"""
    EdgeCollapsedRips{I, T} <: AbstractRipsFiltration{I, T}

Perform a sequence of edge collapses on a filtration. This may significantly reduce
computation time and does not change the result. The speedup is especially apparent with
datasets that have a boundary, and with high-dimensional persistent homology computation.

The drawback is that doing the collapses themselves can be time-consuming. The construction
does not require a lot of memory (``\mathcal{O}(n^2)`` for ``n`` vertices). This might still
be a good choice for large inputs if you are willing to wait but don't have enough memory
to compute persistent homology with `Rips`.

Setting the `eps` keyword argument approximates the collapses instead. For small `eps`, the
bottleneck distance between the exact and approximate diagrams is guaranteed to be at most
`eps`. For large `eps`, new infinite intervals may appear.

See the reference below for a description of the algorithm.

# Constructors

* `EdgeCollapsedRips(::AbstractRipsFiltration; progress=false, threshold=nothing, eps=0)`:
  Collapse a given filtration. Setting `progress` shows a progress bar.

* `EdgeCollapsedRips(::EdgeCollapsedRips; progress=false, threshold=nothing, eps=0)`: Allows
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

julia> ripserer(collapsed; dim_max=5)
6-element Vector{PersistenceDiagramsBase.PersistenceDiagram}:
 100-element 0-dimensional PersistenceDiagram
 75-element 1-dimensional PersistenceDiagram
 37-element 2-dimensional PersistenceDiagram
 14-element 3-dimensional PersistenceDiagram
 1-element 4-dimensional PersistenceDiagram
 0-element 5-dimensional PersistenceDiagram

```

# Reference

Boissonnat, J. D., & Pritam, S. (2019). [Edge Collapse and Persistence of Flag
Complexes](https://hal.inria.fr/hal-02395227/).
"""
struct EdgeCollapsedRips{I,T} <: AbstractRipsFiltration{I,T}
    adj::SparseMatrixCSC{T,Int}
    threshold::T
    eps::T
end

function EdgeCollapsedRips(
    rips::AbstractRipsFiltration{I,T}; progress=false, threshold=nothing, eps=zero(T)
) where {I,T}
    adj = _core_graph(rips; progress=progress, eps=T(eps))
    if isnothing(threshold)
        threshold = maximum(adj)
    end
    return EdgeCollapsedRips{I,T}(adj, T(threshold), T(eps))
end
function EdgeCollapsedRips{I}(dist_or_points; progress=false, eps=0, kwargs...) where {I}
    return EdgeCollapsedRips(Rips{I}(dist_or_points; kwargs...); progress=progress, eps=eps)
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
