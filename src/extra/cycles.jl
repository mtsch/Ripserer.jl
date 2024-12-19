"""
    OneSkeleton

A realization of the one-skeleton of a filtered simplicial complex at time t as an
`AbstractGraph`. Some edges of the graph may be removed by placing them into the `removed`
set. This will ignore said edges when looking for neighbours or shortest paths.
Weights are usually simplex birth times, but can be customized.
"""
struct OneSkeleton{T,F<:AbstractFiltration,S<:AbstractCell{1},A} <: AbstractGraph{Int}
    filtration::F
    threshold::T
    weights::A
    removed::Set{S}
end
function OneSkeleton(
    filtration::F,
    thresh::T=threshold(filtration),
    removed=(),
    weights=distance_matrix(filtration),
) where {F,T}
    S = simplex_type(filtration, 1)
    removed = Set{S}(removed)
    return OneSkeleton{T,F,S,typeof(weights)}(filtration, thresh, weights, removed)
end

function Graphs.SimpleGraph(g::OneSkeleton)
    SimpleGraph(edges(g))
end

_birth_or_value(σ::AbstractCell) = birth(σ)
_birth_or_value(σ) = σ
_in(σ::S, g::OneSkeleton{S}) where {S} = !isnothing(σ) && σ ≤ g.threshold && σ ∉ g.removed
_in(σ, g::OneSkeleton) = !isnothing(σ) && birth(σ) ≤ g.threshold && σ ∉ g.removed

# Convert linear indices to filtration indices. Applicable to Cubical because it has
# vertices of type CartesianIndex.
_linear(g, i) = LinearIndices(vertices(g.filtration))[i]
_inv_linear(g, i) = vertices(g.filtration)[i]

Graphs.edgetype(::OneSkeleton) = Edge{Int}

function Graphs.has_edge(g::OneSkeleton, u::Integer, v::Integer)
    if u ∉ vertices(g) || v ∉ vertices(g)
        return false
    else
        σ = simplex(g.filtration, Val(1), _inv_linear.(Ref(g), (v, u)))
        return _in(σ, g)
    end
end

function Graphs.edges(g::OneSkeleton)
    result = edgetype(g)[]
    for sx in Ripserer.edges(g.filtration)
        if _in(sx, g)
            u, v = _linear.(Ref(g), sx)
            push!(result, Edge(u, v))
        end
    end
    return result
end
function Graphs.outneighbors(g::OneSkeleton, u::Integer)
    root = simplex(g.filtration, Val(0), (_inv_linear(g, u),))
    neighbors = Int[]
    for sx in Ripserer.coboundary(g.filtration, root)
        if _in(sx, g)
            v, w = _linear.(Ref(g), sx)
            if v == u
                push!(neighbors, w)
            else
                push!(neighbors, v)
            end
        end
    end
    return neighbors
end
Graphs.inneighbors(g::OneSkeleton, u::Integer) = outneighbors(g, u)

Base.eltype(::OneSkeleton) = Int
Graphs.has_vertex(g::OneSkeleton, u) = 1 ≤ u ≤ nv(g)
Graphs.ne(g::OneSkeleton) = length(edges(g))
Graphs.nv(g::OneSkeleton) = nv(g.filtration)
Graphs.vertices(g::OneSkeleton) = Base.OneTo(nv(g))
Graphs.is_directed(::OneSkeleton) = false
Graphs.is_directed(::Type{<:OneSkeleton}) = false

# need default dist for filtrations...
Graphs.weights(g::OneSkeleton) = g.weights

_heuristic(filtration, src, distances) = dst -> distances[dst, src]
function _path_length(dists, cyc)
    result = zero(eltype(dists))
    for edge in cyc
        i, j = src(edge), dst(edge)
        result += dists[i, j]
    end
    return result
end

"""
    _find_cycle(g, dists)

Find the shortest cycle in `g` that has exactly one edge from `g.removed`.
"""
function _find_cycle(skeleton, dists)
    # Idea:
    # best_weight is the current best length of shortest cycle.
    # best_path is the current best candidate for shortest cycle.
    # best_simplex completes best_path to a cycle.
    # We go through all edges in skeleton.removed and find the shortest path through its
    # endpoints. The shortest among those is the shortest cycle.
    graph = SimpleGraph(skeleton)

    best_weight = missing
    best_path = edgetype(skeleton)[]
    best_sx = first(skeleton.removed)
    for sx in skeleton.removed
        u, v = _linear.(Ref(skeleton), sx)
        path = a_star(graph, u, v, dists, _heuristic(skeleton.filtration, v, dists))
        weight = _path_length(dists, path) + dists[u, v]
        if !isempty(path) && isless(weight, best_weight)
            best_weight = weight
            best_path = path
            best_sx = sx
        end
    end
    if ismissing(best_weight)
        error("no cycle found!")
    else
        result = [best_sx]
        # Convert the type of best_path from `Edge`s to the simplex type of the filtration.
        for e in best_path
            u, v = _inv_linear.(Ref(skeleton), (src(e), dst(e)))
            push!(result, simplex(skeleton.filtration, Val(1), (u, v)))
        end
        return result
    end
end

"""
    reconstruct_cycle(filtration, interval[, t]; distances=distance_matrix(filtration))

Reconstruct the shortest representative cycle for the first homology group of given
`interval`. The optional argument `t` sets the time at which the cycle is to be computed. It
defaults to interval birth time, which gives a cycle similar to a representative cycle
computed from homology. In general, higher times will yield nicer cycles. `t` can be a
simplex or a number.

The optional `distances` keyword argument can be used to change the distance matrix used for
determining edge lengths.

This method uses the representative _co_cycle to compute the cycle. As such, the interval
must be computed with the default cohomology algorithm and must include a representative. To
get such an interval, run `ripserer` with the keyword argument `reps=true` or `reps=1`.

Note that this method only works in the first dimension, as it is based on finding shortest
paths in a graph.

!!! warning
    This feature is still experimental.
"""
function reconstruct_cycle(
    filtration::AbstractFiltration{<:Any,T},
    interval,
    r=birth_simplex(interval);
    distances=distance_matrix(filtration),
) where {T}
    if !hasproperty(interval, :representative)
        throw(
            ArgumentError("interval has no representative! Run `ripserer` with `reps=true`")
        )
    elseif !(eltype(interval.representative) <: AbstractChainElement{<:AbstractCell{1}})
        throw(
            ArgumentError("cycles can only be reconstructed for 1-dimensional intervals.")
        )
    elseif !(birth(interval) ≤ _birth_or_value(r) < death(interval))
        return simplex_type(filtration, 1)[]
    else
        reps = filter!(simplex.(representative(interval))) do sx
            birth(sx) ≤ _birth_or_value(r)
        end
        g = OneSkeleton(filtration, r, reps)
        return _find_cycle(g, distances)
    end
end
