struct OneSkeleton{T, F<:AbstractFiltration, S<:AbstractSimplex{1}} <: AbstractGraph{Int}
    filtration::F
    threshold::T
    removed::Set{S}
end
function OneSkeleton(
    filtration::F, thresh::T=threshold(filtration), removed=()
) where {F, T}
    S = simplex_type(filtration, 1)
    removed = Set{S}(removed)
    return OneSkeleton{T, F, S}(filtration, thresh, removed)
end

_birth_or_value(σ::AbstractSimplex) = birth(σ)
_birth_or_value(σ) = σ
_in(σ::S, g::OneSkeleton{S}) where S = !isnothing(σ) && σ ≤ g.threshold && σ ∉ g.removed
_in(σ, g::OneSkeleton) = !isnothing(σ) && birth(σ) ≤ g.threshold && σ ∉ g.removed

LightGraphs.edgetype(::OneSkeleton) = Edge{Int}

function LightGraphs.has_edge(g::OneSkeleton, u::Integer, v::Integer)
    if u ∉ vertices(g) || v ∉ vertices(g)
        return false
    else
        σ = simplex(g.filtration, Val(1), _inv_linear.(Ref(g), (v, u)))
        return _in(σ, g)
    end
end

function LightGraphs.edges(g::OneSkeleton)
    result = edgetype(g)[]
    for sx in Ripserer.edges(g.filtration)
        if _in(sx, g)
            u, v = _linear.(Ref(g), sx)
            push!(result, Edge(u, v))
        end
    end
    return result
end
function LightGraphs.outneighbors(g::OneSkeleton, u::Integer)
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
LightGraphs.inneighbors(g::OneSkeleton, u::Integer) = outneighbors(g, u)

Base.eltype(::OneSkeleton) = Int
LightGraphs.has_vertex(g::OneSkeleton, u) = 1 ≤ u ≤ nv(g)
LightGraphs.ne(g::OneSkeleton) = length(edges(g))
LightGraphs.nv(g::OneSkeleton) = nv(g.filtration)
LightGraphs.vertices(g::OneSkeleton) = Base.OneTo(nv(g))
LightGraphs.is_directed(::OneSkeleton) = false
LightGraphs.is_directed(::Type{<:OneSkeleton}) = false

# need default dist for filtrations...
LightGraphs.weights(g::OneSkeleton) = dist(g.filtration)

_heuristic(filtration, src) = dst -> dist(filtration, dst, src)
function _path_length(dists, cyc)
    result = zero(eltype(dists))
    for edge in cyc
        i, j = src(edge), dst(edge)
        result += dists[i, j]
    end
    return result
end

# Convert linear indices to filtration indices. Applicable to Cubical.
_linear(g, i) = LinearIndices(vertices(g.filtration))[i]
_inv_linear(g, i) = vertices(g.filtration)[i]

function _find_cycle(g)
    best_weight = missing
    best_path = edgetype(g)[]
    best_sx = first(g.removed)
    dists = dist(g.filtration)
    for sx in g.removed
        u, v = _linear.(Ref(g), sx)
        path = a_star(g, u, v, dists, _heuristic(g.filtration, u))
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
        for e in best_path
            u, v = _inv_linear.(Ref(g), (src(e), dst(e)))
            push!(result, simplex(g.filtration, Val(1), (u, v)))
        end
        return result
    end
end

"""
    reconstruct_cycle(filtration, interval[, t])

Reconstruct the shortest representative cycle for the first homology group of given
`interval`. The optional argument `t` sets the time at which the cycle is to be computed. It
defaults to interval birth time, which gives a cycle similar to a representative cycle
computed from homology. In general, higher times will yield nicer cycles.

This method uses the representative cocycle to compute the cycle. As such, the interval must
include a representative. To get such an interval, run `ripserer` with the keyword argument
`reps=true`.
"""
function reconstruct_cycle(
    filtration::AbstractFiltration{<:Any, T}, interval, r=birth_simplex(interval)
) where T
    if !hasproperty(interval, :representative)
        throw(ArgumentError(
            "interval has no representative! Run `ripserer` with `reps=true`"
        ))
    elseif !(birth(interval) ≤ _birth_or_value(r) < death(interval))
        return simplex_type(filtration, 1)[]
    else
        reps = filter!(simplex.(representative(interval))) do sx
            birth(sx) ≤ _birth_or_value(r)
        end
        g = OneSkeleton(filtration, r, reps)
        return _find_cycle(g)
    end
end
