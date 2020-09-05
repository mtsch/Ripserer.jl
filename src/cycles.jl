struct OneSkeleton{T, F<:AbstractFiltration, S<:AbstractSimplex{1, T}} <: AbstractGraph{Int}
    filtration::F
    threshold::T
    removed::Set{S}
end
function OneSkeleton(filtration, threshold, removed=())
    S = simplex_type(filtration, 1)
    removed = Set{S}(removed)
    return OneSkeleton(filtration, threshold, removed)
end

_in(σ::S, g::OneSkeleton{S}) where S = !isnothing(σ) && σ ≤ g.threshold && σ ∉ g.removed
_in(σ, g::OneSkeleton) = !isnothing(σ) && birth(σ) ≤ g.threshold && σ ∉ g.removed

LightGraphs.edgetype(::OneSkeleton) = Edge{Int}

function LightGraphs.has_edge(g::OneSkeleton, u::Integer, v::Integer)
    σ = unsafe_simplex(g.filtration, Val(1), TupleTools.sort((v, u), rev=true))
    return _in(σ, g)
end

function LightGraphs.edges(g::OneSkeleton)
    result = edgetype(g)[]
    for sx in Ripserer.edges(g.filtration)
        if birth(sx) ≤ g.threshold && sx ∉ g.repset
            u, v = sx
            push!(result, Edge(u, v))
        end
    end
    return result
end
function LightGraphs.outneighbors(g::OneSkeleton, u::Integer)
    root = simplex(g.filtration, Val(0), (u,))
    neighbors = Int[]
    for sx in Ripserer.coboundary(g.filtration, root)
        if _in(sx, g)
            v, w = sx
            if v == u
                push!(neighbors, w)
            else
                push!(neighbors, v)
            end
        end
    end
    return neighbors
end

Base.eltype(::OneSkeleton) = Int
LightGraphs.has_vertex(g::OneSkeleton, u) = 1 ≤ u ≤ nv(g)
LightGraphs.ne(g::OneSkeleton) = length(edges(g))
LightGraphs.nv(g::OneSkeleton) = length(nv(g.filtration))
LightGraphs.vertices(g::OneSkeleton) = vertices(g)
LightGraphs.is_directed(::OneSkeleton) = false
LightGraphs.is_directed(::Type{<:OneSkeleton}) = false

# need default dist for filtrations...
LightGraphs.weights(g::OneSkeleton) = dist(g.filtration)

function reconstruct_cycle(
    filtration::AbstractFiltration{<:Any, T}, interval, r=birth_simplex(interval)
) where T
    g = OneSkeleton(filtration, interval, r)
    fw = floyd_warshall_shortest_paths(g)
    parents = fw.parents
    dists = fw.dists
    best_u, best_v = 0, 0
    best_sx = first(representative(interval))
    best_dist = typemax(T)
    for sx in Iterators.filter(σ -> birth(σ) ≤ g.threshold, g.repset)
        u, v = linear.(Ref(g), sx)
        if dists[u, v] < best_dist
            best_dist = dists[u, v]
            best_u, best_v = u, v
            best_sx = sx
        end
    end
    if best_dist == typemax(best_dist)
        error("no cycle found!")
    end
    u, v = best_u, best_v

    result = simplex_type(filtration, 1)[]
    while v ≠ u
        w = parents[u, v]
        @assert w ≠ v
        @assert w ≠ 0
        push!(result, simplex(filtration, Val(1), (u, v)))
        v = w
    end
    @assert maximum(birth, result) ≤ g.threshold
    push!(result, best_sx)
    return result
end
