struct CircumSimplex{D, T, I<:Integer} <: AbstractSimplex{D, T, I}
    index::I
    birth::T
    circum::T

    function CircumSimplex{D, T, I}(index::Integer, birth, circum) where {D, T, I}
        D ≥ 0 || throw(DomainError(D, "dimension must be a non-negative integer"))
        return new{D, T, I}(I(index), T(birth), T(circum))
    end
end

function CircumSimplex{D}(vertices, birth::T, circ::T) where {D, T}
    idx = index(TupleTools.sort(Tuple(vertices), rev=true))
    return CircumSimplex{D, T, eltype(vertices)}(idx, birth, circ)
end
function CircumSimplex{1}(index_or_vertices, birth)
    return CircumSimplex{1}(index_or_vertices, birth, birth)
end

birth(sx::CircumSimplex) = sx.birth
index(sx::CircumSimplex) = abs(sx.index)
Base.sign(sx::CircumSimplex) = sign(sx.index)
function Base.:-(sx::CircumSimplex{D, T, I}) where {D, T, I}
    return CircumSimplex{D, T, I}(-sx.index, birth(sx), circum(sx))
end
circum(sx::CircumSimplex) = sx.circum

function Base.isless(sx1::S, sx2::S) where S<:CircumSimplex
    (birth(sx1), circum(sx1), -index(sx1)) < (birth(sx2), circum(sx2), -index(sx2))
end

# Dense version.
@inline @propagate_inbounds function unsafe_cofacet(
    ::Type{S},
    rips::AbstractRipsFiltration,
    simplex::CircumSimplex,
    cofacet_vertices,
    new_vertex,
    sign,
) where {I, T, D, S<:CircumSimplex{D, T, I}}
    diam = birth(simplex)
    circumference = circum(simplex)
    for v in cofacet_vertices
        v == new_vertex && continue
        d = dist(rips, new_vertex, v)
        if ismissing(d) || d > threshold(rips)
            return nothing
        else
            _d::T = d
            diam = ifelse(_d > diam, _d, diam)
            circumference += _d
        end
    end
    return S(I(sign) * index(cofacet_vertices), diam, circumference)
end

# Sparse version.
@inline @propagate_inbounds function unsafe_cofacet(
    ::Type{S},
    rips::AbstractRipsFiltration,
    simplex::CircumSimplex,
    cofacet_vertices,
    ::Any,
    sign,
    new_edges,
) where {I, T, D, S<:CircumSimplex{D, T, I}}
    diameter = birth(simplex)
    circumference = circum(simplex)
    for e in new_edges
        e > threshold(rips) && return nothing
        diameter = ifelse(e > diameter, e, diameter)
        circumference += e
    end
    return S(I(sign) * index(cofacet_vertices), diameter, circumference)
end

function Ripserer.unsafe_simplex(
    ::Type{S}, rips::AbstractRipsFiltration{I, T}, vertices, sign
) where {I, T, S<:CircumSimplex{<:Any, T, I}}
    if dim(S) == 0
        S(I(sign * vertices[1]), birth(rips, vertices[1]), birth(rips, vertices[1]))
    else
        circumference = zero(T)
        diameter = typemin(T)
        n = length(vertices)
        for i in 1:n, j in i+1:n
            d = dist(rips, vertices[i], vertices[j])
            if ismissing(d) || d > threshold(rips)
                return nothing
            else
                _d::T = d
                circumference += _d
                diameter = max(diameter, _d)
            end
        end
        return S(I(sign) * index(vertices), diameter, circumference)
    end
end

# TODO: remove/replace with FiltrationWrapper
struct EqRips{I, T, R<:AbstractRipsFiltration{I, T}} <: AbstractRipsFiltration{I, T}
    filtration::R
end

dist(e::EqRips) = dist(e.filtration)
dist(e::EqRips, u, v) = dist(e.filtration, u, v)
threshold(e::EqRips) = Ripserer.threshold(e.filtration)
adjacency_matrix(e::EqRips) = adjacency_matrix(e.filtration)
simplex_type(::Type{<:EqRips{I, T}}, D) where {I, T} = CircumSimplex{D, T, I}

emergent_pairs(::EqRips) = false

struct GeodesicRips{I, T, G, A<:AbstractMatrix{T}} <: AbstractRipsFiltration{I, T}
    graph::G
    dist::A
    paths::Matrix{Int}
    threshold::T
end

function GeodesicRips{I}(graph, distmx=weights(graph); threshold=nothing) where I
    fw = floyd_warshall_shortest_paths(graph, distmx)
    if isnothing(threshold)
        threshold = radius(fw.dists)
    end
    T = eltype(fw.dists)
    return GeodesicRips{I, T, typeof(graph), typeof(fw.dists)}(
        graph, fw.dists, fw.parents, T(threshold)
    )
end
GeodesicRips(args...; kwargs...) = GeodesicRips{Int}(args...; kwargs...)

dist(gr::GeodesicRips) = gr.dist
threshold(gr::GeodesicRips) = gr.threshold
simplex_type(::Type{<:GeodesicRips{I, T}}, D) where {I, T} = CircumSimplex{D, T, I}

emergent_pairs(::GeodesicRips) = false

function _edges(u, vs, dist, r)
    res = map(vs) do v
        dist[u, v]
    end
    if maximum(res) ≤ r && minimum(res) > 0
        return res
    else
        return nothing
    end
end

function _visit!(stack, vs, visited)
    for v in vs
        if !visited[v]
            push!(stack, v)
            visited[v] = true
        end
    end
end

function _signed_insert(vertices, vertex)
    n = length(vertices)
    sign = iseven(n) ? 1 : -1
    for i in 0:n-1
        if vertices[i+1] < vertex
            return TupleTools.insertafter(vertices, i, (vertex,)), sign
        end
        sign *= -1
    end
    return TupleTools.insertafter(vertices, n, (vertex,)), sign
end

function _min_cofacet(grips::GeodesicRips, σ)
    diam = birth(σ)
    vxs = vertices(σ)
    starting_vertex = vxs[1]
    τ = nothing

    # DFS from one of the vertices, only taking vertices that are less than birth(σ)
    # away from all the others.
    visited = falses(nv(grips))
    stack = Int[]
    _visit!(stack, neighbors(grips.graph, starting_vertex), visited)

    while !isempty(stack)
        v = pop!(stack)
        if v in vxs
            _visit!(stack, neighbors(grips.graph, v), visited)
        else
            weights = _edges(v, vxs, grips.dist, diam)
            if !isnothing(weights)
                _visit!(stack, neighbors(grips.graph, v), visited)
                #new_vxs = TupleTools.sort(TupleTools.insertafter(Tuple(vxs), 0, (v,)), rev=true)
                new_vxs, sign = _signed_insert(Tuple(vxs), v)
                # Ignore the sign from now, we will get it from the call to boundary.
                candidate = unsafe_cofacet(grips, σ, new_vxs, v, sign, weights)
                if isnothing(τ) || !isnothing(candidate) && τ > candidate
                    τ = candidate
                end
            end
        end
    end
    return τ
end

function find_apparent_pairs(grips::GeodesicRips{<:Any, T}, cols, progress) where T
    # TODO only works fine for dense grips
    S = eltype(cols)
    C = simplex_type(grips, dim(S) + 1)
    cols_left = S[]
    apparent = Tuple{S, C}[]

    if progress
        progbar = Progress(length(cols); desc="Finding apparent pairs... ")
    end
    for σ in cols
        τ = _min_cofacet(grips, σ)
        if isnothing(τ)
            push!(cols_left, σ)
        else
            if σ == maximum(boundary(grips, τ))
                push!(apparent, (σ, τ))
            else
                push!(cols_left, σ)
            end
        end
        progress && next!(progbar)
    end
    progress && printstyled(stderr, "$(length(apparent)) apparent pairs.\n", color=:green)
    return cols_left, apparent
end

function filling(grips, element::AbstractChainElement)
    return filling(grips, simplex(element))
end

function filling(grips, sx::AbstractSimplex)
    if dim(sx) > 2
        throw(ArgumentError("currently only dims up to 2 are supported"))
    end
    vxs = vertices(sx)
    result = simplex_type(grips, 1)[]
    for (u, v) in IterTools.subsets(vxs, Val(2))
        while u ≠ v
            u′ = grips.paths[v, u]
            push!(result, simplex(grips, Val(1), (u, u′)))
            u = u′
        end
    end
    return result
end
