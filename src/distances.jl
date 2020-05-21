# Tried to approximate the approach from
# https://mrzv.org/publications/geometry-helps-distances-persistence-diagrams/alenex/ and
# https://www2.cs.arizona.edu/~alon/papers/match.pdf with NearestNeighbors.jl but the
# allocations were so high, it was slower for moderately sized inputs. This approach is also
# much much simpler.
# TODO: try it again some time.

"""
    Matching

A matching between two persistence diagrams.

# Methods

* [`distance(::Matching)`](@ref)
* [`matching(::Matching)`](@ref)
"""
struct Matching{T, L, R}
    left::L
    right::R
    weight::T
    match_left::Vector{Int}
    match_right::Vector{Int}
end

function Matching(left, right, weight, match)
    n = length(left)
    m = length(right)
    match_left = fill(0, n + m)
    match_right = fill(0, m + n)
    for (l, r) in match
        match_left[l] = r
        match_right[r] = l
    end
    return Matching(left, right, weight, match_left, match_right)
end

"""
    distance(::Matching)

Get the weight of a `Matching` object.
"""
distance(match::Matching) = match.weight

Base.length(match::Matching) = length(matching(match))

"""
    matching(::Matching)

Get the matching of a `Matching` object represented by a vector of pairs of intervals.
"""
function matching(match::Matching{T}) where T
    n = length(match.left)
    m = length(match.right)

    # We convert both sides to `PersistenceInterval{T}` in case their types don't match.
    P = PersistenceInterval{T, Nothing}
    result = Pair{P, P}[]
    for i in 1:n
        l = P(match.left[i]...)
        j = match.match_left[i]
        if j > m
            r = P(birth(l), birth(l))
        else
            r = P(match.right[j]...)
        end
        push!(result, l => r)
    end
    # Collect rights matched to diagonals.
    for j in 1:m
        r = P(match.right[j]...)
        i = match.match_right[j]
        if i > n
            l = P(birth(r), birth(r))
        else
            continue
        end
        push!(result, l => r)
    end
    sort!(result)
end

Base.show(io::IO, match::Matching) =
    print(io, "$(length(match))-element Matching with weight $(match.weight)")

function Base.show(io::IO, ::MIME"text/plain", match::Matching)
    print(io, match)
    if length(match) > 0
        print(io, ":")
        show_intervals(io, matching(match))
    end
end


"""
    adj_matrix(diag1::PersistenceDiagram, diag2::PersistenceDiagram, inf)

Get the adjacency matrix of the matching between `diag1` and `diag2`. Distances between
diagonal points and values that should not be matched with them are set to `inf`.

For `length(diag1) == n` and `length(diag2) == m`, it returns a ``(n m) × (m n)`` matrix.

# Example

```jldoctest
diag1 = PersistenceDiagram(0, [(0.0, 1.0), (3.0, 4.5)])
diag2 = PersistenceDiagram(0, [(0.0, 1.0), (4.0, 5.0), (4.0, 7.0)])

adj_matrix(diag1, diag2, Inf)

# output

5×5 Array{Float64,2}:
  0.0   3.5   1.0  Inf   Inf
  4.0   1.0  Inf    1.0  Inf
  6.0   2.5  Inf   Inf    3.0
  1.0  Inf    0.0   0.0   0.0
 Inf    1.5   0.0   0.0   0.0
```
"""
function adj_matrix(diag1::PersistenceDiagram, diag2::PersistenceDiagram, inf)
    function _to_matrix(diag, T)
        pts = [(T(birth(i)), T(death(i))) for i in diag]
        return reshape(reinterpret(T, pts), (2, length(pts)))
    end
    T = promote_type(dist_type(diag1), dist_type(diag2))

    n = length(diag1)
    m = length(diag2)
    adj = Matrix{promote_type(T, typeof(inf))}(undef, n + m, m + n)
    adj .= inf

    dists = pairwise(Chebyshev(), _to_matrix(diag2, T), _to_matrix(diag1, T), dims=2)
    adj[1:m, 1:n] .= dists
    for i in 1:n
        adj[i + m, i] = persistence(diag1[i])
    end
    for j in 1:m
        adj[j, j + n] = persistence(diag2[j])
    end
    adj[m + 1:m + n, n + 1:n + m] .= zero(T)
    return adj
end

"""
    BottleneckGraph{T}

Representation of the bipartite graph used for computing bottleneck distance via the
Hopcroft-Karp algorithm. In all the following functions, `left` and `right` refer to the
vertex sets of the graph. The graph has `n + m` vertices in each set corresponding to the
numbers of points in the diagrams plus the diagonals.

# Fields

* `adj::Matrix{T}`: the adjacency matrix.
* `match_left::Vector{Int}`: matches of left vertices.
* `match_right::Vector{Int}`: matches of right vertices.
* `edges::Vector{T}`: edge lengths, unique and sorted.
* `n::Int`: number of intervals in left diagram.
* `m::Int`: number of intervals in right diagram.
"""
struct BottleneckGraph{T}
    adj::Matrix{T}

    match_left::Vector{Int}
    match_right::Vector{Int}

    edges::Vector{T}

    n::Int
    m::Int
end

function BottleneckGraph(diag1::PersistenceDiagram, diag2::PersistenceDiagram)
    diag1 = filter(isfinite, diag1)
    diag2 = filter(isfinite, diag2)

    n = length(diag1)
    m = length(diag2)
    T = promote_type(dist_type(diag1), dist_type(diag2))
    adj = adj_matrix(diag1, diag2, typemax(T))

    edges = vcat(vec(adj), persistence.(diag1), persistence.(diag2))
    edges = filter!(e -> e < typemax(T), edges) |> unique! |> sort!

    return BottleneckGraph(adj, fill(0, n + m), fill(0, m + n), edges, n, m)
end

function left_neighbors!(buff, graph::BottleneckGraph, vertices, ε, pred)
    empty!(buff)
    for l in vertices
        for r in axes(graph.adj, 1)
            graph.adj[r, l] ≤ ε && pred(r) && push!(buff, r)
        end
    end
    return unique!(buff)
end

function right_neighbors!(buff, graph::BottleneckGraph, vertices)
    empty!(buff)
    for r in vertices
        push!(buff, graph.match_right[r])
    end
    return buff
end

is_exposed_left(graph::BottleneckGraph, l) = graph.match_left[l] == 0

is_exposed_right(graph::BottleneckGraph, r) = graph.match_right[r] == 0

exposed_left(graph::BottleneckGraph) =
    findall(iszero, graph.match_left)

"""
    layer_graph(graph::BottleneckGraph, ε)

Split `graph` into layers by how deep they are from a bfs starting at exposed left
vertices in `graph` only taking into account edges of length smaller than or equal to `ε`.
Return depts of right vertices and maximum depth reached.
"""
function layer_graph(graph::BottleneckGraph, ε)
    depths = fill(0, graph.m + graph.n)
    visited = falses(graph.m + graph.n)
    lefts = exposed_left(graph)
    rights = Int[]
    i = 1
    while true
        left_neighbors!(rights, graph, lefts, ε, r -> !visited[r])
        visited[rights] .= true
        depths[rights] .= i
        if isempty(rights)
            # no augmenting path exists
            return nothing, nothing
        elseif any(r -> is_exposed_right(graph, r), rights)
            return depths, i
        else
            right_neighbors!(lefts, graph, rights)
        end
        i += 1
    end
end

"""
    augmenting_paths(graph::BottleneckGraph, ε)

find a maximal set of augmenting paths in graph, taking only edges with weight less than or
equal to `ε` into account.
"""
function augmenting_paths(graph::BottleneckGraph, ε)
    depths, max_depth = layer_graph(graph, ε)
    paths = Vector{Int}[]
    isnothing(depths) && return paths

    prev_left = fill(0, graph.n + graph.m)
    prev_right = fill(0, graph.m + graph.n)

    rights = Int[]
    lefts = Int[]
    stack = Tuple{Int, Int}[]

    for l_start in exposed_left(graph)
        empty!(stack)
        push!(stack, (l_start, 1))

        while !isempty(stack)
            l, i = pop!(stack)
            left_neighbors!(rights, graph, l, ε, r -> depths[r] == i)
            if i < max_depth
                prev_right[rights] .= l
                right_neighbors!(lefts, graph, rights)
                append!(stack, (l, i + 1) for l in lefts)
            else
                found_path = false
                for r in rights
                    if is_exposed_right(graph, r)
                        prev_right[r] = l
                        path = Int[r]
                        depths[r] = 0
                        while (l = prev_right[r]) ≠ l_start
                            i -= 1
                            r = graph.match_left[l]
                            depths[r] = 0
                            append!(path, (l, r))
                        end
                        push!(path, l_start)
                        reverse!(path)
                        push!(paths, path)
                        found_path = true
                        break
                    end
                end
                found_path && break
            end
        end
    end

    return paths
end

function match!(graph::BottleneckGraph, l, r)
    graph.match_left[l] = r
    graph.match_right[r] = l
    return (l, r)
end

function matched(graph::BottleneckGraph, l, r)
    return graph.match_left[l] == r && graph.match_right[r] == l
end

function unmatch!(graph::BottleneckGraph, l, r)
    @assert matched(graph, l, r)
    graph.match_left[l] = 0
    graph.match_right[r] = 0
    return (0, 0)
end

function unmatch_all!(graph::BottleneckGraph)
    graph.match_left .= 0
    graph.match_right .= 0
    return (0, 0)
end

function hopcroft_karp!(graph, ε)
    unmatch_all!(graph)
    paths = augmenting_paths(graph, ε)
    while !isempty(paths)
        for p in paths
            for i in 1:length(p)-1
                if matched(graph, p[i], p[i + 1])
                    unmatch!(graph, p[i], p[i+1])
                else
                    match!(graph, p[i], p[i+1])
                end
            end
        end
        paths = augmenting_paths(graph, ε)
    end
    matching = [(i, graph.match_left[i])
                for i in 1:graph.n + graph.m if graph.match_left[i] ≠ 0]
    is_maximum = length(matching) == graph.n + graph.m
    return matching, is_maximum
end

"""
    Bottleneck

Use this object to find the bottleneck distance or matching between persistence diagrams.
The distance value is equal to

```math
W_\\infty(X, Y) = \\inf_{\\eta:X\\rightarrow Y} \\sup_{x\\in X} ||x-\\eta(x)||_\\infty,
```

where ``X`` and ``Y`` are the persistence diagrams and ``\\eta`` is a perfect matching
between the intervals. Note the ``X`` and ``Y`` don't need to have the same number of
points, as the diagonal points are considered in the matching as well.

# Methods

* [`matching(::Bottleneck, ::Any, ::Any)`](@ref): construct a bottleneck [`Matching`](@ref).
* [`distance(::Bottleneck, ::Any, ::Any)`](@ref): find the bottleneck distance.
"""
struct Bottleneck end

"""
    matching(::Bottleneck, diag1, diag2)

Find the bottleneck matching between persistence diagrams `diag1` and `diag2`. Infinite
intervals are ignored.

```jldoctest
diag1 = PersistenceDiagram(0, [(1.0, 2.0), (5.0, 8.0)])
diag2 = PersistenceDiagram(0, [(1.0, 2.0), (3.0, 4.0), (5.0, 10.0)])
matching(Bottleneck(), diag1, diag2)

# Example

# output

3-element Matching with weight 2.0:
 [1.0, 2.0) => [1.0, 2.0)
 [3.0, 3.0) => [3.0, 4.0)
 [5.0, 8.0) => [5.0, 10.0)
```

# See also

* [`Bottleneck`](@ref)
* [`distance`](@ref)
"""
function matching(::Bottleneck, diag1, diag2)
    diag1 = filter(isfinite, diag1)
    diag2 = filter(isfinite, diag2)

    graph = BottleneckGraph(diag1, diag2)
    edges = graph.edges

    lo = 1
    hi = length(edges)
    while lo < hi - 1
        m = lo + ((hi - lo) >>> 0x01)
        _, succ = hopcroft_karp!(graph, edges[m])
        if succ
            hi = m
        else
            lo = m
        end
    end
    match, succ = hopcroft_karp!(graph, edges[lo])
    distance = edges[lo]
    if !succ
        distance = edges[hi]
        match, _ = hopcroft_karp!(graph, edges[hi])
    end
    @assert length(match) == length(diag1) + length(diag2)

    return Matching(diag1, diag2, distance, match)
end

"""
    distance(::Bottleneck, diag1, diag2)

Compute the bottleneck distance between persistence diagrams `diag1` and `diag2`. Infinite
intervals are ignored.

# Example

```jldoctest
diag1 = PersistenceDiagram(0, [(1.0, 2.0), (5.0, 8.0)])
diag2 = PersistenceDiagram(0, [(1.0, 2.0), (3.0, 4.0), (5.0, 10.0)])
distance(Bottleneck(), diag1, diag2)

# output

2.0
```

# See also

* [`Bottleneck`](@ref)
* [`matching`](@ref)
"""
distance(::Bottleneck, diag1::PersistenceDiagram, diag2::PersistenceDiagram) =
    matching(Bottleneck(), diag1, diag2).weight

"""
    Wasserstein(q=1)

Use this object to find the Wasserstein distance or matching between persistence diagrams.
The distance value is equal to

```math
W_q(X,Y)=\\left[\\inf_{\\eta:X\\rightarrow Y}\\sum_{x\\in X}||x-\\eta(x)||_\\infty^q\\right],
```

where ``X`` and ``Y`` are the persistence diagrams and ``\\eta`` is a perfect matching
between the intervals. Note the ``X`` and ``Y`` don't need to have the same number of
points, as the diagonal points are considered in the matching as well.

# Methods

* [`matching(::Wasserstein, ::Any, ::Any)`](@ref): construct a bottleneck [`Matching`](@ref).
* [`distance(::Wasserstein, ::Any, ::Any)`](@ref): find the bottleneck distance.
"""
struct Wasserstein{T}
    q::T

    Wasserstein(q=1) = new{typeof(q)}(q)
end

"""
    matching(::Wasserstein, diag1, diag2)

Find the Wasserstein matching between persistence diagrams `diag1` and `diag2`. Infinite
intervals are ignored.

# Example

```jldoctest
diag1 = PersistenceDiagram(0, [(1.0, 2.0), (5.0, 8.0)])
diag2 = PersistenceDiagram(0, [(1.0, 2.0), (3.0, 4.0), (5.0, 10.0)])
matching(Wasserstein(), diag1, diag2)

# output

3-element Matching with weight 3.0:
 [1.0, 2.0) => [1.0, 2.0)
 [3.0, 3.0) => [3.0, 4.0)
 [5.0, 8.0) => [5.0, 10.0)
```

# See also

* [`Wasserstein`](@ref)
* [`distance`](@ref)
"""
function matching(w::Wasserstein, diag1::PersistenceDiagram, diag2::PersistenceDiagram)
    adj = adj_matrix(diag2, diag1, missing).^w.q
    match = collect(enumerate(hungarian(adj)[1]))
    distance = sum(adj[i, j] for (i, j) in match)^(w.q == 1 ? 1 : 1/w.q)

    return Matching(diag1, diag2, distance, match)
end

"""
    distance(::Wasserstein, diag1, diag2)

Compute the Wasserstein distance between persistence diagrams `diag1` and `diag2`. Infinite
intervals are ignored.

# Example

```jldoctest
diag1 = PersistenceDiagram(0, [(1.0, 2.0), (5.0, 8.0)])
diag2 = PersistenceDiagram(0, [(1.0, 2.0), (3.0, 4.0), (5.0, 10.0)])
distance(Wasserstein(), diag1, diag2)

# output

3.0
```

# See also

* [`Wasserstein`](@ref)
* [`matching`](@ref)
"""
distance(w::Wasserstein, diag1::PersistenceDiagram, diag2::PersistenceDiagram) =
    matching(w, diag1, diag2).weight
