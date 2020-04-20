# binomial table ========================================================================= #
"""
    Binomial(n_max, k_max)

Table of precomputed binomial coefficients up to `n_max` and `k_max`. Can be called like a
function and should be identical to [`Base.binomial`](@ref) for values of `0 ≤ n ≤ n_max`
and `0 ≤ k ≤ k_max`.
"""
struct Binomial
    table::Matrix{Int64}
end

function Binomial(n, k)
    table = zeros(Int, n+1, k+1)
    for i in 1:n+1
        table[i, 1] = 1;
        for j in 2:min(i, k+1)
            table[i, j] = table[i-1, j-1] + table[i-1, j];
            if (i <= k)
                table[i, i] = 1
            end
        end
    end
    Binomial(table)
end

Base.show(io::IO, bin::Binomial) =
    print(io, "Binomial$(size(bin.table) .- 1)")
@propagate_inbounds (bin::Binomial)(n, k) =
    bin.table[n+1, k+1]

# coboundary enumerator ================================================================== #
"""
    Coboundary{S<:AbstractSimplex, F<:AbstractFiltration}

As the name suggests, this type is used to enumerate coboundaries of simplices via the
`coboundary` function.

# Constructor

    Coboundary(::AbstractFiltration, dim_max::Integer)

# Usage

TODO
"""
struct Coboundary{S, F<:AbstractFiltration{<:Any, S}}
    filtration   ::F
    binomial     ::Binomial
    dim_max      ::Int
    vertex_cache ::Vector{Int}
end

Coboundary(filtration::AbstractFiltration, dim_max::Integer) =
    Coboundary(
        filtration, Binomial(n_vertices(filtration), dim_max+2), dim_max, Int[])

Base.binomial(cob::Coboundary, n, k) =
    cob.binomial(n, k)
dim_max(cob::Coboundary) =
    cob.dim_max
diam(cob::Coboundary, args...) =
    diam(cob.filtration, args...)
n_vertices(cob::Coboundary) =
    n_vertices(cob.filtration)

"""
    index(cob::Coboundary, vertices)

Compute the index from a collection of `vertices`. Vertices must be in descending order.
"""
@propagate_inbounds function index(cob::Coboundary, vertices)
    res = 0
    for l in eachindex(vertices)
        res += binomial(cob, vertices[end - l + 1] - 1, l)
    end
    res + 1
end

"""
    find_max_vertex(filtration, idx, k)

Find largest vertex index of vertex for which `binomial(i, k) ≤ idx` holds. This finds the
index of the first vertex in simplex.
"""
@propagate_inbounds function find_max_vertex(cob::Coboundary, idx, k, n_max)
    hi = n_max + 1
    lo = k - 1
    while lo < hi - 1
        m = lo + ((hi - lo) >>> 0x01)
        if binomial(cob, m, k) ≤ idx
            lo = m
        else
            hi = m
        end
    end
    lo
end

"""
    get_vertices!(filtration, index)

Find and copy vertices of simplex with `index` to `filtration`'s vertex cache.
"""
function get_vertices!(cob::Coboundary, index, dim)
    resize!(cob.vertex_cache, dim + 1)
    index = index - 1
    v = n_vertices(cob.filtration)
    @inbounds for (i, k) in enumerate(dim+1:-1:1)
        v = find_max_vertex(cob, index, k, v-1)

        cob.vertex_cache[i] = v + 1
        index -= binomial(cob, v, k)
    end
    cob.vertex_cache
end

"""
    vertices(coboundary_enumerator, simplex, dim)

    vertices(coboundary_enumerator, index, dim)

Get vertices of `simplex`. Vertices are only recomputed when the vertex cache in
`coboundary_enumerator` is invalid.
"""
vertices(cob::Coboundary{S}, sx::S, dim) where S =
    vertices(cob, index(sx), dim)
function vertices(cob::Coboundary, idx::Integer, dim)
    # If we already have the correct vertices there is no need to recalculate them.
    # Checking if we do is orders of magnitude faster than getting the vertices again.
    if length(cob.vertex_cache) != dim+1 || index(cob, cob.vertex_cache) != idx
        get_vertices!(cob, idx, dim)
    end
    cob.vertex_cache
end

# coboundary ============================================================================= #
struct CoboundaryIterator{A, S<:AbstractSimplex, C<:Coboundary{S}}
    coboundary ::C
    simplex    ::S
    dim        ::Int

    function CoboundaryIterator{A}(cob, sx, dim) where A
        dim ≤ dim_max(cob) + 1 ||
            throw(ArgumentError("`dim` should be smaller or equal than `dim_max + 1`"))
        new{A, typeof(sx), typeof(cob)}(cob, sx, dim)
    end
end

Base.IteratorSize(::Type{CoboundaryIterator}) =
    Base.SizeUnknown()
Base.IteratorEltype(::Type{CoboundaryIterator}) =
    Base.HasEltype()
Base.eltype(::Type{CoboundaryIterator{<:Any, S}}) where S =
    S

"""
    coboundary(filtration, simplex, dim, [all_cofaces])

Return an iterator that iterates over all cofaces of `simplex` of dimension `dim + 1` in
decreasing order by index.

If `all_cofaces` is `::Val{false}`, only return cofaces where the added vertex has an index
higher than all other vertices. This is used in sparse `assemble_columns!` to find all
simplices.
"""
(cob::Coboundary)(simplex::AbstractSimplex, dim, ::Val{false}) =
    CoboundaryIterator{false}(cob, simplex, dim)
(cob::Coboundary)(simplex::AbstractSimplex, dim) =
    CoboundaryIterator{true}(cob, simplex, dim)

function Base.iterate(
    ci::CoboundaryIterator{all_cofaces, S},
    st = (n_vertices(ci.coboundary), ci.dim + 1, index(ci.simplex) - 1, 0)
) where {all_cofaces, S}
    # bounds are checked in constructor
    diameter = ∞
    v, k, idx_below, idx_above = st

    @inbounds while diameter == ∞
        v -= 1
        if !all_cofaces && binomial(ci.coboundary, v, k) ≤ idx_below
            break
        end
        while v > 0 && v ≥ k && binomial(ci.coboundary, v, k) ≤ idx_below
            idx_below -= binomial(ci.coboundary, v, k)
            idx_above += binomial(ci.coboundary, v, k + 1)
            v -= 1; k -= 1
        end
        if v < k
            break
        else
            # We could avoid this call, but that would make nested coboundary iteration
            # unsafe.
            vxs = vertices(ci.coboundary, ci.simplex, ci.dim)
            diameter = diam(ci.coboundary, ci.simplex, vxs, v+1)
        end
    end
    if diameter != ∞
        coefficient = ifelse(k % 2 == 1, -coef(ci.simplex), coef(ci.simplex))
        new_index = idx_above + binomial(ci.coboundary, v, k + 1) + idx_below + 1

        S(diameter, new_index, coefficient), (v, k, idx_below, idx_above)
    else
        nothing
    end
end
