# cofaces ================================================================================ #
"""
    coface(reduction_state, simplex, new_vertex)

Get coface by adding `new_vertex` to `simplex`.
"""
function coface(st::ReductionState{M, T}, simplex::DiameterSimplex{M, T},
                new_vertex) where {M, T}
    vxs = vertices(st, simplex)

    index = 0
    new_vertex_done = false
    k = dim(st) + 2
    i = 1
    diameter = diam(simplex)
    coefficient = coef(simplex)
    while k > 0
        if !new_vertex_done && (k == 1 || new_vertex > vxs[i])
            new_vertex_done = true
            index += binomial(st, new_vertex - 1, k)
            diameter = max(diameter, maximum(dist(st, new_vertex, j) for j in vxs))
            coefficient = (k % 1 == 1 ? M - 1 : 1) * coefficient % M
            k -= 1
        else
            index += binomial(st, vxs[i] - 1, k)
            i += 1
            k -= 1
        end
    end
    DiameterSimplex{M}(diameter, index + 1, coefficient)
end

struct CoboundaryIterator{M, T, R<:ReductionState{M, T}}
    state   ::R
    simplex ::DiameterSimplex{M, T}
end

Base.IteratorSize(::Type{CoboundaryIterator}) =
    Base.SizeUnknown()
Base.IteratorEltype(::Type{CoboundaryIterator}) =
    Base.HasEltype()
Base.eltype(::Type{CoboundaryIterator{M, T}}) where {M, T} =
    DiameterSimplex{M, T}

function Base.iterate(ni::CoboundaryIterator, i = 1) where {M, T}
    vxs = vertices(ni.state, ni.simplex)
    while i <= n_vertices(ni.state) && !is_connected(ni.state, i, vxs)
        i += 1
    end
    if i > n_vertices(ni.state)
        nothing
    else
        coface(ni.state, ni.simplex, i), i + 1
    end
end

coboundary(st, simplex) =
    CoboundaryIterator(st, simplex)

# columns ================================================================================ #
const Column{M, T} =
    BinaryHeap{DiameterSimplex{M, T}, DiameterSimplexComparer}

"""
    initialize!(column, vertex_buffer, simplex, dim, dist, binomial)

Initialize column by putting all cofaces of `dim`-dimensional `simplex` on `column` heap.

# TODO emergent pairs.
"""
function initialize!(column::Column, #=vertex_buffer,=# simplex, dim, dist, binomial)
    vertices = get_vertices!(#=vertex_buffer=#Int[], simplex, dim, size(dist, 1), binomial)
    common = get_common_neighbors(vertices, dist)
    for v in common
        push!(column, coface(simplex, vertices, v, dist, binomial))
    end
    column
end

function pop_pivot!(column::Column)
    isempty(column) && return nothing

    pivot = pop!(column)
    while !isempty(column)
        if coef(pivot) == 0
            pivot = pop!(column)
        elseif index(top(column)) == index(pivot)
            pivot += pop!(column)
        else
            break
        end
    end
    coef(pivot) == 0 ? nothing : pivot
end

function pivot(column::Column)
    pivot = pop_pivot!(column)
    push!(column, pivot)
    pivot
end

struct ReductionMatrices{M, T}
    reduction_matrix ::CompressedSparseMatrix{DiameterSimplex{M, T}}
    pivot_index      ::Dict{T, Simplex{M, T}}
    working_column   ::Column{M, T}
end

function reconstruct!(pivot, reduction_matrix, index, dim, dist, binomial)
    # greš po indeksih v matriki, dodaš vse koboundarije, sproti addaš
    # cofacets of simplex
    for simplex in reduction_matrix[index]
    end
end

# za add:
# * naštimaš drug stolpec
# * vzameš vn pivot od obeh -> se izničta
# * ponavljaš dokler se izničujeta
# * fukneš vse iz col2 na col1

# boljš:
# * greš po matriki in mečeš vse na col1
# * popneš pivot
# * če še obstaja je treba več dodajat
function add!(col1, col2, reduction_matrix, index, dim, dist, binomial)
    pivot = pop_pivot!(col)
    for other in reduction_matrix[index]
        initialize!(col2, other, dim, dist, binomial)
    end
end
