
const DiameterSimplexHeap{M} =
    BinaryHeap{DiameterSimplex{M}, DiameterSimplexComparer}

struct CurrentColumn{M}
    heap            ::DiameterSimplexHeap{M}
    vertex_buffer   ::Vector{Int}
    #neighbor_buffer ::Set{Int}
end
CurrentColumn{M}() where M =
    CurrentColumn(DiameterSimplexHeap{M}(), Int[])

function diam(vertices, dist)
    maximum(dist[i, j] for i in vertices, j in vertices)
end

"""
    coface(simplex::DiameterSimplex{M}, vertices, new_vertex, dist, binomial)

Creat coface by adding `new_vertex` to `simplex` with `vertices`. `dist` is used to
determine the new diameter.
"""
function coface(simplex::DiameterSimplex{M}, vertices, new_vertex, dist, binomial) where M
    index = 0
    new_vertex_done = false
    k = length(vertices) + 1
    i = 1
    diameter = diam(simplex)
    coefficient = coef(simplex)
    while k > 0
        if !new_vertex_done && (k == 1 || new_vertex > vertices[i])
            new_vertex_done = true
            index += binomial(new_vertex - 1, k)
            diameter = max(diameter, maximum(dist[new_vertex, j] for j in vertices))
            coefficient = (k % 1 == 1 ? M - 1 : 1) * coefficient % M
            k -= 1
        else
            index += binomial(vertices[i] - 1, k)
            i += 1
            k -= 1
        end
    end
    DiameterSimplex{M}(diameter, index + 1, coefficient)
end

"""
    get_common_neighbors(vertices, dist::AbstractMatrix)

Get neighbors `vertices` have in common. If `dist` is sparse, it is treated as a graph,
otherwise, all indices except for `vertices` are returned.
"""
function get_common_neighbors(vertices, dist::SparseMatrixCSC)
    colptr = dist.colptr
    rowval = dist.rowval

    common = Set(rowval[colptr[vertices[1]]:colptr[vertices[1]+1]-1])
    for v in vertices
        intersect!(common, Set(rowval[colptr[v]:colptr[v+1]-1]))
    end
    common
end

function get_common_neighbors(vertices, dist::Matrix)
    setdiff!(BitSet(1:size(dist, 1)), vertices)
end

"""
    initialize!(column, simplex, dim, dist, binomial)

Initialize column by putting all cofaces of `dim`-dimensional `simplex` on `column`'s heap.

# TODO emergent pairs.
"""
function initialize!(col::CurrentColumn{M}, simplex, dim, dist, binomial) where M
    heap = col.heap
    vertices = get_vertices!(col.vertex_buffer, simplex, dim, size(dist, 1), binomial)
    common = get_common_neighbors(vertices, dist)
    for v in common
        push!(heap, coface(simplex, vertices, v, dist, binomial))
    end
    col
end

function pop_pivot!(col::CurrentColumn{M}) where M
    heap = col.heap
    isempty(heap) && return nothing

    pivot = pop!(heap)
    while !isempty(heap)
        if coef(pivot) == 0
            pivot = pop!(heap)
        elseif index(top(heap)) == index(pivot)
            pivot += pop!(heap)
        else
            break
        end
    end
    coef(pivot) == 0 ? nothing : pivot
end
