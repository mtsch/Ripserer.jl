const Column{M, T} =
    BinaryHeap{DiameterSimplex{M, T}, DiameterSimplexComparer}

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
    initialize!(column, vertex_buffer, simplex, dim, dist, binomial)

Initialize column by putting all cofaces of `dim`-dimensional `simplex` on `column` heap.

# TODO emergent pairs.
# TODO do we even need the vertex buffer?
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
