# since this is the first step, maybe dont pass edges and n and just pass distance matrix.
# maybe dont convert simplex -> vertices, but rather do it the other way around.
# docs
"""
    compute_0_dim_pairs!(columns, edges, n, binomial)

Compute 0-dimensional persistent homology using Kruskal's Algorithm.
"""
function compute_0_dim_pairs!(columns::AbstractVector{DiameterSimplex{M, T}},
                              edges, n, binomial) where {M, T}
    dset = IntDisjointSets(n)
    buff = [0, 0]
    res = Tuple{T, T}[]

    for e in edges
        u, v = get_vertices!(buff, e, 1, n, binomial)
        u_i = find_root(dset, u)
        v_i = find_root(dset, v)
        if u_i ≠ v_i
            union!(dset, u_i, v_i)
            if diam(e) > 0
                push!(res, (0, diam(e)))
            end
        else
            push!(columns, e)
        end
    end
    for _ in 1:num_groups(dset)
        push!(res, (0, typemax(T)))
    end
    reverse!(columns)
    res
end

function compute_0_dim_pairs!(columns::AbstractVector{DiameterSimplex{M, T}},
                              dist, binomial) where {M, T}
    dset = IntDisjointSets(size(dist, 1))
    res = Tuple{T, T}[]

    for (l, (u, v)) in edges(dist)
        i = find_root(dset, u)
        j = find_root(dset, v)
        if i ≠ j
            union!(dset, i, j)
            if l > 0
                push!(res, (zero(T), T(l)))
            end
        else
            push!(columns, DiameterSimplex{M}(T(l), index((u, v), binomial), 1))
        end
    end
    for _ in 1:num_groups(dset)
        push!(res, (zero(T), typemax(T)))
    end
    reverse!(columns)
    res
end

function compute_pairs!(columns::AbstractVector{DiameterSimplex{M, T}},
                        pivot_column_index, dim, binomial) where {M, T}
end
