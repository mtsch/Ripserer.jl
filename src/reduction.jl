"""
    compute_0_dim_pairs!(columns, dist, binomial)

Compute 0-dimensional persistent homology using Kruskal's Algorithm.
"""
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
