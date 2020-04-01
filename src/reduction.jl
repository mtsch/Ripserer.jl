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

#assemble_columns_to_reduce

function compute_pairs!(columns::AbstractVector{DiameterSimplex{M, T}},
                        pivot_dict, dim, dist, binomial) where {M, T}
    reduction_matrix = CompressedSparseMatrix{DiameterSimplex{M, T}}()

    for i in eachindex(columns)
        column = columns[i]
        push!(reduction_matrix, column)

        working_reduction_column = CurrentColumn()
        working_coboundary = CurrentColumn()

        initialize!(working_coboundary, column, dim, dist)
        pivot = pop_pivot!(working_coboundary)

        while true
            if isnothing(pivot)
                # pogledaš v dict
                # dobiš drug pivot (key)
                # indeks stolpca (value)
                if true # je v mapu
                # factor = -pivot / other_pivot
                # add!(reduction_matrix, columns, indeks stolpca, factor, dim, oba stolpca)
                # pivot = get_pivot(coboundary)
                else
                # pivot_dict[simplex(pivot)] = i
                # while true
                #   e = pop_pivot!(working reduction)
                #   isnothing(e) && break
                #   @assert coef(e) > 0
                #   reduction_matrix.push_back(e)
                end
            else
                break
            end
        end
    end

        # get column and its diameter
        # append column to reduction matrix
        # initialize working column

        # while true
        #  if get_index(pivot == -1), break
        #  pair = find pivot in pivot_column_index
        #  if pair ≠ pivot_column_index[end]
        #    sx = get simplex from pivot
        #    index_col_to_reduce = index of hash map
        #    factor = modulus - coef(pivot) * inv(coef(other_pivot)) % M
        #    add_coboundary(reduction_matrix, columns, index_col_to_reduce, factor, dim, working_column, working_coboundary)
        #    pivot = get_pivot(working_coboundary)
        #  else
        #    pivot_column_index.insert(simplex(pivot), index_column_to_reduce
        #    while true
        #      e = pop_pivot
        #      if index(e) == -1 break
        #      reduction_matrix.push_back(e)
        #   end
        # end
end
