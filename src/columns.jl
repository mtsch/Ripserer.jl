# cofaces ================================================================================ #
"""
    coface(reduction_state, simplex, new_vertex)

Get coface by adding `new_vertex` to `simplex`.
"""
function coface(st::ReductionState{M, T}, simplex::DiameterSimplex{M, T},
                new_vertex, dim) where {M, T}
    vxs = vertices(st, simplex, dim)

    index = 0
    new_vertex_done = false
    k = dim + 2
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
    dim     ::Int
end

Base.IteratorSize(::Type{CoboundaryIterator}) =
    Base.SizeUnknown()
Base.IteratorEltype(::Type{CoboundaryIterator}) =
    Base.HasEltype()
Base.eltype(::Type{CoboundaryIterator{M, T}}) where {M, T} =
    DiameterSimplex{M, T}

function Base.iterate(ni::CoboundaryIterator, i = 1) where {M, T}
    vxs = vertices(ni.state, ni.simplex, ni.dim)
    while i <= n_vertices(ni.state) && !is_connected(ni.state, i, vxs)
        i += 1
    end
    if i > n_vertices(ni.state)
        nothing
    else
        coface(ni.state, ni.simplex, i, ni.dim), i + 1
    end
end

coboundary(st, simplex, dim) =
    CoboundaryIterator(st, simplex, dim)

# columns ================================================================================ #
const Column{M, T} =
    BinaryHeap{DiameterSimplex{M, T}, DiameterSimplexComparer}

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
    if !isnothing(pivot)
        push!(column, pivot)
    end
    pivot
end

# main algo stuff ======================================================================== #
struct ReductionMatrices{M, T, R<:ReductionState{M, T}}
    state             ::R
    reduction_matrix  ::CompressedSparseMatrix{DiameterSimplex{M, T}}
    column_index      ::Dict{Int, Int} # index(sx) => column index of reduction_matrix
    working_column    ::Column{M, T}
    reduction_entries ::Column{M, T}
    current_dim       ::Ref{Int}
end

ReductionMatrices(st::ReductionState{M, T}, dim) where {M, T} =
    ReductionMatrices(st, CompressedSparseMatrix{DiameterSimplex{M, T}}(),
                      Dict{Int, Int}(), Column{M, T}(), Column{M, T}(), Ref(dim))

"""
    add!(rm::ReductionMatrices, index)

Add column with column index `index` multiplied by the correct factor to `working_column`.
Also record the addition in `reduction_matrix`.
"""
# add mora vedet katere simplekse je dodal, ne njihovih coboundaryjev
function add!(rm::ReductionMatrices, index)
    inv_pivot = inv(pivot(rm.working_column))
    for simplex in rm.reduction_matrix[index]
        push!(rm.reduction_entries, -simplex * inv_pivot)
        for coface in coboundary(rm.state, simplex, rm.current_dim[])
            push!(rm.working_column, -coface * inv_pivot)
        end
    end
    pivot(rm.working_column)
end

function reduce_working_column!(rm::ReductionMatrices, res, column_simplex)
    # initialize!(rm, column_simplex)
    while !isempty(rm.working_column)
        pop!(rm.working_column)
    end
    while !isempty(rm.reduction_entries)
        pop!(rm.reduction_entries)
    end

    for coface in coboundary(rm.state, column_simplex, rm.current_dim[])
        push!(rm.working_column, coface)
    end
    # end initialize

    add_column!(rm.reduction_matrix)
    push!(rm.reduction_matrix, column_simplex)

    current_pivot = pivot(rm.working_column)
    while !isnothing(current_pivot) && haskey(rm.column_index, index(current_pivot))
        current_pivot = add!(rm, rm.column_index[index(current_pivot)])
    end
    if isnothing(current_pivot)
        #println("reduced")
    else
        death = diam(current_pivot)
        birth = diam(column_simplex)
        #println("int. ($birth, $death)")
        if death > birth
            push!(res, (birth, death))
        end

        rm.column_index[index(current_pivot)] = length(rm.reduction_matrix)
        current_entry = pop_pivot!(rm.reduction_entries)
        while !isnothing(current_entry)
            # tukaj se pusha simplekse IZ KATERIH SI RAČUNAL COBOUNDARY!!
            push!(rm.reduction_matrix, current_entry)
            current_entry = pop_pivot!(rm.reduction_entries)
        end
    end
end

"""
    compute_0_dim_pairs!(reduction_state, columns)

Compute 0-dimensional persistent homology using Kruskal's Algorithm.
"""
function compute_0_dim_pairs!(st::ReductionState{M, T},
                              columns::AbstractVector{DiameterSimplex{M, T}}) where {M, T}
    dset = IntDisjointSets(n_vertices(st))
    res = Tuple{T, T}[]

    for (l, (u, v)) in edges(st)
        i = find_root(dset, u)
        j = find_root(dset, v)
        if i ≠ j
            union!(dset, i, j)
            if l > 0
                push!(res, (zero(T), T(l)))
            end
        else
            push!(columns, DiameterSimplex{M}(st, T(l), (u, v), 1))
        end
    end
    for _ in 1:num_groups(dset)
        push!(res, (zero(T), typemax(T)))
    end
    reverse!(columns)
    res
end

function compute_pairs!(rm::ReductionMatrices{M, T}, columns_to_reduce, dim) where {M, T}
    res = Tuple{T, T}[]
    rm.current_dim[] = dim
    for column in columns_to_reduce
        reduce_working_column!(rm, res, column)
    end
    res
end

function ripserer(dist, dim_max, modulus) # todo thresholding
    st = ReductionState{modulus}(dist, dim_max)
    push!(res, compute_0_dim_pairs!(st, columns))

    if dim_max > 0
        rm = ReductionMatrices(st, 1)
        for dim in 1:dim_max
            push!(res, compute_pairs!(rm, columns, dim))
            if dim < dim_max
                error("not implemented for dim > 1")
            end
        end
    end
end
