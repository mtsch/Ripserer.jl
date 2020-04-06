# cofaces ================================================================================ #
"""
    coface(reduction_state, simplex, new_vertex)

Get coface by adding `new_vertex` to `simplex`.
"""
function coface(st::ReductionState{M, T}, simplex::Simplex{M, T},
                new_vertex, dim) where {M, T}
    vxs = vertices(st, simplex, dim)
    diameter = max(diam(simplex), max_dist(st, vxs, new_vertex))

    index = 0
    new_vertex_done = false
    k = dim + 2
    i = 1
    coefficient = coef(simplex)
    while k > 0
        if !new_vertex_done && (k == 1 || new_vertex > vxs[i])
            new_vertex_done = true
            index += binomial(st, new_vertex - 1, k)
            coefficient = (k % 2 == 1 ? M - 1 : 1) * coefficient % M
            k -= 1
        else
            index += binomial(st, vxs[i] - 1, k)
            i += 1
            k -= 1
        end
    end
    Simplex{M}(diameter, index + 1, coefficient)
end

struct CoboundaryIterator{M, T, R<:ReductionState{M, T}}
    state   ::R
    simplex ::Simplex{M, T}
    dim     ::Int
end

Base.IteratorSize(::Type{CoboundaryIterator}) =
    Base.SizeUnknown()
Base.IteratorEltype(::Type{CoboundaryIterator}) =
    Base.HasEltype()
Base.eltype(::Type{CoboundaryIterator{M, T}}) where {M, T} =
    Simplex{M, T}

function Base.iterate(ni::CoboundaryIterator, i = n_vertices(ni.state))
    vxs = vertices(ni.state, ni.simplex, ni.dim)
    while i > 0 && !is_connected(ni.state, vxs, i)
        i -= 1
    end
    if i < 1
        nothing
    else
        coface(ni.state, ni.simplex, i, ni.dim), i - 1
    end
end

coboundary(st, simplex, dim) =
    CoboundaryIterator(st, simplex, dim)

# columns ================================================================================ #
const Column{M, T} =
    BinaryHeap{Simplex{M, T}, SimplexComparer}

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
    reduction_matrix  ::CompressedSparseMatrix{Simplex{M, T}}
    column_index      ::Dict{Int, Tuple{Int, Int}} # index(sx) => column index of reduction_matrix
                                       # should be sx => column index of reduction_matrix
    working_column    ::Column{M, T}
    reduction_entries ::Column{M, T}
    dim               ::Int
end

ReductionMatrices(st::ReductionState{M, T}, dim) where {M, T} =
    ReductionMatrices(st, CompressedSparseMatrix{Simplex{M, T}}(),
                      Dict{Int, Tuple{Int, Int}}(), Column{M, T}(), Column{M, T}(), dim)

"""
    add!(rm::ReductionMatrices, index)

Add column with column index `index` multiplied by the correct factor to `working_column`.
Also record the addition in `reduction_matrix`.
"""
# add mora vedet katere simplekse je dodal, ne njihovih coboundaryjev
function add!(rm::ReductionMatrices, idx, other_coef)
    λ = coef(pivot(rm.working_column) / other_coef)
    for simplex in rm.reduction_matrix[idx]
        push!(rm.reduction_entries, -simplex * λ)
        for coface in coboundary(rm.state, simplex, rm.dim)
            push!(rm.working_column, -coface * λ)
        end
    end
    pivot(rm.working_column)
end

function initialize!(rm::ReductionMatrices, column_simplex)
    empty!(rm.working_column.valtree)
    empty!(rm.reduction_entries.valtree)

    for coface in coboundary(rm.state, column_simplex, rm.dim)
        if diam(coface) == diam(column_simplex) && !haskey(rm.column_index, index(coface))
            return coface
        end
    end
    for coface in coboundary(rm.state, column_simplex, rm.dim)
        push!(rm.working_column, coface)
    end
    pivot(rm.working_column)
end

function reduce_working_column!(rm::ReductionMatrices, res, column_simplex)
    current_pivot = initialize!(rm, column_simplex)

    add_column!(rm.reduction_matrix)
    push!(rm.reduction_matrix, column_simplex)

    while !isnothing(current_pivot) && haskey(rm.column_index, index(current_pivot))
        old_pivot = current_pivot
        current_pivot = add!(rm, rm.column_index[index(current_pivot)]...)
        @assert old_pivot != current_pivot
    end
    if isnothing(current_pivot)
        push!(res, (birth, typemax(T)))
    else
        death = diam(current_pivot)
        birth = diam(column_simplex)
        if death > birth
            push!(res, (birth, death))
        end

        rm.column_index[index(current_pivot)] = (length(rm.reduction_matrix),
                                                 coef(current_pivot))
        current_entry = pop_pivot!(rm.reduction_entries)
        while !isnothing(current_entry)
            push!(rm.reduction_matrix, current_entry)
            current_entry = pop_pivot!(rm.reduction_entries)
        end
    end
end

"""
    compute_0_dim_pairs!(reduction_state, columns)

Compute 0-dimensional persistent homology using Kruskal's Algorithm.
"""
function compute_0_dim_pairs!(st::ReductionState{M, T}, simplices, columns) where {M, T}
    dset = IntDisjointSets(n_vertices(st))
    res = Tuple{T, T}[]

    for (l, (u, v)) in edges(st)
        push!(simplices, Simplex{M}(st, l, (u, v), 1))
        i = find_root(dset, u)
        j = find_root(dset, v)
        if i ≠ j
            union!(dset, i, j)
            if l > 0
                push!(res, (zero(T), T(l)))
            end
        else
            push!(columns, Simplex{M}(st, T(l), (u, v), 1))
        end
    end
    for _ in 1:num_groups(dset)
        push!(res, (zero(T), typemax(T)))
    end
    reverse!(columns)
    res
end

function compute_pairs!(rm::ReductionMatrices{M, T}, columns) where {M, T}
    res = Tuple{T, T}[]
    for column in columns
        reduce_working_column!(rm, res, column)
    end
    res
end

# move me.
function diam(st::ReductionState{M, T}, vertices) where {M, T}
    n = length(vertices)
    res = typemin(T)
    for i in 1:n, j in i+1:n
        d = dist(st, vertices[j], vertices[i])
        if d == 0
            return typemax(T)
        else
            res = max(res, d)
        end
    end
    res
end

function assemble_columns!(rm::ReductionMatrices{M, T}, simplices, columns) where {M, T}
    empty!(columns)
    new_simplices = Simplex{M, T}[]

    for simplex in simplices
        for coface in coboundary(rm.state, simplex, rm.dim)
            coface = set_coef(coface, 1)
            push!(new_simplices, coface)
            if !haskey(rm.column_index, index(coface))
                push!(columns, coface)
            end
        end
    end
    copy!(simplices, unique!(new_simplices))
    sort!(unique!(columns), lt=SimplexComparer(), rev=true)
    columns
end

ripserer(dists, dim_max=1, modulus=2) =
    ripserer(dists, dim_max, Val(modulus))

function ripserer(dists::AbstractMatrix{T}, dim_max, ::Val{M}) where {M, T}
    st = ReductionState{M}(dists, dim_max)
    res = Vector{Tuple{T, T}}[]
    simplices = Simplex{M, T}[]
    columns = Simplex{M, T}[]

    push!(res, compute_0_dim_pairs!(st, simplices, columns))

    for dim in 1:dim_max
        rm = ReductionMatrices(st, dim)
        push!(res, compute_pairs!(rm, columns))
        if dim < dim_max
            assemble_columns!(rm, simplices, columns)
        end
    end
    res
end
