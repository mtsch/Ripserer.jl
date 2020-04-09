# compressed sparse matrix =============================================================== #
"""
    CompressedSparseMatrix{T}

Compressed immutable sparse matrix data structure that supports efficient column insertion,
pushing to the last column via [`push!`](@ref) and iterating over columns.

It's up to the value type `T` to know about its row position.
"""
struct CompressedSparseMatrix{T}
    colptr::Vector{Int}
    nzval::Vector{T}
end

CompressedSparseMatrix{T}() where T =
    CompressedSparseMatrix(Int[1], T[])

function Base.show(io::IO, csm::CompressedSparseMatrix{T}) where T
    println(io, "CompressedSparseMatrix{$T}[")
    for i in 1:length(csm)
        println(io, "  $i: ", collect(csm[i]))
    end
    print(io, "]")
end

function Base.push!(csm::CompressedSparseMatrix, value)
    push!(csm.nzval, value)
    csm.colptr[end] += 1
    value
end

add_column!(csm::CompressedSparseMatrix) =
    push!(csm.colptr, csm.colptr[end])
Base.eltype(csm::CompressedSparseMatrix{T}) where T =
    T
Base.length(csm::CompressedSparseMatrix) =
    length(csm.colptr) - 1
Base.getindex(csm::CompressedSparseMatrix, i) =
    CSMColumnIterator(csm, i)

struct CSMColumnIterator{T}
    csm ::CompressedSparseMatrix{T}
    idx ::Int
end

Base.IteratorSize(::Type{CSMColumnIterator}) =
    Base.HasLength()
Base.IteratorEltype(::Type{CSMColumnIterator{T}}) where T =
    Base.HasEltype()
Base.eltype(::Type{CSMColumnIterator{T}}) where T =
    T
Base.length(ci::CSMColumnIterator) =
    ci.csm.colptr[ci.idx + 1] - ci.csm.colptr[ci.idx]

function Base.iterate(ci::CSMColumnIterator, i=1)
    colptr = ci.csm.colptr
    index = i + colptr[ci.idx] - 1
    if index ≥ colptr[ci.idx + 1]
        nothing
    else
        (ci.csm.nzval[index], i + 1)
    end
end

# columns ================================================================================ #
const Column{M, T} =
    BinaryMinHeap{Simplex{M, T}}

"""
    pop_pivot!(column)

Pop the pivot from `column`. If there are multiple simplices of with the same index on the
top of the column, sum them together. If they sum to 0 pop the next column.
"""
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

"""
    pivot(column)

Return the pivot of the column.
"""
function pivot(column::Column)
    pivot = pop_pivot!(column)
    if !isnothing(pivot)
        push!(column, pivot)
    end
    pivot
end

# reduction matrix ======================================================================= #
struct ReductionMatrix{M, T, C<:SimplicialComplex{M, T}}
    complex           ::C
    reduction_matrix  ::CompressedSparseMatrix{Simplex{M, T}}
    column_index      ::Dict{Int, Tuple{Int, Int}}
    working_column    ::Column{M, T}
    reduction_entries ::Column{M, T}
    dim               ::Int
end

ReductionMatrix(scx::SimplicialComplex{M, T}, dim) where {M, T} =
    ReductionMatrix(scx, CompressedSparseMatrix{Simplex{M, T}}(),
                    Dict{Int, Tuple{Int, Int}}(),
                    Column{M, T}(), Column{M, T}(), dim)

"""
    add!(rm::ReductionMatrix, index)

Add column with column `index` multiplied by the correct factor to `rm.working_column`.
Also record the addition in `rm.reduction_entries`.
"""
function add!(rm::ReductionMatrix, idx, other_coef)
    λ = coef(pivot(rm.working_column) / other_coef)
    for simplex in rm.reduction_matrix[idx]
        push!(rm.reduction_entries, -simplex * λ)
        for coface in coboundary(rm.complex, simplex, rm.dim)
            if diam(coface) ≤ threshold(rm.complex)
                push!(rm.working_column, -coface * λ)
            end
        end
    end
    pivot(rm.working_column)
end

"""
    initialize!(rm::ReductionMatrix, column_simplex)

Initialize `rm.working_column` by emptying it and `reduction_entries` and pushing the
coboundary of `column_simplex` to `rm.working_column`.
"""
function initialize!(rm::ReductionMatrix, column_simplex)
    empty!(rm.working_column.valtree)
    empty!(rm.reduction_entries.valtree)

    for coface in coboundary(rm.complex, column_simplex, rm.dim)
        if diam(coface) ≤ threshold(rm.complex)
            if diam(coface) == diam(column_simplex) && !haskey(rm.column_index, index(coface))
                empty!(rm.working_column.valtree)
                return coface
            end
            push!(rm.working_column, coface)
        end
    end

    pivot(rm.working_column)
end

"""
    reduce_working_column!(rm::ReductionMatrix, res, column_simplex)

Reduce the working column by adding other columns to it until it has the lowest pivot.
Record resulting persistence intervals in `res`.
"""
function reduce_working_column!(rm::ReductionMatrix{M, T}, res, column_simplex) where {M, T}
    current_pivot = initialize!(rm, column_simplex)

    add_column!(rm.reduction_matrix)
    push!(rm.reduction_matrix, column_simplex)

    while !isnothing(current_pivot) && haskey(rm.column_index, index(current_pivot))
        current_pivot = add!(rm, rm.column_index[index(current_pivot)]...)
    end
    if isnothing(current_pivot)
        push!(res, (diam(column_simplex), typemax(T)))
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
function compute_0_dim_pairs!(scx::SimplicialComplex{M, T}, columns) where {M, T}
    dset = IntDisjointSets(length(scx))
    res = Tuple{T, T}[]

    for (l, (u, v)) in edges(scx)
        i = find_root!(dset, u)
        j = find_root!(dset, v)
        if i ≠ j
            union!(dset, i, j)
            if l > 0
                push!(res, (zero(T), T(l)))
            end
        else
            push!(columns, Simplex{M}(scx, T(l), (u, v), 1))
        end
    end
    for _ in 1:num_groups(dset)
        push!(res, (zero(T), typemax(T)))
    end
    reverse!(columns)
    res
end

"""
    compute_pairs!(rm::ReductionMatrix, columns)

Compute persistence pairs by reducing `columns` (list of simplices).
"""
function compute_pairs!(rm::ReductionMatrix{M, T}, columns) where {M, T}
    res = Tuple{T, T}[]
    for column in columns
        reduce_working_column!(rm, res, column)
    end
    res
end

"""
    assemble_columns!(rm::ReductionMatrix, columns)

Assemble columns that need to be reduced in the next dimension. Apply clearing optimization.
"""
function assemble_columns!(rm::ReductionMatrix{M, T}, columns) where {M, T}
    empty!(columns)
    n_simplices = binomial(rm.complex, length(rm.complex), rm.dim + 2)
    S = eltype(rm.complex)

    for idx in 1:n_simplices
        if !haskey(rm.column_index, idx)
            sx = S(diam(rm.complex, vertices(rm.complex, idx, rm.dim+1)), idx, 1)
            if diam(sx) ≤ threshold(rm.complex)
                push!(columns, sx)
            end
        end
    end
    sort!(columns, rev=true)
    columns
end

"""
    ripserer(dists::AbstractMatrix{T}; dim_max=1, modulus=2, threshold=typemax(T))

    ripserer(scx::SimplicialComplex)

Compute the persistent homology of metric space represented by `dists` or filtration
represented by `scx`.

# Settings

`dim_max`: compute persistent homology up to this dimension.
`modulus`: compute persistent homology with coefficients in the prime field of inether
           mod `modulus`.
`threshold`: compute persistent homology up to diameter smaller than threshold.
"""
ripserer(dists::AbstractMatrix{T}; dim_max=1, modulus=2, threshold=typemax(T)) where T =
    ripserer(RipsComplex{modulus}(dists, dim_max, threshold))

function ripserer(scx::SimplicialComplex{M, T}) where {M, T}
    res = Vector{Tuple{T, T}}[]
    columns = eltype(scx)[]

    push!(res, compute_0_dim_pairs!(scx, columns))

    for dim in 1:dim_max(scx)
        rm = ReductionMatrix(scx, dim)
        push!(res, compute_pairs!(rm, columns))
        if dim < dim_max(scx)
            assemble_columns!(rm, columns)
        end
    end
    res
end
