# compressed sparse matrix =============================================================== #
"""
    CompressedSparseMatrix{T}

Compressed immutable sparse matrix data structure that supports efficient column insertion,
pushing to the last column via [`Base.push!`](@ref) and iterating over columns.

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
"""
    Column{S<:AbstractSimplex}

Wrapper around `BinaryMinHeap{S}`. Acts like a heap of simplices, where simplices with the
same index are summed together and simplices with coefficient value `0` are ignored.
"""
struct Column{S<:AbstractSimplex}
    heap::BinaryMinHeap{S}

    Column{S}() where S<:AbstractSimplex =
        new{S}(BinaryMinHeap{S}())
end

Base.empty!(col::Column) =
    empty!(col.heap.valtree)
Base.isempty(col::Column) =
    isempty(col.heap)

"""
    pop_pivot!(column::Column)

Pop the pivot from `column`. If there are multiple simplices with the same index on the top
of the column, sum them together. If they sum to 0, pop the next column. Return
`nothing` when column is empty.
"""
function pop_pivot!(column::Column)
    isempty(column) && return nothing
    heap = column.heap

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

"""
    pivot(column)

Return the pivot of the column.
"""
function pivot(column::Column)
    heap = column.heap
    pivot = pop_pivot!(column)
    if !isnothing(pivot)
        push!(column, pivot)
    end
    pivot
end

function Base.push!(column::Column{S}, sx::S) where S
    heap = column.heap
    if !isempty(heap) && index(top(heap)) == index(sx)
        heap.valtree[1] += sx
    else
        push!(heap, sx)
    end
end

# reduction matrix ======================================================================= #
"""
    ReductionMatrix

This structure represents the reduction matrix in the current dimension. A new one is
created for every dimension.

# Fields:

* `filtration`: the filtration we are analyzing.
* `coboundary`: a `Coboundary` object, used to find coboundaries and vertices.
* `reduction_matrix`: the reduction matrix. Each column of the matrix records the operations
  that were performed when reducing the column.
* `column_index`: a `Dict` that maps pivot index to its position in `reduction_matrix` and
  coefficient.
* `working_column`: the current working column, the column we are currently reducing.
* `reduction_entries`: this is where we record which simplices we added to the working
  column.
* `dim`: the current dimension.
"""
struct ReductionMatrix{
    T, I, S<:AbstractSimplex{I, T}, F<:AbstractFiltration{T, S}, C<:Coboundary{S, F}
}
    filtration        ::F
    coboundary        ::C
    reduction_matrix  ::CompressedSparseMatrix{S}
    column_index      ::Dict{Int, Tuple{Int, I}}
    working_column    ::Column{S}
    reduction_entries ::Column{S}
    dim               ::Int
end

function ReductionMatrix(
    coboundary::Coboundary{S, F},
    filtration::F,
    dim,
) where {T, I, S<:AbstractSimplex{I, T}, F<:AbstractFiltration{T, S}}
    ReductionMatrix(
        filtration,
        coboundary,
        CompressedSparseMatrix{S}(),
        Dict{Int, Tuple{Int, I}}(),
        Column{S}(),
        Column{S}(),
        dim,
    )
end

"""
    add!(rm::ReductionMatrix, index)

Add column with column `index` multiplied by the correct factor to `rm.working_column`.
Also record the addition in `rm.reduction_entries`.
"""
function add!(rm::ReductionMatrix, current_pivot, idx, other_coef)
    λ = -coef(current_pivot / other_coef)
    for simplex in rm.reduction_matrix[idx]
        push!(rm.reduction_entries, simplex * λ)
        for coface in rm.coboundary(simplex, rm.dim)
            push!(rm.working_column, coface * λ)
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
    empty!(rm.working_column)
    empty!(rm.reduction_entries)

    for coface in rm.coboundary(column_simplex, rm.dim)
        if diam(coface) == diam(column_simplex) && !haskey(rm.column_index, index(coface))
            empty!(rm.working_column)
            return coface
        end
        push!(rm.working_column, coface)
    end
    pivot(rm.working_column)
end

"""
    reduce_working_column!(rm::ReductionMatrix, res, column_simplex)

Reduce the working column by adding other columns to it until it has the lowest pivot or is
reduced. Record resulting persistence intervals in `res`.
"""
function reduce_working_column!(rm::ReductionMatrix, res, column_simplex)
    current_pivot = initialize!(rm, column_simplex)

    add_column!(rm.reduction_matrix)
    push!(rm.reduction_matrix, column_simplex)

    while !isnothing(current_pivot) && haskey(rm.column_index, index(current_pivot))
        current_pivot = add!(rm, current_pivot, rm.column_index[index(current_pivot)]...)
    end
    if isnothing(current_pivot)
        push!(res, (diam(column_simplex), ∞))
    else
        birth = diam(column_simplex)
        death = diam(current_pivot)
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
    rm
end

"""
    compute_0_dim_pairs!(filtration, columns)

Compute 0-dimensional persistent homology using Kruskal's Algorithm.
If `filtration` is sparse, also return a vector of all 1-simplices with diameter below
threshold.
"""
function compute_0_dim_pairs!(coboundary::Coboundary, columns)
    filtration = coboundary.filtration
    T = dist_type(filtration)
    dset = IntDisjointSets(n_vertices(filtration))
    res = Tuple{T, Union{T, Infinity}}[]
    # We only collect simplices if the filtration is sparse.
    simplices = issparse(filtration) ? eltype(filtration)[] : nothing

    for (l, (u, v)) in edges(filtration)
        i = find_root!(dset, u)
        j = find_root!(dset, v)
        issparse(filtration) &&
            push!(simplices, eltype(filtration)(T(l), index(coboundary, (u, v)), 1))
        if i ≠ j
            union!(dset, i, j)
            if l > 0
                push!(res, (zero(T), T(l)))
            end
        else
            push!(columns, eltype(filtration)(T(l), index(coboundary, (u, v)), 1))
        end
    end
    for _ in 1:num_groups(dset)
        push!(res, (zero(T), ∞))
    end
    reverse!(columns)
    res, simplices
end

"""
    compute_pairs!(rm::ReductionMatrix, columns)

Compute persistence pairs by reducing `columns` (list of simplices).
"""
function compute_pairs!(rm::ReductionMatrix, columns)
    T = dist_type(rm.filtration)
    res = Tuple{T, Union{T, Infinity}}[]
    for column in columns
        reduce_working_column!(rm, res, column)
    end
    res
end

"""
    assemble_columns!(rm::ReductionMatrix, columns, simplices)

Assemble columns that need to be reduced in the next dimension. Apply clearing optimization.
The algorithm used depends on whether the filtration is sparse or not. When it is, we
collect columns by only looking through the cofaces of simplices from the previous
dimension. When it's not, we go through all valid simplex indices.
"""
# This method is used when filtration is _not_ sparse.
function assemble_columns!(rm::ReductionMatrix, columns, ::Nothing)
    empty!(columns)
    n_simplices = binomial(rm.coboundary, n_vertices(rm.filtration), rm.dim + 2)
    S = eltype(rm.filtration)

    simplices = trues(n_simplices)
    for k in keys(rm.column_index)
        @inbounds simplices[k] = false
    end
    sizehint!(columns, sum(simplices))
    for idx in 1:n_simplices
        if simplices[idx]
            diameter = diam(rm.filtration, vertices(rm.coboundary, idx, rm.dim + 1))
            #sx = S(diam(rm.filtration, vertices(rm.filtration, idx, rm.dim + 1)), idx, 1)
            if diameter < ∞
                push!(columns, S(diameter, idx, 1))
            end
        end
    end
    sort!(columns, rev=true)
    columns
end

function assemble_columns!(rm::ReductionMatrix, columns, simplices)
    empty!(columns)
    new_simplices = eltype(simplices)[]

    for simplex in simplices
        for coface in rm.coboundary(simplex, rm.dim, Val(false))
            push!(new_simplices, coface)
            if !haskey(rm.column_index, index(coface))
                push!(columns, coface)
            end
        end
    end
    copy!(simplices, new_simplices)
    sort!(columns, rev=true)
    columns
end

"""
    ripserer(dists::AbstractMatrix{T}; dim_max=1, modulus=2, threshold=typemax(T))
    ripserer(points, metric; dim_max=1, modulus=2, threshold=typemax(T))

Compute the persistent homology of metric space represented by `dists` or `points` and
`metric`.

# Keyoword Arguments

* `dim_max`: compute persistent homology up to this dimension.
* `modulus`: compute persistent homology with coefficients in the prime field of integers
             mod `modulus`.
* `threshold`: compute persistent homology up to diameter smaller than threshold.
  Defaults to radius of input space.
* `sparse`: if `true`, use `SparseRipsFiltration`. Defaults to `issparse(dists)`.
"""
function ripserer(dists::AbstractMatrix; sparse=issparse(dists), dim_max=1, kwargs...)
    if issparse(dists)
        ripserer(SparseRipsFiltration(dists; kwargs...), dim_max=dim_max)
    else
        ripserer(RipsFiltration(dists; kwargs...), dim_max=dim_max)
    end
end
function ripserer(points, metric; sparse=false, dim_max=1, kwargs...)
    if sparse
        ripserer(SparseRipsFiltration(points, metric; kwargs...), dim_max=dim_max)
    else
        ripserer(RipsFiltration(points, metric; kwargs...), dim_max=dim_max)
    end
end

"""
    ripserer(filtration::AbstractFiltration)

Compute persistent homology from `filtration` object.
"""
function ripserer(filtration::AbstractFiltration{T}; dim_max=1) where T
    res = Vector{Tuple{T, Union{T, Infinity}}}[]
    columns = eltype(filtration)[]
    coboundary = Coboundary(filtration, dim_max)

    res_0, simplices = compute_0_dim_pairs!(coboundary, columns)
    push!(res, res_0)

    for dim in 1:dim_max
        rm = ReductionMatrix(coboundary, filtration, dim)
        push!(res, compute_pairs!(rm, columns))
        if dim < dim_max
            assemble_columns!(rm, columns, simplices)
        end
    end
    res
end
