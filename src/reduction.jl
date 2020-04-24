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
    Ds, Dc, I, S<:AbstractSimplex{Ds, I}, C<:AbstractSimplex{Dc, I}, F<:AbstractFiltration
}
    filtration        ::F
    reduction_matrix  ::CompressedSparseMatrix{S}
    column_index      ::Dict{Int, Tuple{Int, I}}
    working_column    ::Column{C}
    reduction_entries ::Column{S}

    function ReductionMatrix(
        filtration::F,
        reduction_matrix::CompressedSparseMatrix{S},
        column_index::Dict{Int, Tuple{Int, I}},
        working_column::Column{C},
        reduction_entries::Column{S},
    ) where {
        Ds, Dc, I,
        S<:AbstractSimplex{Ds, I},
        C<:AbstractSimplex{Dc, I},
        F<:AbstractFiltration,
    }
        Dc == Ds + 1 || error("dimension mismatch Dc=$Dc, Ds=$Ds")
        new{Ds, Dc, I, S, C, F}(filtration, reduction_matrix, column_index, working_column, reduction_entries)
    end
end

function ReductionMatrix(
    filtration::F,
    ::Type{S},
) where {D, I, S<:AbstractSimplex{D, I}, F<:AbstractFiltration}
    C = coface_type(S)
    ReductionMatrix(
        filtration,
        CompressedSparseMatrix{S}(),
        Dict{Int, Tuple{Int, I}}(),
        Column{C}(),
        Column{S}(),
    )
end

simplex_type(rm::ReductionMatrix{<:Any, <:Any, <:Any, S}) where S =
    S

"""
    add!(rm::ReductionMatrix, index)

Add column with column `index` multiplied by the correct factor to `rm.working_column`.
Also record the addition in `rm.reduction_entries`.
"""
function add!(rm::ReductionMatrix, current_pivot, idx, other_coef)
    λ = -coef(current_pivot / other_coef)
    for simplex in rm.reduction_matrix[idx]
        push!(rm.reduction_entries, simplex * λ)
        for coface in coboundary(rm.filtration, simplex)
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
function initialize!(rm::ReductionMatrix, column_simplex::AbstractSimplex)
    empty!(rm.working_column)
    empty!(rm.reduction_entries)

    for coface in coboundary(rm.filtration, column_simplex)
        if diam(coface) == diam(column_simplex) && !haskey(rm.column_index, index(coface))
            empty!(rm.working_column)
            return coface
        end
        push!(rm.working_column, coface)
    end
    # Pivot always exists, so we use :: to make sure the correct type is inferred.
    pivot(rm.working_column)::coface_type(column_simplex)
end

"""
    reduce_working_column!(rm::ReductionMatrix, res, column_simplex)

Reduce the working column by adding other columns to it until it has the lowest pivot or is
reduced. Record it in the reduction matrix and return the persistence interval.
"""
function reduce_working_column!(rm::ReductionMatrix, column_simplex::AbstractSimplex)
    current_pivot = initialize!(rm, column_simplex)

    add_column!(rm.reduction_matrix)
    push!(rm.reduction_matrix, column_simplex)

    while !isnothing(current_pivot) && haskey(rm.column_index, index(current_pivot))
        current_pivot = add!(rm, current_pivot, rm.column_index[index(current_pivot)]...)
    end
    if isnothing(current_pivot)
        death = ∞
    else
        rm.column_index[index(current_pivot)] = (length(rm.reduction_matrix),
                                                 coef(current_pivot))
        current_entry = pop_pivot!(rm.reduction_entries)
        while !isnothing(current_entry)
            push!(rm.reduction_matrix, current_entry)
            current_entry = pop_pivot!(rm.reduction_entries)
        end
        death = diam(current_pivot)

    end
    birth = diam(column_simplex)
    (birth, death)
end

"""
    compute_pairs!(rm::ReductionMatrix, columns)

Compute persistence intervals by reducing `columns` (list of simplices).
"""
function compute_intervals!(rm::ReductionMatrix, columns)
    T = dist_type(rm.filtration)
    intervals = Tuple{T, Union{T, Infinity}}[]
    for column in columns
        birth, death = reduce_working_column!(rm, column)
        if death > birth
            push!(intervals, (birth, death))
        end
    end
    intervals
end

"""
    assemble_columns!(rm::ReductionMatrix, columns, simplices)

Assemble columns that need to be reduced in the next dimension. Apply clearing optimization.
The algorithm used depends on whether the filtration is sparse or not. When it is, we
collect columns by only looking through the cofaces of simplices from the previous
dimension. When it's not, we go through all valid simplex indices.
"""
# This method is used when filtration is _not_ sparse.
function assemble_columns!(rm::ReductionMatrix{D}, ::Nothing) where D
    S = coface_type(simplex_type(rm))

    n_simplices = binomial(n_vertices(rm.filtration), D + 2)
    simplices = trues(n_simplices)
    for k in keys(rm.column_index)
        simplices[k] = false
    end

    columns = S[]
    sizehint!(columns, sum(simplices))
    for idx in 1:n_simplices
        if simplices[idx]
            diameter = diam(rm.filtration, vertices(idx, Val(D+1)))
            if diameter < ∞
                push!(columns, S(diameter, idx, 1))
            end
        end
    end
    sort!(columns, rev=true)
    columns, nothing
end

function assemble_columns!(rm::ReductionMatrix, simplices)
    S = coface_type(simplex_type(rm))
    columns = S[]
    new_simplices = S[]

    for simplex in simplices
        for coface in coboundary(rm.filtration, simplex, Val(false))
            push!(new_simplices, coface)
            if !haskey(rm.column_index, index(coface))
                push!(columns, coface)
            end
        end
    end
    sort!(columns, rev=true)
    columns, new_simplices
end

"""
    zeroth_intervals(filtration)

Compute 0-dimensional persistent homology using Kruskal's Algorithm.
If `filtration` is sparse, also return a vector of all 1-simplices with diameter below
threshold.
"""
function zeroth_intervals(filtration)
    T = dist_type(filtration)
    dset = IntDisjointSets(n_vertices(filtration))
    intervals = Tuple{T, Union{T, Infinity}}[]
    # We only collect simplices if the filtration is sparse.
    simplices = issparse(filtration) ? edge_type(filtration)[] : nothing
    columns = edge_type(filtration)[]

    for sx in edges(filtration)
        u, v = vertices(sx)
        i = find_root!(dset, u)
        j = find_root!(dset, v)
        issparse(filtration) &&
            push!(simplices, sx)
        if i ≠ j
            union!(dset, i, j)
            if diam(sx) > 0
                push!(intervals, (zero(T), diam(sx)))
            end
        else
            push!(columns, sx)
        end
    end
    for _ in 1:num_groups(dset)
        push!(intervals, (zero(T), ∞))
    end
    reverse!(columns)
    intervals, columns, simplices
end

function nth_intervals(filtration,
                       columns::Vector{S},
                       simplices;
                       next=true,
                       ) where S<:AbstractSimplex

    rm = ReductionMatrix(filtration, S)
    intervals = compute_intervals!(rm, columns)
    if next
        (intervals, assemble_columns!(rm, simplices)...)
    else
        (intervals, nothing, nothing)
    end
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
ripserer(filtration::AbstractFiltration; dim_max=1) =
    ripserer(filtration, Val(dim_max))

@generated function ripserer(filtration::AbstractFiltration{T}, ::Val{D}) where {T, D}
    # We unroll the loop over 1:D to ensure type stability.
    # Generated code looks something like:
    #
    # ints_0, cols_1, sxs_1 = zeroth_intervals(filtration)
    # ints_1, cols_2, sxs_2 = nth_itervals(filtration, cols_1, sxs_1)
    # ...
    # ints_D, _, _ = nth_itervals(filtration, cols_D-1, sxs_D-1, next=false)
    #
    # (ints_0, ints_1, ..., ints_D)

    ints = [Symbol("ints_", i) for i in 1:D]
    cols = [Symbol("cols_", i) for i in 1:D]
    sxs = [Symbol("sxs_", i) for i in 1:D]

    expr = quote
        ints_0, $(cols[1]), $(sxs[1]) = zeroth_intervals(filtration)
    end

    for i in 1:D-1
        expr = quote
            $expr
            $(ints[i]), $(cols[i+1]), $(sxs[i+1]) =
                nth_intervals(filtration, $(cols[i]), $(sxs[i]))
        end
    end

    quote
        $expr
        $(ints[D]), _, _ =
            nth_intervals(filtration, $(cols[D]), $(sxs[D]), next=false)

        Array{Tuple{T, Union{T, Infinity}}}[ints_0, $(ints...)]
    end
end
