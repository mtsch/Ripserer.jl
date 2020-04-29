"""
    ReductionMatrix{S}

TODO
"""
struct ReductionMatrix{S}
    column_index ::Dict{Int, Int}
    colptr       ::Vector{Int}
    nzval        ::Vector{S}
end

function Base.show(io::IO, ::MIME"text/plain", rm::ReductionMatrix{S}) where S
    println(io, "ReductionMatrix{$S}[")
    for i in keys(rm.column_index)
        println(io, "  $i: ", collect(rm[i]))
    end
    print(io, "]")
end

ReductionMatrix{S}() where S =
    ReductionMatrix(Dict{Int, Int}(), Int[1], S[])

has_column(rm::ReductionMatrix, i) =
    haskey(rm.column_index, i)

function insert_column!(rm::ReductionMatrix, i)
    rm.column_index[i] = length(rm.colptr)
    push!(rm.colptr, rm.colptr[end])
    rm
end

function Base.push!(rm::ReductionMatrix, value)
    push!(rm.nzval, value)
    rm.colptr[end] += 1
    value
end

function Base.sizehint!(rm::ReductionMatrix, n)
    sizehint!(rm.column_index, n)
    sizehint!(rm.colptr, n)
    sizehint!(rm.nzval, n)
end

Base.eltype(rm::ReductionMatrix{T}) where T =
    T
Base.length(rm::ReductionMatrix) =
    length(rm.colptr) - 1
Base.lastindex(rm::ReductionMatrix) =
    length(rm.colptr) - 1
Base.getindex(rm::ReductionMatrix, i) =
    RMColumnIterator(rm, rm.column_index[i])
"""
    RMColumnIterator{S}

An iterator over a column of a `ReductionMatrix{S}`.
"""
struct RMColumnIterator{S}
    rm  ::ReductionMatrix{S}
    idx ::Int
end

Base.IteratorSize(::Type{RMColumnIterator}) =
    Base.HasLength()
Base.IteratorEltype(::Type{RMColumnIterator{T}}) where T =
    Base.HasEltype()
Base.eltype(::Type{RMColumnIterator{T}}) where T =
    T
Base.length(ci::RMColumnIterator) =
    ci.rm.colptr[ci.idx + 1] - ci.rm.colptr[ci.idx]

function Base.iterate(ci::RMColumnIterator, i=1)
    colptr = ci.rm.colptr
    index = i + colptr[ci.idx] - 1
    if index ≥ colptr[ci.idx + 1]
        nothing
    else
        (ci.rm.nzval[index], i + 1)
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
    move!([dst, ]col::Column; times=1)

Move contents of column into `dst` by repeatedly calling `push!`. `dst` defaults to `S[]`.
Multipy all elements that are moved by `times`.
"""
move!(col::Column{S}) where S =
    move!(S[], col)
function move!(dst, col::Column{S}; times=one(S)) where S
    pivot = pop_pivot!(col)
    while !isnothing(pivot)
        push!(dst, times * pivot)
        pivot = pop_pivot!(col)
    end
    dst
end

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
    column
end

# disjointset with birth ================================================================= #
"""
    DisjointSetsWithBirth{T}

Almost identical to `DataStructures.IntDisjointSets`, but keeps track of vertex birth times.
Has no `num_groups` method.
"""
struct DisjointSetsWithBirth{T}
    parents ::Vector{Int}
    ranks   ::Vector{Int}
    births  ::Vector{T}

    function DisjointSetsWithBirth(births::AbstractVector{T}) where T
        n = length(births)
        new{T}(collect(1:n), fill(0, n), copy(births))
    end
end

function DataStructures.find_root!(s::DisjointSetsWithBirth, x)
    parents = s.parents
    p = parents[x]
    @inbounds if parents[p] != p
        parents[x] = p = find_root!(s, p)
    end
    p
end

function Base.union!(s::DisjointSetsWithBirth, x, y)
    parents = s.parents
    xroot = find_root!(s, x)
    yroot = find_root!(s, y)
    xroot != yroot ? root_union!(s, xroot, yroot) : xroot
end

function DataStructures.root_union!(s::DisjointSetsWithBirth, x, y)
    parents = s.parents
    rks = s.ranks
    births = s.births
    @inbounds xrank = rks[x]
    @inbounds yrank = rks[y]

    if xrank < yrank
        x, y = y, x
    elseif xrank == yrank
        rks[x] += 1
    end
    @inbounds parents[y] = x
    @inbounds births[x] = min(births[x], births[y])
    x
end

birth(dset::DisjointSetsWithBirth, i) =
    dset.births[i]

# reduction matrix ======================================================================= #
"""
    ReductionState

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
struct ReductionState{
    Ds, Dc, I, S<:AbstractSimplex{Ds, I}, C<:AbstractSimplex{Dc, I}, F<:AbstractFiltration
}
    filtration        ::F
    reduction_matrix  ::ReductionMatrix{S}
    working_column    ::Column{C}
    reduction_entries ::Column{S}

    function ReductionState(
        filtration        ::F,
        reduction_matrix  ::ReductionMatrix{S},
        working_column    ::Column{C},
        reduction_entries ::Column{S},
    ) where {
        Ds, Dc, I,
        S<:AbstractSimplex{Ds, I},
        C<:AbstractSimplex{Dc, I},
        F<:AbstractFiltration,
    }
        Dc == Ds + 1 || error("dimension mismatch Dc=$Dc, Ds=$Ds")
        new{Ds, Dc, I, S, C, F}(
            filtration,
            reduction_matrix,
            working_column,
            reduction_entries
        )
    end
end

function ReductionState(
    filtration::F,
    ::Type{S},
) where {D, I, S<:AbstractSimplex{D, I}, F<:AbstractFiltration}
    C = coface_type(S)
    ReductionState(
        filtration,
        ReductionMatrix{S}(),
        Column{C}(),
        Column{S}(),
    )
end

simplex_type(rs::ReductionState{<:Any, <:Any, <:Any, S}) where S =
    S

"""
    add!(rs::ReductionState, index)

Add column with column `index` multiplied by the correct factor to `rs.working_column`.
Also record the addition in `rs.reduction_entries`.
"""
function add!(rs::ReductionState, current_pivot)
    λ = -coef(current_pivot)
    for simplex in rs.reduction_matrix[index(current_pivot)]
        push!(rs.reduction_entries, λ * simplex)
        for coface in coboundary(rs.filtration, simplex)
            push!(rs.working_column, λ * coface)
        end
    end
    pivot(rs.working_column)
end

"""
    initialize!(rs::ReductionState, column_simplex)

Initialize `rs.working_column` by emptying it and `reduction_entries` and pushing the
coboundary of `column_simplex` to `rs.working_column`.
"""
function initialize!(rs::ReductionState, column_simplex::AbstractSimplex)
    empty!(rs.working_column)
    empty!(rs.reduction_entries)

    for coface in coboundary(rs.filtration, column_simplex)
        if diam(coface) == diam(column_simplex) && !has_column(rs.reduction_matrix, index(coface))
            empty!(rs.working_column)
            return coface
        end
        push!(rs.working_column, coface)
    end
    pivot(rs.working_column)
end

"""
    reduce_working_column!(rs::ReductionState, res, column_simplex)

Reduce the working column by adding other columns to it until it has the lowest pivot or is
reduced. Record it in the reduction matrix and return the persistence interval.
"""
function reduce_working_column!(
    rs::ReductionState,
    column_simplex::AbstractSimplex,
    ::Val{cocycles},
) where cocycles
    current_pivot = initialize!(rs, column_simplex)

    while !isnothing(current_pivot) && has_column(rs.reduction_matrix, index(current_pivot))
        current_pivot = add!(rs, current_pivot)
    end
    if isnothing(current_pivot)
        death = ∞
    else
        insert_column!(rs.reduction_matrix, index(current_pivot))
        push!(rs.reduction_entries, column_simplex)
        move!(rs.reduction_matrix, rs.reduction_entries, times=inv(coef(current_pivot)))
        death = diam(current_pivot)
    end
    birth = diam(column_simplex)
    if cocycles
        cocycle = collect(rs.reduction_matrix[index(current_pivot)])
        PersistenceInterval(birth, death, cocycle)
    else
        PersistenceInterval(birth, death)
    end
end

"""
    compute_pairs!(rs::ReductionState, columns)

Compute persistence intervals by reducing `columns`, a collection of simplices.
"""
function compute_intervals!(
    rs::ReductionState, columns, ratio, ::Val{cocycles},
) where cocycles
    T = dist_type(rs.filtration)
    if cocycles
        C = eltype(columns)
        intervals = PersistenceInterval{T, Vector{C}}[]
    else
        intervals = PersistenceInterval{T, Nothing}[]
    end
    for column in columns
        interval = reduce_working_column!(rs, column, Val(cocycles))
        if death(interval) > birth(interval) * ratio
            push!(intervals, interval)
        end
    end
    sort!(PersistenceDiagram(dim(eltype(columns)), intervals))
end

"""
    assemble_columns!(rs::ReductionState, columns, simplices)

Assemble columns that need to be reduced in the next dimension. Apply clearing optimization.
"""
function assemble_columns!(rs::ReductionState, simplices)
    S = coface_type(simplex_type(rs))
    columns = S[]
    new_simplices = S[]

    for simplex in simplices
        for coface in coboundary(rs.filtration, simplex, Val(false))
            push!(new_simplices, coface)
            if !has_column(rs.reduction_matrix, index(coface))
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
function zeroth_intervals(filtration, ratio=1)
    T = dist_type(filtration)
    dset = DisjointSetsWithBirth([birth(filtration, v) for v in 1:n_vertices(filtration)])
    intervals = PersistenceInterval{T, Nothing}[]
    simplices = edge_type(filtration)[]
    columns = edge_type(filtration)[]

    for sx in edges(filtration)
        u, v = vertices(sx)
        i = find_root!(dset, u)
        j = find_root!(dset, v)
        push!(simplices, sx)
        if i ≠ j
            # According to the elder rule, the vertex with the lower birth will fall
            # into a later interval.
            interval = PersistenceInterval(max(birth(dset, i), birth(dset, j)), diam(sx))
            if death(interval) > birth(interval) * ratio
                push!(intervals, interval)
            end
            union!(dset, i, j)
        else
            push!(columns, sx)
        end
    end
    for v in 1:n_vertices(filtration)
        if find_root!(dset, v) == v
            push!(intervals, PersistenceInterval(birth(dset, v), ∞))
        end
    end
    reverse!(columns)
    sort!(PersistenceDiagram(0, intervals)), columns, simplices
end

"""
    nth_intervals(filtration, columns, simplices; next=true)

Compute the ``n``-th intervals of persistent cohomology. The ``n`` is determined from the
`eltype` of `columns`. If `next` is `true`, assemble columns for the next dimension.
"""
function nth_intervals(
    filtration, columns::Vector{S}, simplices, ratio, ::Val{cocycles}; next=true,
) where {S<:AbstractSimplex, cocycles}

    rs = ReductionState(filtration, S)
    sizehint!(rs.reduction_matrix, length(columns))
    intervals = compute_intervals!(rs, columns, ratio, Val(cocycles))
    if next
        (intervals, assemble_columns!(rs, simplices)...)
    else
        (intervals, nothing, nothing)
    end
end

"""
    ripserer(dists::AbstractMatrix{T}; kwargs...)
    ripserer(points; metric=Euclidean(), births, kwargs...)

Compute the persistent homology of metric space represented by `dists` or `points` and
`metric`. `points` must be an array of bitstypes, such as `NTuple`s or `SVectors`.

# Keyoword Arguments

* `dim_max`: compute persistent homology up to this dimension. Defaults to `1`.
* `modulus`: compute persistent homology with coefficients in the prime field of integers
  mod `modulus`. Defaults to `2`.
* `threshold`: compute persistent homology up to diameter smaller than threshold.
  For Rips filtrations, it defaults to radius of input space.
* `sparse`: if `true`, use `SparseRipsFiltration`. Defaults to `false || issparse(dists)`.
* `ratio`: only keep intervals with `death(interval) > birth(interval) * ratio`.
  Defaults to `1`.
* `cocycles`: if `true`, return representative cocycles along with persistence
  intervals. Defaults to `false`.
* `metric`: when calculating persistent homology from points, any metric from
  [`Distances.jl`](https://github.com/JuliaStats/Distances.jl) can be used. Defaults to
  `Euclidean()`.
* `births`: when calculating persistent homology from points, births can be used to add
  birth times to vertices. Defaults to all births equal to `0`.
"""
function ripserer(
    dists::AbstractMatrix;
    dim_max=1,
    sparse=false || issparse(dists),
    ratio=1,
    cocycles=false,
    kwargs...
)
    if sparse
        filtration = SparseRipsFiltration(dists; kwargs...)
    else
        filtration = RipsFiltration(dists; kwargs...)
    end
    ripserer(filtration; dim_max=dim_max, cocycles=cocycles, ratio=ratio)
end

function ripserer(points; metric=Euclidean(), births=nothing, kwargs...)
    dists = distances(metric, points, births)
    ripserer(dists; kwargs...)
end

"""
    ripserer(filtration::AbstractFiltration; dim_max=1)

Compute persistent homology from `filtration` object.
"""
ripserer(filtration::AbstractFiltration; dim_max=1, cocycles=false, ratio=1) =
    ripserer(filtration, ratio, Val(dim_max), Val(cocycles))

function ripserer(filtration::AbstractFiltration, ratio, ::Val{0}, ::Val{<:Any})
    diagram, _, _ = zeroth_intervals(filtration, ratio)
    [diagram]
end

@generated function ripserer(
    filtration::AbstractFiltration, ratio, ::Val{D}, ::Val{C},
) where {D, C}
    # We unroll the loop over 1:D to ensure type stability.
    # Generated code looks something like:
    #
    # ints_0, cols_1, sxs_1 = zeroth_intervals(filtration)
    # ints_1, cols_2, sxs_2 = nth_itervals(filtration, cols_1, sxs_1)
    # ...
    # ints_D, _, _ = nth_itervals(filtration, cols_D-1, sxs_D-1, next=false)
    #
    # [ints_0, ints_1, ..., ints_D]
    ints = [Symbol("ints_", i) for i in 1:D]
    cols = [Symbol("cols_", i) for i in 1:D]
    sxs = [Symbol("sxs_", i) for i in 1:D]

    expr = quote
        ints_0, $(cols[1]), $(sxs[1]) = zeroth_intervals(filtration, ratio)
    end
    for i in 1:D-1
        expr = quote
            $expr
            $(ints[i]), $(cols[i+1]), $(sxs[i+1]) =
                nth_intervals(filtration, $(cols[i]), $(sxs[i]), ratio, Val(C))
        end
    end
    quote
        $expr
        $(ints[D]), _, _ =
            nth_intervals(filtration, $(cols[D]), $(sxs[D]), ratio, Val(C), next=false)

        [ints_0, $(ints...)]
    end
end
