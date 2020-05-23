"""
    ReductionMatrix{C, SE}

A representation of the reduction matrix. It is indexed by simplices of type `C` (or chain
element types) and holds chain elements of one dimension lower. Supports the following
operations.

* `insert_column!(::ReductionMatrix, i::K)`: add a new column with column index `i`.
* `has_column(::ReductionMatrix, i::K)`: return `true` if the matrix has any entries in the
  `i`-th column.
* `push!(::ReductionMatrix, val::V)`: push `val` to *the last* column that was added to the
  matrix.
"""
struct ReductionMatrix{C<:AbstractSimplex, SE<:AbstractChainElement}
    column_index ::Dict{C, Int}
    colptr       ::Vector{Int}
    nzval        ::Vector{SE}

    function ReductionMatrix{C, SE}() where {C, SE}
        new{C, SE}(Dict{C, Int}(), Int[1], SE[])
    end
end

has_column(rm::ReductionMatrix{C}, sx::C) where C =
    haskey(rm.column_index, abs(sx))
has_column(rm::ReductionMatrix{C}, ce::AbstractChainElement{C}) where C =
    haskey(rm.column_index, simplex(ce))

insert_column!(rm::ReductionMatrix{C}, ce::AbstractChainElement{C}) where C =
    insert_column!(rm, simplex(ce))
function insert_column!(rm::ReductionMatrix{C}, sx::C) where C
    rm.column_index[abs(sx)] = length(rm.colptr)
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

Base.eltype(::Type{<:ReductionMatrix{<:Any, SE}}) where SE =
    SE
Base.length(rm::ReductionMatrix) =
    length(rm.colptr) - 1
Base.lastindex(rm::ReductionMatrix) =
    length(rm.colptr) - 1
Base.getindex(rm::ReductionMatrix{C}, ce::AbstractChainElement{C}) where C =
    RMColumnIterator{typeof(rm)}(rm, rm.column_index[simplex(ce)])
Base.getindex(rm::ReductionMatrix{C}, sx::C) where C =
    RMColumnIterator{typeof(rm)}(rm, rm.column_index[sx])

"""
    RMColumnIterator{R}

An iterator over a column of a `ReductionMatrix`.
"""
struct RMColumnIterator{R<:ReductionMatrix}
    rm  ::R
    idx ::Int
end

Base.IteratorSize(::Type{RMColumnIterator}) =
    Base.HasLength()
Base.IteratorEltype(::Type{RMColumnIterator}) =
    Base.HasEltype()
Base.eltype(::Type{<:RMColumnIterator{R}}) where R =
    eltype(R)
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
    Column{CE<:AbstractChainElement}

Wrapper around `BinaryMinHeap{CE}`. Acts like a heap of chain elements, where elements with
the same simplex are summed together and elements with coefficient value `0` are ignored.
Support `push!`ing chain elements or simplices. Unlike a regular heap, it returns nothing
when `pop!` is used on an empty column.
"""
struct Column{CE<:AbstractChainElement}
    heap::BinaryMinHeap{CE}

    Column{CE}() where CE =
        new{CE}(BinaryMinHeap{CE}())
end

Base.empty!(col::Column) =
    empty!(col.heap.valtree)
Base.isempty(col::Column) =
    isempty(col.heap)

function Base.sizehint!(col::Column, size)
    sizehint!(col.heap, size)
    col
end

"""
    move_mul!(dst, col::Column{CE}, times=one(CE))

Move contents of column into `dst` by repeatedly calling `pop_pivot!` on the column and
`push!`ing the pivot to `dst`. Multipy all elements that are moved are multiplied by
`times`.
"""
function move_mul!(dst, col::Column{CE}, times=one(CE)) where CE
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
`nothing` if column is empty or if all pivots summed to 0.
"""
function pop_pivot!(column::Column)
    isempty(column) && return nothing
    heap = column.heap

    pivot = pop!(heap)
    while !isempty(heap)
        if iszero(pivot)
            pivot = pop!(heap)
        elseif top(heap) == pivot
            pivot += pop!(heap)
        else
            break
        end
    end
    iszero(pivot) ? nothing : pivot
end

"""
    pivot(column)

Return the pivot of the column - the element with the lowest diameter.
"""
function pivot(column::Column)
    heap = column.heap
    pivot = pop_pivot!(column)
    if !isnothing(pivot)
        push!(column, pivot)
    end
    pivot
end

Base.push!(column::Column{CE}, simplex) where CE =
    push!(column, CE(simplex))
function Base.push!(column::Column{CE}, element::CE) where CE
    heap = column.heap
    if !isempty(heap) && top(heap) == element
        heap.valtree[1] += element
    else
        push!(heap, element)
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

Base.length(s::DisjointSetsWithBirth) =
    length(s.parents)

function DataStructures.find_root!(s::DisjointSetsWithBirth, x)
    parents = s.parents
    p = parents[x]
    @inbounds if parents[p] != p
        parents[x] = p = find_root!(s, p)
    end
    p
end

"""
    find_leaves!(s::DisjointSetsWithBirth, x)

Find all leaves below `x`, i.e. vertices that have `x` as root.
"""
function find_leaves!(s::DisjointSetsWithBirth, x)
    leaves = Int[]
    for i in 1:length(s)
        find_root!(s, i) == x && push!(leaves, i)
    end
    leaves
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
