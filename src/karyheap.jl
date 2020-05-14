# TODO: propose this to DataStructures once
# https://github.com/JuliaCollections/DataStructures.jl/pull/547 is merged.
import Base.Order: Forward, Ordering, lt

# K-Heap operations on flat arrays
# --------------------------------

function heapchildren(i::Integer, ::Val{K}) where K
    K * (i-1) + 2:K * (i-1) + K + 1
end
firstchild(i::Integer, ::Val{K}) where K = first(heapchildren(i, Val(K)))
heapparent(i::Integer, ::Val{K}) where K = ifelse(i == 1, 0, div(i-2, K) + 1)

function percolate_down!(xs::AbstractArray, i::Integer, ::Val{K},
                         x=xs[i], o::Ordering=Forward, len::Integer=length(xs)) where K
    @inbounds while (j = firstchild(i, Val(K))) <= len
        for c in heapchildren(i, Val(K))
            c > len && break
            j = ifelse(lt(o, xs[j], xs[c]), j, c)
        end
        if lt(o, xs[j], x)
            xs[i] = xs[j]
            i = j
        else
            break
        end
    end
    xs[i] = x
end

# Binary min-heap percolate down.
percolate_down!(xs::AbstractArray, i::Integer, ::Val{K}, o::Ordering, len::Integer=length(xs)) where K =
    percolate_down!(xs, i, Val(K), xs[i], o, len)

# Binary min-heap percolate up.
function percolate_up!(xs::AbstractArray, i::Integer, ::Val{K}, x=xs[i], o::Ordering=Forward) where K
    @inbounds while (j = heapparent(i, Val(K))) >= 1
        if lt(o, x, xs[j])
            xs[i] = xs[j]
            i = j
        else
            break
        end
    end
    xs[i] = x
end

percolate_up!(xs::AbstractArray, i::Integer, ::Val{K}, o::Ordering) where K =
    percolate_up!(xs, i, Val(K), xs[i], o)

# Turn an arbitrary array into a K-ary min-heap (by default) in linear time.
"""
    heapify!(v, ord::Ordering=Forward)
In-place [`heapify`](@ref).
"""
function DataStructures.heapify!(xs::AbstractArray, ::Val{K}, o::Ordering=Forward) where K
    for i in heapparent(length(xs), Val(K)):-1:1
        percolate_down!(xs, i, o, Val(K))
    end
    xs
end

"""
    heappop!(v, [ord])
Given a binary heap-ordered array, remove and return the lowest ordered element.
For efficiency, this function does not check that the array is indeed heap-ordered.
"""
function DataStructures.heappop!(xs::AbstractArray, ::Val{K}, o::Ordering=Forward) where K
    x = xs[1]
    y = pop!(xs)
    if !isempty(xs)
        percolate_down!(xs, 1, Val(K), y, o)
    end
    x
end

"""
    heappush!(v, x, [ord])
Given a binary heap-ordered array, push a new element `x`, preserving the heap property.
For efficiency, this function does not check that the array is indeed heap-ordered.
"""
function DataStructures.heappush!(xs::AbstractArray, x, ::Val{K}, o::Ordering=Forward) where K
    push!(xs, x)
    percolate_up!(xs, length(xs), Val(K), o)
    xs
end

struct KAryHeap{K, T, O<:Base.Ordering} <: AbstractHeap{T}
    ordering::O
    valtree::Vector{T}

    function KAryHeap{K, T, O}() where {K, T, O}
        new{K, T, O}(O(), T[])
    end

    function KAryHeap{K, T, O}(xs::AbstractVector{T}) where {K, T, O}
        ordering = O()
        valtree = heapify(xs, ordering, Val(K))
        new{K, T, O}(ordering, valtree)
    end
end

const KAryMinHeap{K, T} = KAryHeap{K, T, Base.Order.ForwardOrdering}
const KAryMaxHeap{K, T} = KAryHeap{K, T, Base.Order.ReverseOrdering}

#################################################
#
#   interfaces
#
#################################################

"""
    length(h::KAryHeap)

Returns the number of elements in heap `h`.
"""
Base.length(h::KAryHeap) = length(h.valtree)

"""
    isempty(h::KAryHeap)

Returns whether the heap `h` is empty.
"""
Base.isempty(h::KAryHeap) = isempty(h.valtree)

"""
    push!(h::KAryHeap, value)

Adds the `value` element to the heap `h`.
"""
function Base.push!(h::KAryHeap{K}, v) where K
    heappush!(h.valtree, v, Val(K), h.ordering)
    h
end

"""
    top(h::KAryHeap)

Returns the element at the top of the heap `h`.
"""
@inline DataStructures.top(h::KAryHeap) = h.valtree[1]

"""
    pop!(h::KAryHeap)

Removes and returns the element at the top of the heap `h`.
"""
Base.pop!(h::KAryHeap{K}) where K = heappop!(h.valtree, Val(K), h.ordering)

"""
    sizehint!(h::KAryHeap, n::Integer)

Suggest that heap `h` reserve capacity for at least `n` elements. This can improve
performance.
"""
function Base.sizehint!(h::KAryHeap, n::Integer)
    sizehint!(h.valtree, n)
    h
end
