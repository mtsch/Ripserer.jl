"""
    WorkingChain{E<:AbstractChainElement, O<:Base.Ordering}

This structure hold the chain that is currently being reduced. Elements are added to the
chain with `push!`. The pivot element is the minimal element with respect to `O` and can be
extracted with `pop!`. `nonheap_push!` followed by `repair!` can also be used to add
elemetns to the chain.

`move!` is used to empty the chain and transfer its elements to a new array.
"""
struct WorkingChain{E<:AbstractChainElement, O<:Base.Ordering}
    heap::Vector{E}
    ordering::O

    function WorkingChain{E}(ordering::O) where {E, O}
        new{E, O}(E[], ordering)
    end
end

Base.empty!(col::WorkingChain) = empty!(col.heap)
Base.isempty(col::WorkingChain) = isempty(col.heap)

function Base.sizehint!(col::WorkingChain, size)
    sizehint!(col.heap, size)
    return col
end

function Base.pop!(column::WorkingChain)
    isempty(column) && return nothing
    heap = column.heap

    pivot = heappop!(heap, column.ordering)
    while !isempty(heap)
        if iszero(pivot)
            pivot = heappop!(heap, column.ordering)
        elseif first(heap) == pivot
            pivot += heappop!(heap, column.ordering)
        else
            break
        end
    end
    return iszero(pivot) ? nothing : pivot
end

function Base.push!(column::WorkingChain{E}, element::E) where E
    heappush!(column.heap, element, column.ordering)
end

function nonheap_push!(column::WorkingChain{E}, simplex::AbstractSimplex) where E
    push!(column.heap, E(simplex))
end
function nonheap_push!(column::WorkingChain{E}, element::E) where E
    push!(column.heap, element)
end

repair!(column::WorkingChain) = heapify!(column.heap, column.ordering)

function Base.sort!(column::WorkingChain; kwargs...)
    return sort!(column.heap; order=column.ordering, kwargs...)
end

Base.first(column::WorkingChain) = first(column.heap)
