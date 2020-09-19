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

Base.empty!(chain::WorkingChain) = empty!(chain.heap)
Base.isempty(chain::WorkingChain) = isempty(chain.heap)

Base.eltype(chain::WorkingChain{E}) where E = E

function Base.sizehint!(chain::WorkingChain, size)
    sizehint!(chain.heap, size)
    return chain
end

function Base.pop!(chain::WorkingChain)
    isempty(chain) && return nothing
    heap = chain.heap

    pivot = heappop!(heap, chain.ordering)
    while !isempty(heap)
        if iszero(pivot)
            pivot = heappop!(heap, chain.ordering)
        elseif first(heap) == pivot
            pivot += heappop!(heap, chain.ordering)
        else
            break
        end
    end
    return iszero(pivot) ? nothing : pivot
end

function Base.push!(chain::WorkingChain{E}, element::E) where E
    heappush!(chain.heap, element, chain.ordering)
end

function nonheap_push!(chain::WorkingChain{E}, simplex::AbstractSimplex) where E
    push!(chain.heap, E(simplex))
end
function nonheap_push!(chain::WorkingChain{E}, element::E) where E
    push!(chain.heap, element)
end

repair!(chain::WorkingChain) = heapify!(chain.heap, chain.ordering)

function Base.sort!(chain::WorkingChain; kwargs...)
    return sort!(chain.heap; order=chain.ordering, kwargs...)
end

Base.first(chain::WorkingChain) = first(chain.heap)

function move!(chain::WorkingChain)
    result = eltype(chain)[]
    while (pivot = pop!(chain)) ≢ nothing
        push!(result, pivot)
    end
    return result
end
