export Chain

abstract type AbstractChain{F,S,O<:Base.Ordering} <: AbstractVector{ChainElement{S,F}} end

struct Chain{F,S,E,O} <: AbstractChain{F,S,O}
    elements::Vector{E}
    ordering::O
end

_internal_eltype(::Type{Mod{2}}, ::Type{S}) where {S} = S
_internal_eltype(::Type{F}, ::Type{S}) where {F,S} = chain_element_type(S, F)
_internal_eltype(chain::Chain) = eltype(chain.elements)

function Chain{F,S}(ordering::O=Base.Order.Forward) where {F,S,O<:Base.Ordering}
    E = _internal_eltype(F,S)
    return Chain{F,S,E,O}(_internal_eltype(F,S)[], ordering)
end
function Chain{F,S}(elements, ordering::O=Base.Order.Forward) where {F,S,O<:Base.Ordering}
    E = _internal_eltype(F,S)
    return Chain{F,S,E,O}(convert.(E, elements), ordering)
end

# Array stuff
Base.size(chain::Chain) = size(chain.elements)
@propagate_inbounds function Base.getindex(chain::Chain{Mod{2}}, i::Int)
    return eltype(chain)(chain.elements[i])
end
@propagate_inbounds function Base.getindex(chain::Chain, i::Int)
    el = chain.elements[i]
    return eltype(chain)(simplex(el), coefficient(el))
end
@propagate_inbounds function Base.getindex(chain::Chain{F,S,E,O}, is) where {F,S,E,O}
    return Chain{F,S,E,O}(chain.elements[is], chain.ordering)
end

function Base.setindex!(chain::Chain, val, i::Int)
    return chain.elements[i] = convert(_internal_eltype(chain), val)
end
function Base.setindex!(chain::Chain{Mod{2}}, val, i::Int)
    return chain.elements[i] = convert(_internal_eltype(chain), val)
end

function Base.resize!(chain::Chain, n)
    resize!(chain.elements, n)
    return chain
end
function Base.empty!(chain::Chain)
    empty!(chain.elements)
    return chain
end
function Base.pop!(chain::Chain{F,S}) where {F,S}
    return convert(ChainElement{S,F}, pop!(chain.elements))
end

function Base.copy(chain::Chain{F,S,E,O}) where {F,S,E,O}
    return Chain{F,S,E,O}(copy(chain.elements), chain.ordering)
end

# Heap stuff
function DataStructures.heapify!(chain::Chain)
    heapify!(chain.elements, chain.ordering)
    return chain
end
function DataStructures.heappush!(chain::Chain, el)
    return heappush!(chain.elements, convert(_internal_eltype(chain), el), chain.ordering)
end
function DataStructures.heappop!(chain::Chain)
    if isempty(chain)
        return nothing
    else
        heap = chain.elements
        E = eltype(chain)

        top = convert(E, heappop!(heap, chain.ordering))
        while !isempty(heap)
            if iszero(top)
                top = convert(E, heappop!(heap, chain.ordering))
            elseif convert(E, first(heap)) == top
                top += convert(E, heappop!(heap, chain.ordering))
            else
                break
            end
        end
        return iszero(top) ? nothing : top
    end
end
function heapmove!(chain::Chain{F,S}) where {F,S}
    result = Chain{F,S}(chain.ordering)
    while (pivot = heappop!(chain)) ≢ nothing
        push!(result, pivot)
    end
    return result
end

# Other stuff
function clean!(chain::Chain{F}, factor=one(F)) where F
    @inbounds if !isempty(chain)
        sort!(chain.elements; alg=QuickSort, order=chain.ordering)
        i = 0
        current = chain[1]
        for j in 2:length(chain)
            top = chain[j]
            if top == current
                current += top
            else
                if !iszero(current)
                    i += 1
                    chain[i] = current * factor
                end
                current = top
            end
        end
        if !iszero(current)
            i += 1
            chain[i] = current * factor
        end
        resize!(chain, i)
    end
    return chain
end
