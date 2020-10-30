"""
    Chain{F,S,E} <: AbstractVector{ChainElement{S,F}}

A chain should behave exactly like an array of `ChainElements`. Internally, it may store its
elements in a more efficient way.

It can be used as a heap with `heapify!`, `heappush!`, and `heapop!`. When used this way, it
will only return unique elements, summing elements with the same simplex together.

`clean!` can be used to sum duplicates together and sort the chain.
"""
struct Chain{F,S,E} <: AbstractVector{ChainElement{S,F}}
    elements::Vector{E}
end

_internal_eltype(::Type{Mod{2}}, ::Type{S}) where {S} = S
_internal_eltype(::Type{F}, ::Type{S}) where {F,S} = chain_element_type(S, F)
_internal_eltype(chain::Chain) = eltype(chain.elements)

function Chain{F,S}() where {F,S}
    E = _internal_eltype(F, S)
    return Chain{F,S,E}(_internal_eltype(F, S)[])
end
function Chain{F,S}(elements) where {F,S}
    E = _internal_eltype(F, S)
    return Chain{F,S,E}(convert.(E, elements))
end

function Base.summary(io::IO, chain::Chain{F,S}) where {F,S}
    return print(io, length(chain), "-element Chain{$F,$S}")
end

simplex_type(chain::Chain) = simplex_type(eltype(chain))
field_type(chain::Chain) = field_type(eltype(chain))

# Array stuff
Base.size(chain::Chain) = size(chain.elements)
@propagate_inbounds function Base.getindex(chain::Chain{Mod{2}}, i::Int)
    return eltype(chain)(chain.elements[i])
end
@propagate_inbounds function Base.getindex(chain::Chain, i::Int)
    el = chain.elements[i]
    return eltype(chain)(simplex(el), coefficient(el))
end
@propagate_inbounds function Base.getindex(chain::Chain{F,S,E}, is) where {F,S,E}
    return Chain{F,S,E}(chain.elements[is])
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

function Base.copy(chain::Chain{F,S,E}) where {F,S,E}
    return Chain{F,S,E}(copy(chain.elements))
end
function Base.sizehint!(chain::Chain, n)
    sizehint!(chain.elements, n)
    return chain
end

# Heap stuff
function DataStructures.heapify!(chain::Chain, ordering::Base.Ordering)
    heapify!(chain.elements, ordering)
    return chain
end
function DataStructures.heappush!(chain::Chain, el, ordering::Base.Ordering)
    return heappush!(chain.elements, convert(_internal_eltype(chain), el), ordering)
end
function DataStructures.heappop!(chain::Chain, ordering::Base.Ordering)
    if isempty(chain)
        return nothing
    else
        heap = chain.elements
        E = eltype(chain)

        top = convert(E, heappop!(heap, ordering))
        while !isempty(heap)
            if iszero(top)
                top = convert(E, heappop!(heap, ordering))
            elseif convert(E, first(heap)) == top
                top += convert(E, heappop!(heap, ordering))
            else
                break
            end
        end
        return iszero(top) ? nothing : top
    end
end
function heapmove!(chain::Chain{F,S}, ordering::Base.Ordering) where {F,S}
    result = Chain{F,S}()
    while (pivot = heappop!(chain, ordering)) ≢ nothing
        push!(result, pivot)
    end
    return result
end

# Other stuff
function clean!(chain::Chain{F}, ordering::Base.Ordering, factor=one(F)) where {F}
    @inbounds if !isempty(chain)
        sort!(chain.elements; alg=QuickSort, order=ordering)
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
