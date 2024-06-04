# A chain should behave exactly like an array of `ChainElements`. Internally, it may store
# its elements in a more efficient way.

# It can be used as a heap with `heapify!`, `heappush!`, and `heapop!`. When used this way,
# it will only return unique elements, summing elements with the same simplex together.

# `clean!` can be used to sum duplicates together and sort the chain.
"""
    Chain{F,S,E} <: AbstractVector{ChainElement{S,F}}

An internal representation of a chain. Behaves like an array of pairs `S => F`, where `S` is
the simplex type, and `F` is the coefficient type (probably a subtype of [`Mod`](@ref)).

Most functions that can be called on `AbstractCell`s can also be called on the elements of a
`Chain`.

# Examples

```jldoctest; setup=(using Random; Random.seed!(1337))
julia> data = [(rand(), rand(), rand()) for _ in 1:100];

julia> result = ripserer(data; reps=true, modulus=7);

julia> chain = result[end][end].representative
28-element Chain{Mod{7},Simplex{1, Float64, Int64}}:
 +Simplex{1}((68, 54), 0.25178097927369875) => 6 mod 7
 +Simplex{1}((54, 46), 0.2575262844682746) => 1 mod 7
 +Simplex{1}((88, 56), 0.25936896586973557) => 1 mod 7
 +Simplex{1}((79, 22), 0.26690101517910964) => 6 mod 7
 +Simplex{1}((53, 13), 0.27161939814693414) => 6 mod 7
 +Simplex{1}((24, 13), 0.28686687771884056) => 6 mod 7
 +Simplex{1}((56, 54), 0.2869400047523196) => 6 mod 7
 +Simplex{1}((50, 42), 0.29274772369750884) => 6 mod 7
 +Simplex{1}((93, 50), 0.3022995886824861) => 1 mod 7
 +Simplex{1}((73, 54), 0.30771273548368216) => 6 mod 7
 ⋮
 +Simplex{1}((70, 34), 0.3412250751087395) => 6 mod 7
 +Simplex{1}((49, 42), 0.3414959312908324) => 6 mod 7
 +Simplex{1}((49, 7), 0.35568373080560794) => 6 mod 7
 +Simplex{1}((88, 46), 0.35629002982624475) => 1 mod 7
 +Simplex{1}((95, 88), 0.3566986102316938) => 6 mod 7
 +Simplex{1}((95, 54), 0.363029725360828) => 6 mod 7
 +Simplex{1}((50, 13), 0.3648936036309352) => 6 mod 7
 +Simplex{1}((24, 3), 0.3653704317993549) => 6 mod 7
 +Simplex{1}((85, 42), 0.369953146904745) => 6 mod 7

julia> simplex.(chain)
28-element Vector{Simplex{1, Float64, Int64}}:
 +Simplex{1}((68, 54), 0.25178097927369875)
 +Simplex{1}((54, 46), 0.2575262844682746)
 +Simplex{1}((88, 56), 0.25936896586973557)
 +Simplex{1}((79, 22), 0.26690101517910964)
 +Simplex{1}((53, 13), 0.27161939814693414)
 +Simplex{1}((24, 13), 0.28686687771884056)
 +Simplex{1}((56, 54), 0.2869400047523196)
 +Simplex{1}((50, 42), 0.29274772369750884)
 +Simplex{1}((93, 50), 0.3022995886824861)
 +Simplex{1}((73, 54), 0.30771273548368216)
 ⋮
 +Simplex{1}((70, 34), 0.3412250751087395)
 +Simplex{1}((49, 42), 0.3414959312908324)
 +Simplex{1}((49, 7), 0.35568373080560794)
 +Simplex{1}((88, 46), 0.35629002982624475)
 +Simplex{1}((95, 88), 0.3566986102316938)
 +Simplex{1}((95, 54), 0.363029725360828)
 +Simplex{1}((50, 13), 0.3648936036309352)
 +Simplex{1}((24, 3), 0.3653704317993549)
 +Simplex{1}((85, 42), 0.369953146904745)

julia> vertices.(chain)
28-element Vector{Tuple{Int64, Int64}}:
 (68, 54)
 (54, 46)
 (88, 56)
 (79, 22)
 (53, 13)
 (24, 13)
 (56, 54)
 (50, 42)
 (93, 50)
 (73, 54)
 ⋮
 (70, 34)
 (49, 42)
 (49, 7)
 (88, 46)
 (95, 88)
 (95, 54)
 (50, 13)
 (24, 3)
 (85, 42)

julia> coefficient.(chain)
28-element Vector{Mod{7}}:
 6 mod 7
 1 mod 7
 1 mod 7
 6 mod 7
 6 mod 7
 6 mod 7
 6 mod 7
 6 mod 7
 1 mod 7
 6 mod 7
       ⋮
 6 mod 7
 6 mod 7
 6 mod 7
 1 mod 7
 6 mod 7
 6 mod 7
 6 mod 7
 6 mod 7
 6 mod 7

```
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
function Chain(elements::AbstractVector{ChainElement{S,F}}) where {S,F}
    return Chain{F,S}(collect(elements))
end

function Base.summary(io::IO, chain::Chain{F,S}) where {F,S}
    return print(io, length(chain), "-element Chain{$F,$S}")
end

simplex_type(chain::Chain) = simplex_type(eltype(chain))
coefficient_type(chain::Chain) = coefficient_type(eltype(chain))

# Array overloads
function Base.size(chain::Chain)
    return size(chain.elements)
end
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
function Base.setindex!(chain::Chain{F,S}, (s, c)::Tuple{S,F}, i::Int) where {F,S}
    element = ChainElement{S,F}(s, c)
    return chain[i] = element
end
function Base.IndexStyle(::Chain)
    return Base.IndexLinear()
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
# Compat: these methods are here for compatibility reasons.
# TODO: Remove when 1.6 becomes LTS
function Base.append!(chain::Chain, elements)
    n = length(chain)
    resize!(chain, n + length(elements))
    @inbounds for (i, e) in enumerate(elements)
        chain[i + n] = e
    end
    return chain
end
function Base.append!(chain1::C, chain2::C) where {C<:Chain}
    return append!(chain1.elements, chain2.elements)
end
function Base.push!(chain::Chain, element)
    resize!(chain, length(chain) + 1)
    return @inbounds chain[end] = element
end

# Heap stuff
function DataStructures.heapify!(chain::Chain, ordering::Base.Ordering)
    heapify!(chain.elements, ordering)
    return chain
end
function DataStructures.heappush!(chain::Chain, val, ordering::Base.Ordering)
    return heappush!(chain.elements, convert(_internal_eltype(chain), val), ordering)
end
function DataStructures.heappush!(
    chain::Chain{F,S}, (s, c)::Tuple{S,F}, ordering::Base.Ordering
) where {F,S}
    element = ChainElement{S,F}(s, c)
    return heappush!(chain, element, ordering)
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
function heapmove!(dst::Chain{F,S}, chain::Chain{F,S}, ordering::Base.Ordering) where {F,S}
    while (pivot = heappop!(chain, ordering)) ≢ nothing
        push!(dst, pivot)
    end
    return dst
end
function heapmove!(chain::Chain{F,S}, ordering::Base.Ordering) where {F,S}
    return heapmove!(Chain{F,S}(), chain, ordering)
end

_no_ordering_error() = error("no ordering given")
DataStructures.heapify!(::Chain) = _no_ordering_error()
DataStructures.heappop!(::Chain) = _no_ordering_error()
DataStructures.heappush!(::Chain, _) = _no_ordering_error()
heapmove!(::Chain) = _no_ordering_error()
heapmove!(::Chain, ::Chain) = _no_ordering_error()

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

"""
    is_cocycle(filtration, chain::Chain, time)

Check whether a `chain` is a cocycle in `filtration` at specified `time`.
"""
function is_cocycle(filtration, chain::Chain{F,S}, time) where {F,S}
    buffer = Chain{F,simplex_type(filtration, dim(S) + 1)}()
    for (simplex, coefficient) in chain
        for cofacet in coboundary(filtration, simplex)
            if birth(cofacet) < time
                heappush!(buffer, (cofacet, coefficient), Base.Order.Forward)
            end
        end
    end
    return isnothing(heappop!(buffer, Base.Order.Forward))
end
