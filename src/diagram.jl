# interval =============================================================================== #
"""
    PersistenceInterval{T, C}

The type that represents a persistence interval. It behaves exactly like a
`Tuple{T, Union{T, Infinity}}`, but may have a representative cocycle attached to it.
"""
struct PersistenceInterval{T, R}
    birth          ::T
    death          ::Union{T, Infinity}
    representative ::R

    function PersistenceInterval{T, R}(birth, death, rep::R=nothing) where {T, R}
        if death ≡ ∞
            return new{T, R}(T(birth), death, rep)
        else
            return new{T, R}(T(birth), T(death), rep)
        end
    end
    function PersistenceInterval{T, R}((birth, death), rep::R=nothing) where {T, R}
        if death ≡ ∞
            return new{T, R}(T(birth), death, rep)
        else
            return new{T, R}(T(birth), T(death), rep)
        end
    end
end

function PersistenceInterval(birth::T, death::Union{T, Infinity}, rep::R=nothing) where {T, R}
    return PersistenceInterval{T, R}(birth, death, rep)
end
function PersistenceInterval(t::Tuple{<:Any, <:Any})
    return PersistenceInterval(t...)
end
function Base.convert(::Type{P}, tp::Tuple{<:Any, <:Any}) where P<:PersistenceInterval
    return P(tp)
end

function Base.show(io::IO, int::PersistenceInterval)
    print(io, "[", birth(int), ", ", death(int), ")")
end
function Base.show(io::IO, int::PersistenceInterval{<:AbstractFloat})
    b = round(birth(int), sigdigits=3)
    d = isfinite(death(int)) ? round(death(int), sigdigits=3) : "∞"
    print(io, "[$b, $d)")
end
function Base.show(io::IO, ::MIME"text/plain", int::PersistenceInterval{T}) where T
    print(io, "PersistenceInterval{", T, "}", (birth(int), death(int)))
    if !isnothing(int.representative)
        println(io, " with representative:")
        show(io, MIME"text/plain"(), representative(int))
    end
end

"""
    birth(interval::PersistenceInterval)

Get the birth time of `interval`.
"""
birth(int::PersistenceInterval) = int.birth

"""
    death(interval::PersistenceInterval)

Get the death time of `interval`.
"""
death(int::PersistenceInterval) = int.death

"""
    death(interval::PersistenceInterval)

Get the persistence of `interval`, which is equal to `death - birth`. When
`T<:AbstractFloat`, `Inf` is returned instead of `∞`.
"""
persistence(int::PersistenceInterval) = isfinite(death(int)) ? death(int) - birth(int) : ∞

Base.isfinite(int::PersistenceInterval) = isfinite(death(int))

"""
    representative(interval::PersistenceInterval)

Get the representative cocycle attached to `interval`. If representatives were not computed,
throw an error.
"""
function representative(int::PersistenceInterval)
    if !isnothing(int.representative)
        return int.representative
    else
        error("$int has no representative. Run ripserer with `representatives=true`")
    end
end

function Base.iterate(int::PersistenceInterval, i=1)
    if i == 1
        return birth(int), i+1
    elseif i == 2
        return death(int), i+1
    else
        return nothing
    end
end

Base.length(::PersistenceInterval) = 2
Base.IteratorSize(::Type{<:PersistenceInterval}) = Base.HasLength()
Base.IteratorEltype(::Type{<:PersistenceInterval}) = Base.HasEltype()
Base.eltype(::Type{<:PersistenceInterval{T}}) where T = Union{T, Infinity}

dist_type(::Type{<:PersistenceInterval{T}}) where T = T
dist_type(::PersistenceInterval{T}) where T = T

function Base.getindex(int, i)
    if i == 1
        return birth(int)
    elseif i == 2
        return death(int)
    else
        throw(BoundsError(int, i))
    end
end

Base.firstindex(int::PersistenceInterval) = 1
Base.lastindex(int::PersistenceInterval) = 2

Base.:(==)(int::PersistenceInterval, (b, d)::Tuple) = birth(int) == b && death(int) == d

function Base.isless(int1::PersistenceInterval, int2::PersistenceInterval)
    if birth(int1) ≠ birth(int2)
        return isless(birth(int1), birth(int2))
    else
        return isless(death(int1), death(int2))
    end
end

# diagram ================================================================================ #
"""
    PersistenceDiagram{P<:PersistenceInterval} <: AbstractVector{P}

Type for representing persistence diagrams. Behaves exactly like an array of
`PersistenceInterval`s, but is aware of its dimension and supports pretty printing and
plotting.
"""
struct PersistenceDiagram{T, P<:PersistenceInterval{T}} <: AbstractVector{P}
    dim       ::Int
    threshold ::Union{T, Infinity}
    intervals ::Vector{P}
end

function PersistenceDiagram(dim, threshold, intervals::AbstractVector{<:Tuple})
    return PersistenceDiagram(dim, threshold, PersistenceInterval.(intervals))
end
function PersistenceDiagram(dim, intervals)
    return PersistenceDiagram(dim, ∞, intervals)
end

function show_intervals(io::IO, pd)
    limit = get(io, :limit, false) ? first(displaysize(io)) : typemax(Int)
    if length(pd) + 1 < limit
        for i in eachindex(pd)
            if isassigned(pd, i)
                print(io, "\n ", pd[i])
            else
                print(io, "\n #undef")
            end
        end
    else
        for i in 1:limit÷2-2
            if isassigned(pd, i)
                print(io, "\n ", pd[i])
            else
                print(io, "\n #undef")
            end
        end
        print(io, "\n ⋮")
        for i in lastindex(pd)-limit÷2+3:lastindex(pd)
            if isassigned(pd, i)
                print(io, "\n ", pd[i])
            else
                print(io, "\n #undef")
            end
        end
    end
end
function Base.show(io::IO, pd::PersistenceDiagram)
    print(io, length(pd), "-element ", dim(pd), "-dimensional PersistenceDiagram")
end
function Base.show(io::IO, ::MIME"text/plain", pd::PersistenceDiagram)
    print(io, pd)
    if length(pd) > 0
        print(io, ":")
        show_intervals(io, pd.intervals)
    end
end

threshold(pd::PersistenceDiagram) = pd.threshold

Base.size(pd::PersistenceDiagram) = size(pd.intervals)
Base.getindex(pd::PersistenceDiagram, i::Integer) = pd.intervals[i]
Base.setindex!(pd::PersistenceDiagram, x, i::Integer) = pd.intervals[i] = x
Base.firstindex(pd::PersistenceDiagram) = 1
Base.lastindex(pd::PersistenceDiagram) = length(pd.intervals)

function Base.similar(pd::PersistenceDiagram)
    return PersistenceDiagram(dim(pd), threshold(pd), similar(pd.intervals))
end
function Base.similar(pd::PersistenceDiagram, dims::Tuple)
    return PersistenceDiagram(dim(pd), threshold(pd), similar(pd.intervals, dims))
end

"""
    dim(::PersistenceDiagram)

Get the dimension of persistence diagram.
"""
dim(pd::PersistenceDiagram) = pd.dim

dist_type(pd::PersistenceDiagram) = dist_type(eltype(pd))
dist_type(::Type{P}) where P<:PersistenceDiagram = dist_type(eltype(P))
