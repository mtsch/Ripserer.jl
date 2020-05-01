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

    PersistenceInterval(birth::T, death::Union{T, Infinity}) where T =
        new{T, Nothing}(birth, death, nothing)
    PersistenceInterval(birth::T, death::Union{T, Infinity}, rep::R) where {T, R} =
        new{T, R}(birth, death, rep)
end

PersistenceInterval(t::Tuple{<:Any, <:Any}) =
    PersistenceInterval(t...)

Base.show(io::IO, int::PersistenceInterval) =
    print(io, "[", int.birth, ", ", int.death, ")")
function Base.show(io::IO, ::MIME"text/plain", int::PersistenceInterval{T}) where T
    print(io, "PersistenceInterval{", T, "}", (int.birth, int.death))
    if !isnothing(int.representative)
        println(io, " with representative:")
        show(io, MIME"text/plain"(), int.representative)
    end
end

"""
    birth(interval::PersistenceInterval)

Get the birth time of `interval`.
"""
birth(int::PersistenceInterval) =
    int.birth

"""
    death(interval::PersistenceInterval)

Get the death time of `interval`. When `T<:AbstractFloat`, `Inf` is returned instead of `∞`.
"""
death(int::PersistenceInterval) =
    int.death
death(int::PersistenceInterval{T}) where T<:AbstractFloat =
    isfinite(int.death) ? int.death : typemax(T)

"""
    death(interval::PersistenceInterval)

Get the persistence of `interval`, which is equal to `death - birth`. When
`T<:AbstractFloat`, `Inf` is returned instead of `∞`.
"""
persistence(int::PersistenceInterval) =
    isfinite(death(int)) ? death(int) - birth(int) : ∞
persistence(int::PersistenceInterval{T}) where T<:AbstractFloat =
    isfinite(death(int)) ? death(int) - birth(int) : typemax(T)

Base.isfinite(int::PersistenceInterval) =
    isfinite(death(int))

"""
    representative(interval::PersistenceInterval)

Get the representative cocycle attached to `interval`. If representatives were not computed,
return `nothing`.
"""
representative(int::PersistenceInterval) =
    int.representative

function Base.iterate(int::PersistenceInterval, i=1)
    if i == 1
        birth(int), i+1
    elseif i == 2
        death(int), i+1
    else
        nothing
    end
end

Base.length(::PersistenceInterval) =
    2

Base.IteratorSize(::Type{<:PersistenceInterval}) =
    Base.HasLength()

Base.IteratorEltype(::Type{<:PersistenceInterval}) =
    Base.HasEltype()

Base.eltype(::Type{<:PersistenceInterval{T}}) where T =
    Union{T, Infinity}

dist_type(::Type{<:PersistenceInterval{T}}) where T =
    T
dist_type(::PersistenceInterval{T}) where T =
    T

function Base.getindex(int, i)
    if i == 1
        birth(int)
    elseif i == 2
        death(int)
    else
        throw(BoundsError(int, i))
    end
end

Base.firstindex(int::PersistenceInterval) =
    1

Base.lastindex(int::PersistenceInterval) =
    2

Base.:(==)(int::PersistenceInterval, tup::NTuple{2}) =
    birth(int) == tup[1] && death(int) == tup[2]

function Base.isless(int1::PersistenceInterval, int2::PersistenceInterval)
    if birth(int1) ≠ birth(int2)
        isless(birth(int1), birth(int2))
    else
        isless(death(int1), death(int2))
    end
end

# diagram ================================================================================ #
"""
    PersistenceDiagram{P<:PersistenceInterval} <: AbstractVector{P}

Type for representing persistence diagrams. Behaves exactly like an array of
`PersistenceInterval`s, but is aware of its dimension and supports pretty printing and
plotting.
"""
struct PersistenceDiagram{P<:PersistenceInterval} <: AbstractVector{P}
    dim       ::Int
    intervals ::Vector{P}
end

PersistenceDiagram(dim, intervals::AbstractVector{<:Tuple}) =
    PersistenceDiagram(dim, PersistenceInterval.(intervals))

Base.show(io::IO, pd::PersistenceDiagram) =
    print(io, length(pd), "-element ", dim(pd), "-dimensional PersistenceDiagram")

function Base.show(io::IO, ::MIME"text/plain", pd::PersistenceDiagram)
    print(io, pd)
    if length(pd) > 0
        print(io, ":")
        limit = get(io, :limit, false) ? first(displaysize(io)) : typemax(Int)
        if length(pd) + 1 < limit
            for i in eachindex(pd)
                if isassigned(pd.intervals, i)
                    print(io, "\n ", pd[i])
                else
                    print(io, "\n #undef")
                end
            end
        else
            for i in 1:limit÷2-2
                if isassigned(pd.intervals, i)
                    print(io, "\n ", pd[i])
                else
                    print(io, "\n #undef")
                end
            end
            print(io, "\n ⋮")
            for i in lastindex(pd)-limit÷2+3:lastindex(pd)
                if isassigned(pd.intervals, i)
                    print(io, "\n ", pd[i])
                else
                    print(io, "\n #undef")
                end
            end
        end
    end
end

Base.size(pd::PersistenceDiagram) =
    size(pd.intervals)

Base.getindex(pd::PersistenceDiagram, i::Integer) =
    pd.intervals[i]

Base.setindex!(pd::PersistenceDiagram, x, i::Integer) =
    pd.intervals[i] = x

Base.firstindex(pd::PersistenceDiagram) =
    1

Base.lastindex(pd::PersistenceDiagram) =
    length(pd.intervals)

Base.similar(pd::PersistenceDiagram) =
    PersistenceDiagram(dim(pd), similar(pd.intervals))

Base.similar(pd::PersistenceDiagram, dims::Tuple) =
    PersistenceDiagram(dim(pd), similar(pd.intervals, dims))

"""
    dim(::PersistenceDiagram)

Get the dimension of persistence diagram.
"""
dim(pd::PersistenceDiagram) =
    pd.dim

dist_type(pd::PersistenceDiagram) =
    dist_type(eltype(pd))
dist_type(::Type{P}) where P<:PersistenceDiagram =
    dist_type(eltype(P))
