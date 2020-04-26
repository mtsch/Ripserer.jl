struct PersistenceInterval{T, C}
    birth   ::T
    death   ::Union{T, Infinity}
    cocycle ::C

    PersistenceInterval(birth::T, death::Union{T, Infinity}) where T =
        new{T, Nothing}(birth, death, nothing)

    PersistenceInterval(birth::T, death::Union{T, Infinity}, cocycle::C) where {T, C} =
        new{T, C}(birth, death, cocycle)
end

birth(int::PersistenceInterval) =
    int.birth
death(int::PersistenceInterval) =
    int.death

# Make PersistenceInterval indistinguishable from tuple.
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
function Base.getindex(int, i)
    if i == 1
        birth(int)
    elseif i == 2
        death(int)
    else
        throw(BoundsError(int, i))
    end
end
Base.lastindex(int::PersistenceInterval) =
    2
Base.:(==)(int::PersistenceInterval, tup::NTuple{2}) =
    birth(int) == tup[1] && death(int) == tup[2]
