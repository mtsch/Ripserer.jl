"""
    PersistenceInterval{T, C}

The type that represents a persistence interval. It behaves exactly like a
`Tuple{T, Union{T, Infinity}}`, but may have a cocycle attached to it.
"""
struct PersistenceInterval{T, C}
    birth   ::T
    death   ::Union{T, Infinity}
    cocycle ::C

    PersistenceInterval(birth::T, death::Union{T, Infinity}) where T =
        new{T, Nothing}(birth, death, nothing)

    PersistenceInterval(birth::T, death::Union{T, Infinity}, cocycle::C) where {T, C} =
        new{T, C}(birth, death, cocycle)
end

Base.show(io::IO, int::PersistenceInterval) =
    print(io, "[", int.birth, ", ", int.death, ")")
function Base.show(io::IO, ::MIME"text/plain", int::PersistenceInterval{T}) where T
    print(io, "PersistenceInterval{", T, "}", (int.birth, int.death))
    if !isnothing(int.cocycle)
        println(io, " with cocycle:")
        show(io, MIME"text/plain"(), int.cocycle)
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

Get the death time of `interval`.
"""
death(int::PersistenceInterval) =
    int.death
"""
    cocycle(interval::PersistenceInterval)

Get the cocycle attached to `interval`. If cocycles were not computed, return `nothing`.
"""
cocycle(int::PersistenceInterval) =
    int.cocycle

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

"""
    PersistenceDiagram{P<:PersistenceInterval} <: AbstractVector{P}

Type for representing persistence diagrams. Behaves exactly like an immutable array of
`PersistenceInterval`s, but is aware of its dimension and supports pretty printing and
plotting.
"""
struct PersistenceDiagram{P<:PersistenceInterval} <: AbstractVector{P}
    dim       ::Int
    intervals ::Vector{P}

    PersistenceDiagram(dim, intervals) =
        new{eltype(intervals)}(dim, sort(intervals))
end

Base.show(io::IO, pd::PersistenceDiagram) =
    print(io, length(pd), "-element ", dim(pd), "-dimensional PersistenceDiagram")

function Base.show(io::IO, ::MIME"text/plain", pd::PersistenceDiagram)
    print(io, pd)
    if length(pd) > 0
        print(io, ":")
        limit = get(io, :limit, false) ? first(displaysize(io)) : typemax(Int)
        if length(pd) + 1 < limit
            for i in pd
                print(io, "\n ", i)
            end
        else
            for i in pd[1:limit÷2-2]
                print(io, "\n ", i)
            end
            print(io, "\n ⋮")
            for i in pd[end-limit÷2+3:end]
                print(io, "\n ", i)
            end
        end
    end
end

Base.size(pd::PersistenceDiagram) =
    size(pd.intervals)

Base.getindex(pd::PersistenceDiagram, i::Integer) =
    pd.intervals[i]
Base.getindex(pd::PersistenceDiagram, is) =
    PersistenceDiagram(dim(pd), pd.intervals[is])

Base.firstindex(pd::PersistenceDiagram) =
    1

Base.lastindex(pd::PersistenceDiagram) =
    length(pd.intervals)

"""
    dim(::PersistenceDiagram)

Get the dimension of persistence diagram.
"""
dim(pd::PersistenceDiagram) =
    pd.dim

# plots ================================================================================== #
"""
    get_max_y_and_inf(diagram::PersistenceDiagram)

Get the time of last death in `diagram`. If it's infinite guess a good number to place `∞`
at.
"""
function get_max_y_and_inf(pd::PersistenceDiagram)
    last_death = max(maximum(death.(pd)),
                     maximum(birth.(pd)))
    last_finite = max(maximum(filter(isfinite, death.(pd))),
                      maximum(birth.(pd)))
    is_finite = last_finite == last_death
    if is_finite
        last_death, false
    else
        rounded = round(last_finite, RoundUp)
        rounded + (last_finite ≥ 1 ? length(digits(Int(rounded))) : 0), true
    end
end
function get_max_y_and_inf(pds)
    last_death = maximum(max(maximum(death.(pd)),
                             maximum(birth.(pd))) for pd in pds)
    last_finite = maximum(max(maximum(filter(isfinite, death.(pd))),
                              maximum(birth.(pd))) for pd in pds)
    is_finite = last_finite == last_death
    if is_finite
        last_death, false
    else
        rounded = round(last_finite, RoundUp)
        rounded + (last_finite ≥ 1 ? length(digits(Int(rounded))) : 0), true
    end
end

@recipe function f(pd::Union{PersistenceDiagram, AbstractVector{<:PersistenceDiagram}};
                   infinity=nothing)
    max_y, has_inf = get_max_y_and_inf(pd)
    if !isnothing(infinity)
        max_y = infinity
    end
    if pd isa PersistenceDiagram
        pds = (pd,)
    else
        pds = pd
    end
    xguide --> "birth"
    yguide --> "death"

    # x = y line
    @series begin
        seriestype := :path
        seriescolor := :black
        label := ""

        [0, max_y], [0, max_y]
    end
    # line at infinity
    if has_inf
        @series begin
            seriestype := :path
            seriescolor := :grey
            line := :dot
            label := "∞"

            [0, max_y], [max_y, max_y]
        end
    end

    for pd in pds
        @series begin
            seriestype := :scatter
            label --> "dim $(dim(pd))"

            births = birth.(pd)
            deaths = map(x -> isfinite(x) ? x : max_y, death.(pd))
            births, deaths
        end
    end
end

@userplot Barcode
@recipe function f(bc::Barcode; infinity=nothing)
    length(bc.args) == 1 ||
        throw(ArgumentError("expected single argument, got $(bc.args)"))
    arg = only(bc.args)
    if arg isa PersistenceDiagram
        pds = [arg]
    elseif !(arg isa Vector{<:PersistenceDiagram})
        throw(ArgumentError("expected `PersistenceDiagram` or " *
                            "`AbstractVector{<:PersistenceDiagram}`, got $(typeof(arg))"))
    else
        pds = arg
    end
    max_y, has_inf = get_max_y_and_inf(pds)
    if !isnothing(infinity)
        max_y = infinity
    end

    yticks --> []
    xguide --> "t"

    if has_inf
        @series begin
            seriestype := :vline
            seriescolor := :grey
            line := :dot
            label := "∞"

            [max_y]
        end
    end

    offset = 0
    for pd in pds
        @series begin
            seriestype := :path
            label --> "dim $(dim(pd))"
            linewidth --> 1

            xs = Float64[]
            ys = Float64[]
            for (i, int) in enumerate(pd)
                b, d = birth(int), min(death(int), max_y)
                append!(xs, (b, d, NaN))
                append!(ys, (offset, offset, NaN))
                offset += 1
            end
            xs, ys
        end
    end
end
