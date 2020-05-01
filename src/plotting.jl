# diarams ================================================================================ #
"""
    t_limits(diagram)

Get the minimum and maximum birth/death times in `diagram` and a boolean specifing if any of
the intervals was infinite. If any of the intervals are infinite, guess a good number to
place `∞` at.
"""
function t_limits(pd::PersistenceDiagram)
    t_min = foldl(min, [birth.(pd); death.(filter(isfinite, pd))], init=0)
    t_max = foldl(max, [birth.(pd); death.(filter(isfinite, pd))], init=0)
    infinite = any(!isfinite, pd)
    if infinite
        rounded = round(t_max, RoundUp)
        t_max = rounded + (t_max ≥ 1 ? ndigits(Int(rounded)) : 0)
    end
    t_min, t_max, infinite
end
function t_limits(pds)
    t_min = 0
    t_max = 0
    infinite = false
    for pd in pds
        t_min = foldl(min, [birth.(pd); death.(filter(isfinite, pd))], init=t_min)
        t_max = foldl(max, [birth.(pd); death.(filter(isfinite, pd))], init=t_max)
        infinite |= any(!isfinite, pd)
    end
    if infinite
        rounded = round(t_max, RoundUp)
        t_max = rounded + (t_max ≥ 1 ? ndigits(Int(rounded)) : 0)
    end
    t_min, t_max, infinite
end

@recipe function f(pd::Union{PersistenceDiagram, AbstractVector{<:PersistenceDiagram}};
                   infinity=nothing)
    t_min, t_max, infinite = t_limits(pd)
    if infinite && !isnothing(infinity)
        t_max = infinity
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
        primary := false

        [t_min, t_max], [t_min, t_max]
    end
    # line at infinity
    if infinite
        @series begin
            seriestype := :hline
            seriescolor := :grey
            line := :dot
            primary := false

            [t_max]
        end
    end

    for pd in pds
        isempty(pd) && continue
        @series begin
            seriestype := :scatter
            label --> "dim $(dim(pd))"

            births = birth.(pd)
            deaths = map(x -> isfinite(x) ? death(x) : t_max, pd)
            births, deaths
        end
    end

    if infinite
        annotations --> (t_max/2, t_max, "∞")
    end
    primary := false
    ()
end

"""
    barcode(diagram; infinity=nothing)

Plot the barcode plot or `AbstractVector` of diagrams. The `infinity` keyword argument
determines where the infinity line is placed. If set to `nothing` the function tries to
guess a good infinity poistion.
"""
barcode(::Union{PersistenceDiagram, AbstractVector{<:PersistenceDiagram}})

@userplot struct Barcode{
    T<:Union{PersistenceDiagram,AbstractVector{<:PersistenceDiagram}},
}
    args::Tuple{T}
end
@recipe function f(bc::Barcode; infinity=nothing)
    arg = only(bc.args)
    if arg isa PersistenceDiagram
        pds = (arg,)
    elseif arg isa Vector{<:PersistenceDiagram}
        pds = arg
    end

    yticks --> []
    xguide --> "t"

    t_min, t_max, infinite = t_limits(pds)
    if infinite && !isnothing(infinity)
        t_max = infinity
    end
    if infinite
        @series begin
            seriestype := :path
            seriescolor := :grey
            line := :dot
            label := "∞"
            primary := false

            [t_max, t_max], [1, sum(length, pds)]
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
                offset += 1
                b, d = birth(int), min(death(int), t_max)
                append!(xs, (b, d, NaN))
                append!(ys, (offset, offset, NaN))
            end
            xs, ys
        end
    end
    if infinite
        annotations --> (t_max, (offset+1)/2, "∞")
    end
    primary := false
    ()
end

# simplex plots ========================================================================== #
const SxVector{D} = AbstractVector{<:AbstractSimplex{D}}

"""
    index_data(indices, args...)

Index into each arg with `indices`, placing `NaN`s where `indices` are zero.
"""
function index_data(indices, pts::AbstractVector{<:NTuple{N}}) where N
    map(indices) do i
        i == 0 ? ntuple(_ -> NaN, N) : Float64.(pts[i])
    end
end

function index_data(indices, x)
    map(indices) do i
        i == 0 ? NaN : Float64(x[i])
    end
end

index_data(indices, x, y) =
    index_data(indices, x), index_data(indices, y)

index_data(indices, x, y, z) =
    index_data(indices, x), index_data(indices, y), index_data(indices, z)

"""
    plottable(sx::AbstractSimplex, args...)
    plottable(sx::PersistenceInterval, args...)
    plottable(sx::Vector{<:AbstractSimplex}, args...)

Turn `sx` into a series that can be plotted. `args...`, should contain the data which tells
the plot where to place simplex vertices.
"""
plottable(sx::AbstractSimplex, args...) =
    plottable([sx], args...)

plottable(int::PersistenceInterval, args...) =
    plottable(representative(int), args...)

plottable(int::PersistenceInterval{<:Any, Nothing}, args...) =
    throw(ArgumentError(
        "interval has no representative. Run `ripserer` with `representatives=true`"
    ))

function plottable(sxs::SxVector{0}, args...)
    indices = only.(vertices.(sxs))
    index_data(indices, args...), [:seriestype => :scatter]
end

function plottable(sxs::SxVector{1}, args...)
    indices = mapreduce(vcat, vertices.(sxs)) do (u, v)
        [u, v, 0]
    end
    index_data(indices, args...), [:seriestype => :path]
end

function plottable(sxs::SxVector, args...)
    indices = mapreduce(vcat, vertices.(sxs)) do vs
        idxs = Int[]
        for (u, v, w) in subsets(vs, Val(3))
            append!(idxs, (u, v, w, u))
        end
        push!(idxs, 0)
        idxs
    end
    index_data(indices, args...), [:seriestype => :path]
end

@recipe function f(
    sx::Union{AbstractSimplex{D}, SxVector{D}, PersistenceInterval}, args...,
) where D
    series, attrs = plottable(sx, args...)
    for (key, value) in attrs
        plotattributes[key] = get(plotattributes, key, value)
    end
    # splat colors over simplices
    for attr in keys(plotattributes)
        value = plotattributes[attr]
        if value isa Vector
            n = D < 2 ? D + 2 : 4 * binomial(D+1, 3) + 1
            plotattributes[attr] = mapreduce(x -> fill(x, n), vcat, value)
        end
    end
    series
end
