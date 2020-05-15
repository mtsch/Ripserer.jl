# diagrams =============================================================================== #
"""
    dim_str(pd)

Get `dim` as subscript string.
"""
function dim_str(pd)
    sub_digits = ("₀", "₁", "₂", "₃", "₄", "₅", "₆", "₇", "₈", "₉")
    join(reverse(sub_digits[digits(dim(pd)) .+ 1]))
end

"""
    t_limits(diagram)

Get the minimum and maximum birth/death times in `diagram` and a boolean specifing if any of
the intervals was infinite. If any of the intervals are infinite, use the value specified in
persistence diagram or try to guess a good number to place `∞` at.
"""
function t_limits(pd::PersistenceDiagram)
    t_min = foldl(min, [birth.(pd); death.(filter(isfinite, pd))], init=0)
    t_max = foldl(max, [birth.(pd); death.(filter(isfinite, pd))], init=0)
    infinite = any(!isfinite, pd)
    if infinite
        if threshold(pd) ≢ ∞
            t_max = threshold(pd)
        else
            rounded = round(t_max, RoundUp)
            t_max = rounded + (t_max ≥ 1 ? ndigits(Int(rounded)) : 0)
        end
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
        if all(isa.(pds, PersistenceDiagram))
            infinity = maximum(threshold.(pds))
        else
            infinity = ∞
        end
        if infinity ≢ ∞
            t_max = infinity
        else
            rounded = round(t_max, RoundUp)
            t_max = rounded + (t_max ≥ 1 ? ndigits(Int(rounded)) : 0)
        end
    end
    t_min, t_max, infinite
end

@recipe function f(pd::Union{PersistenceDiagram, AbstractVector{<:PersistenceDiagram}};
                   infinity=nothing, persistence=false)
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
    yguide --> (persistence ? "persistence" : "death")
    legend --> :bottomright
    title --> "Persistence Diagram"

    # x = y line
    @series begin
        seriestype := :path
        seriescolor := :black
        label := ""
        primary := false

        if persistence
            [t_min, t_max], [t_min, t_min]
        else
            [t_min, t_max], [t_min, t_max]
        end
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
            label --> "H$(dim_str(pd))"
            markercolor --> dim(pd)+1
            markerstrokecolor --> dim(pd)+1
            markeralpha --> 0.7
            markerstrokealpha --> 1
            markershape --> :d

            births = birth.(pd)
            if persistence
                deaths = map(x -> isfinite(x) ? death(x) - birth(x) : t_max, pd)
            else
                deaths = map(x -> isfinite(x) ? death(x) : t_max, pd)
            end
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

    t_min, t_max, infinite = t_limits(pds)
    if infinite && !isnothing(infinity)
        t_max = infinity
    end

    yticks --> []
    xguide --> "t"
    legend --> :outertopright
    title --> "Persistence Barcode"

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
            label --> "H$(dim_str(pd))"
            linewidth --> 1
            seriescolor --> dim(pd)+1

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
const IntervalWithRep{D} = PersistenceInterval{
    <:Any, <:AbstractVector{<:Pair{<:AbstractSimplex{D}, <:Any}}
}

"""
    index_data(indices, args...)

Index into each arg with `indices`, placing `NaN`s where `indices` are zero.
"""
function index_data(indices, pts::AbstractVector{<:NTuple{N}}) where N
    map(indices) do i
        i == 0 ? ntuple(_ -> NaN, N) : Float64.(pts[i])
    end
end

function index_data(indices, x::AbstractVector)
    map(indices) do i
        i == 0 ? NaN : Float64(x[i])
    end
end

# images
function index_data(indices, x::AbstractMatrix)
    cart = CartesianIndices(x)
    map(indices) do i
        i == 0 ? (NaN, NaN) : (cart[i][2], cart[i][1])
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
    plottable(simplex.(representative(int)), args...)

plottable(int::PersistenceInterval{<:Any, Nothing}, args...) =
    error("interval has no representative. Run `ripserer` with `representatives=true`")

function plottable(sxs::SxVector{0}, args...)
    indices = only.(vertices.(sxs))
    index_data(indices, args...), [:seriestype => :scatter], 0
end

function plottable(sxs::SxVector{1}, args...)
    indices = mapreduce(vcat, vertices.(sxs)) do (u, v)
        [u, v, 0]
    end
    index_data(indices, args...), [:seriestype => :path], 1
end

function plottable(sxs::SxVector{D}, args...) where D
    indices = mapreduce(vcat, vertices.(sxs)) do vs
        idxs = Int[]
        for (u, v, w) in subsets(vs, Val(3))
            append!(idxs, (u, v, w, u))
        end
        push!(idxs, 0)
        idxs
    end
    index_data(indices, args...), [:seriestype => :path], D
end

apply_threshold(sx::AbstractSimplex, thresh, thresh_strict) =
    diam(sx) ≤ thresh && diam(sx) < thresh_strict ? sx : nothing
apply_threshold(sxs::SxVector, thresh, thresh_strict) =
    filter(sx -> diam(sx) ≤ thresh && diam(sx) < thresh_strict, sxs)
function apply_threshold(int::PersistenceInterval, thresh, thresh_strict)
    reps = filter(sx -> diam(sx) ≤ thresh && diam(sx) < thresh_strict, representative(int))
    PersistenceInterval(birth(int), death(int), reps)
end

@recipe function f(sx::Union{AbstractSimplex, SxVector, PersistenceInterval}, args...;
                   threshold=∞, threshold_strict=∞)
    sx = apply_threshold(sx, threshold, threshold_strict)
    isnothing(sx) && return ()
    series, attrs, D = plottable(sx, args...)
    for (key, value) in attrs
        plotattributes[key] = get(plotattributes, key, value)
    end
    # splat colors and similar attributes over simplices
    for attr in keys(plotattributes)
        value = plotattributes[attr]
        if value isa Vector
            n = D < 2 ? D + 2 : 4 * binomial(D+1, 3) + 1
            plotattributes[attr] = mapreduce(x -> fill(x, n), vcat, value)
        end
    end
    series
end
