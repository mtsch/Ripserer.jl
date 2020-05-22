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

# This type handles plotting the background things in a persistence diagram, such as the
# infinity line and birth = death line.
struct DiagramStuff
    t_lims::NTuple{2, Float64}
    infinite::Bool
    persistence::Bool
end

function clamp_death(int::PersistenceInterval, t_max)
    return isfinite(int) ? death(int) : t_max
end
function clamp_persistence(int::PersistenceInterval, t_max)
    return isfinite(int) ? persistence(int) : t_max
end

@recipe function f(arg::DiagramStuff)
    t_min, t_max = arg.t_lims
    infinite = arg.infinite
    persistence = arg.persistence

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

    # x = y line
    @series begin
        seriestype := :path
        seriescolor := :black
        primary := false

        if persistence
            [t_min, t_max], [t_min, t_min]
        else
            [t_min, t_max], [t_min, t_max]
        end
    end

    ()
end

@recipe function f(
    pd::Union{PersistenceDiagram, AbstractVector{<:PersistenceDiagram}};
    infinity=nothing, persistence=false
)
    t_min, t_max, infinite = t_limits(pd)
    if infinite && !isnothing(infinity)
        t_max = infinity
    end
    if pd isa PersistenceDiagram
        pds = (pd,)
    else
        pds = pd
    end
    different_dims = allunique(dim.(pds))

    xguide --> "birth"
    yguide --> (persistence ? "persistence" : "death")
    legend --> :bottomright
    title --> "Persistence Diagram"

    @series begin
        primary := false
        DiagramStuff((t_min, t_max), infinite, persistence)
    end

    for pd in pds
        isempty(pd) && continue
        @series begin
            seriestype := :scatter
            label --> "H$(dim_str(pd))"
            if different_dims
                markercolor --> dim(pd)+1
                markerstrokecolor --> plotattributes[:markercolor]
            end
            markeralpha --> 0.7

            births = birth.(pd)
            if persistence
                deaths = clamp_persistence.(pd, t_max)
            else
                deaths = clamp_death.(pd, t_max)
            end
            births, deaths
        end
    end
    if infinite
        annotations --> (t_min + (t_max - t_min)/2, t_max, "∞")
    end
    primary := false
    ()
end

@recipe function f(match::Matching; persistence=false, bottleneck=match.bottleneck)
    pds = [match.left, match.right]
    t_min, t_max, infinite = t_limits(pds)

    xguide --> "birth"
    yguide --> (persistence ? "persistence" : "death")
    legend --> :bottomright
    title --> "Persistence Diagram"

    plotattributes[:infinity] = nothing
    plotattributes[:persistence] = persistence
    @series begin
        label --> "matching"

        xs = Float64[]
        ys = Float64[]
        for (l, r) in matching(match, bottleneck=bottleneck)
            append!(xs, (birth(l), birth(r), NaN))
            if persistence
                append!(ys, (clamp_persistence(l, t_max), clamp_persistence(r, t_max), NaN))
            else
                append!(ys, (clamp_death(l, t_max), clamp_death(r, t_max), NaN))
            end
        end
        xs, ys
    end

    @series begin
        pds
    end

    if infinite
        annotations --> (t_min + (t_max - t_min)/2, t_max, "∞")
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
