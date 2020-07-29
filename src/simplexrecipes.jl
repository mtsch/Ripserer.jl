"""
    index_data(indices, args...)

Index into each arg with `indices`, placing `NaN`s where `indices` are zero.
"""
function index_data(indices, pts::AbstractVector{<:NTuple{N}}) where N
    return map(indices) do i
        i == 0 ? ntuple(_ -> NaN, N) : Float64.(pts[i])
    end
end
function index_data(indices, x::AbstractVector)
    return map(indices) do i
        i == 0 ? NaN : Float64(x[i])
    end
end
# images
function index_data(indices, x::AbstractMatrix)
    cart = CartesianIndices(x)
    return map(indices) do i
        i == 0 ? (NaN, NaN) : (cart[i][2], cart[i][1])
    end
end
function index_data(indices, x, y)
    return index_data(indices, x), index_data(indices, y)
end
function index_data(indices, x, y, z)
    return index_data(indices, x), index_data(indices, y), index_data(indices, z)
end

"""
    plottable(sx::AbstractSimplex, args...)
    plottable(sx::RepresentativeInterval, args...)
    plottable(sx::Vector{<:AbstractSimplex}, args...)

Turn `sx` into a series that can be plotted. `args...`, should contain the data which tells
the plot where to place simplex vertices.
"""
function plottable(sx::AbstractSimplex, args...)
    return plottable([sx], args...)
end
function plottable(int::RepresentativeInterval, args...)
    return plottable(simplex.(representative(int)), args...)
end
function plottable(rep::AbstractVector{<:AbstractChainElement}, args...)
    return plottable(simplex.(rep), args...)
end
function plottable(element::AbstractChainElement, args...)
    return plottable(simplex(element), args...)
end
function plottable(sxs::AbstractVector{<:AbstractSimplex{0}}, args...)
    indices = sort(only.(vertices.(sxs)))
    return index_data(indices, args...), [:seriestype => :scatter], 0
end
function plottable(sxs::AbstractVector{<:AbstractSimplex{1}}, args...)
    indices = mapreduce(vcat, vertices.(sxs)) do (u, v)
        [u, v, 0]
    end
    return index_data(indices, args...), [:seriestype => :path], 1
end
function plottable(sxs::AbstractVector{<:AbstractSimplex{D}}, args...) where D
    indices = mapreduce(vcat, vertices.(sxs)) do vs
        idxs = Int[]
        for (u, v, w) in subsets(vs, Val(3))
            append!(idxs, (u, v, w, u))
        end
        push!(idxs, 0)
        idxs
    end
    return index_data(indices, args...), [:seriestype => :path], D
end

function apply_threshold(sx::AbstractSimplex, thresh, thresh_strict)
    return birth(sx) ≤ thresh && birth(sx) < thresh_strict ? sx : nothing
end
function apply_threshold(sxs::AbstractVector{<:AbstractSimplex}, thresh, thresh_strict)
    return filter(sx -> birth(sx) ≤ thresh && birth(sx) < thresh_strict, sxs)
end
function apply_threshold(els::AbstractVector{<:AbstractChainElement}, thresh, thresh_strict)
    return apply_threshold(simplex.(els), thresh, thresh_strict)
end
function apply_threshold(::PersistenceDiagrams.AbstractInterval, _, _)
    error("interval has no representative. Run `ripserer` with `representatives=true`")
end
function apply_threshold(int::RepresentativeInterval, thresh, thresh_strict)
    return apply_threshold(representative(int), thresh, thresh_strict)
end

@recipe function f(
    sx::Union{
        AbstractSimplex,
        AbstractChainElement,
        AbstractVector{<:AbstractSimplex},
        AbstractVector{<:AbstractChainElement},
        PersistenceDiagrams.AbstractInterval,
        RepresentativeInterval,
    },
    args...;
    threshold=Inf,
    threshold_strict=Inf
)
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
    return series
end
