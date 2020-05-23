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
