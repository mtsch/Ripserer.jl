"""
    module Partition

This submodule contains the following partition functions to be used with
[`CircularCoordinates`](@ref):

* `Partition.linear(r, d) = max(r - d, 0.0)`
* `Partition.quadratic(r, d) = (max(r - d, 0.0))^2`
* `Partition.exponential(r, d) = r - d > 0 ? exp(r^2/(d^2 - r^2)) : 0.0`

"""
module Partition
"""
    Partition.linear(r, d) = max(r - d, 0.0)
"""
linear(r, d) = max(r - d, 0.0)
"""
    Partition.quadratic(r, d) = (max(r - d, 0.0))^2
"""
quadratic(r, d) = (max(r - d, 0.0))^2
"""
    Partition.exponential(r, d) = r - d > 0 ? exp(r^2/(d^2 - r^2)) : 0.0
"""
exponential(r, d) = r - d > 0 ? exp(r^2 / (d^2 - r^2)) : 0.0
end

###
### Landmark selection
###
function _maxmin_sample(points, n, metric=Euclidean())
    l = rand(1:length(points))
    landmarks = zeros(Int, n)
    i = 1
    closest_landmark = fill(Inf, length(points))
    while i ≤ n
        landmarks[i] = l
        for (j, p) in enumerate(points)
            dist = metric(SVector(points[l]), SVector(p))
            closest_landmark[j] = min(closest_landmark[j], dist)
        end
        i += 1
        l = argmax(closest_landmark)
    end
    return SVector.(points[landmarks]), maximum(closest_landmark)
end

function _landmarks_and_radius(points, landmarks::Int, metric)
    return _maxmin_sample(points, landmarks, metric)
end
function _landmarks_and_radius(points, landmarks::AbstractVector{Int}, metric)
    return _landmarks_and_radius(points, points[landmarks], metric)
end
function _landmarks_and_radius(points, landmarks, metric)
    radius = 0.0
    for p in points
        nearest = minimum(metric(SVector(l), SVector(p)) for l in landmarks)
        radius = max(nearest, radius)
    end
    return SVector.(landmarks), radius
end

###
### Cocycle massage
###
function _to_integer_coefficients(chain::Chain{Mod{M},S}) where {M,S}
    result = Chain{Int,S}()
    resize!(result, length(chain))
    for (i, (sx, coef)) in enumerate(chain)
        new_coef = Int(coef) > (M - 1) ÷ 2 ? Int(coef) - M : Int(coef)
        @inbounds result[i] = (sx, Int(new_coef))
    end
    return result
end

function _zero_coboundary_matrix(filtration::AbstractFiltration, time)
    is = Int[]
    js = Int[]
    vs = Float64[]
    for e in edges(filtration)
        if birth(e) < time
            i = index(e)
            u, v = e
            append!(is, (i, i))
            append!(js, (u, v))
            append!(vs, (-1.0, 1.0))
        end
    end
    n = nv(filtration)
    return sparse(is, js, vs, n * (n - 1) ÷ 2, n)
end

function _to_vector(filtration::AbstractFiltration, cocycle::Chain{Int}, time)
    n = nv(filtration)
    vector = zeros(n * (n - 1) ÷ 2)
    for (sx, c) in cocycle
        if birth(sx) < time
            vector[index(sx)] = Float64(c)
        end
    end
    return vector
end

function _mod_z(t)
    t = rem(t, 1.0)
    return t < 0.0 ? 1.0 + t : t
end

function _harmonic_smoothing(filtration, chain; time)
    integral_cocycle = _to_integer_coefficients(chain)
    if !is_cocycle(filtration, integral_cocycle, time)
        error(
            "the cocycle cannot be converted to `Int` coefficients. ",
            "Try using changing the landmarks, or using a different `modulus`.",
        )
    end
    coboundary = _zero_coboundary_matrix(filtration, time)
    real_cocycle = _to_vector(filtration, integral_cocycle, time)

    minimizer = coboundary \ real_cocycle
    real_cocycle .-= coboundary * minimizer
    return minimizer, real_cocycle
end

const PRIMES = [7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79]

struct CircularCoordinateData
    radius::Float64
    coordinate::Vector{Float64}
    cocycle::Vector{Float64}
end

"""
    CircularCoordinates

This `struct` implements Perea's sparse circular coordinates (see reference below).

The idea behind this method is that you pick a subset of your data set (also called
landmarks), compute persistent cohomology on that subset, and use the result to construct
circular coordinates on the whole data set.

For this to work correctly, you need the following.

* A set of landmarks that cover the points well. The default minmax sampling method is
usually a good choice here.
* The persistence diagram on the landmarks must have an interval that is persistent enough.

The circular coordinates returned are on the interval [0, 1). If you are looking for angles,
make sure to multiply them by `2π`.

If you try to compute circular coordinates for a point that is not near any landmark, the
result will be `missing`.

# Constructor

`CircularCoordinates([::AbstractFiltration, ]points, landmarks; kwargs...)`

# Arguments

* `points`: a vector of points. The `eltype` of this vector can be `Tuple`s, `SVector`s, or
  similar.

* `landmarks`: can be an integer, a vector of indices, or a vector of points. If set to an
  integer, it sets the number of landmarks to choose with maxmin sampling. If you are
  looking for non-sparse circular coordinates, use `landmarks=eachindex(points)`.

* `out_dim`: number of most persistent persistence intervals to use to compute coordinates.
  If less than `out_dim` suitable intervals are found, the construction will show a warning
  and construct `CircularCoordinates` with a lower `out_dim`. This warning can be suppressed
  by passing `warn=false`.

* `partition`: a function that defines the partition of unity used when determining the
  coordinates. Should take two arguments, `r`, the radius of the balls, and `d`, the
  distance from the landmark. The partition function should only have support on the ball
  around the landmark (in other words, it should evaluate to 0 when d ≥ r). If that is not
  the case, the circular coordinates are not well-defined and the results may not always
  make sense. See [`Partition`](@ref) for a small collection of predefined functions.

* `metric`: a metric from [Distances.jl](https://github.com/JuliaStats/Distances.jl). Can
  only be set if the filtration argument is `Rips`.

* `kwargs...`: additional keyword arguments passed to `ripserer`. Note that `modulus` is set
  to a random prime between 7 and 79 by default.

# Example

```jldoctest
julia> data = [(sin(t), cos(t)) for t in range(0, 2π, length=101)[1:end-1]];

julia> cc = CircularCoordinates(data, 1:10:100)
CircularCoordinates(
  out_dim=1,
  radius=(0.9510565162951535,),
  n_landmarks=10,
  partition=linear,
  metric=Distances.Euclidean(0.0),
)

julia> summary(cc(data))
"100×1 Array{Union{Missing, Float64},2}"

julia> summary(cc(data, 1))
"100-element Array{Union{Missing, Float64},1}"

```

# Reference

Perea, J. A. (2020, June). Sparse circular coordinates via principal Z-bundles. In
Topological Data Analysis: The Abel Symposium 2018 (Vol. 15, p. 435). Springer Nature.

"""
struct CircularCoordinates{F,P<:SVector,M}
    out_dim::Int
    landmarks::Vector{P}
    partition::F
    metric::M

    coordinate_data::Vector{CircularCoordinateData}
    meta::NamedTuple

    partition_buffer::Vector{Float64}
end

function CircularCoordinates(args...; kwargs...)
    return CircularCoordinates(Rips, args...; kwargs...)
end
function CircularCoordinates(
    ::Type{F},
    points::AbstractVector,
    landmarks;
    out_dim=1,
    partition=Partition.linear,
    metric=Euclidean(),
    coverage=1,
    modulus=rand(PRIMES),
    threshold=nothing,
    verbose=false,
    warn=true,
    kwargs...,
) where {F<:AbstractFiltration}
    @prog_print verbose "Determining radius...   "
    landmarks, min_radius = _landmarks_and_radius(points, landmarks, metric)
    @prog_println verbose "done."

    # TODO with new ripserer interface, change this.
    flt_kwargs = metric == Euclidean() ? NamedTuple() : (metric = Euclidean())
    flt_kwargs = isnothing(threshold) ? flt_kwargs : (; threshold=threshold, flt_kwargs...)

    # compute cohomology
    @prog_print verbose "Computing cohomology... "
    filtration = F(landmarks; flt_kwargs...)
    diagram = ripserer(filtration; modulus=modulus, reps=true, kwargs...)[2]
    @prog_println verbose "done."

    @prog_print verbose "Smoothing cocycles...   "
    coord_data = CircularCoordinateData[]
    for d in 1:out_dim
        if length(diagram) < d
            break
        end
        interval = diagram[end + 1 - d]
        b, d = interval
        if !isfinite(death(interval))
            d = threshold(filtration)
        end
        if max(b, min_radius) ≥ d / 2
            break
        else
            radius = (1 - coverage) * max(b, min_radius) + coverage * d / 2
        end
        coords, cocycle = _harmonic_smoothing(
            filtration, interval.representative; time=2 * radius
        )
        push!(coord_data, CircularCoordinateData(radius, coords, cocycle))
    end
    @prog_println verbose "done."
    n_success = length(coord_data)
    if n_success == 0
        error("no interval is persistent enough to cover the data")
    elseif n_success < out_dim
        s = n_success == 1 ? " was" : "s were"
        warn && @warn "only $n_success interval$s persistent enough to cover the data"
    end

    meta = (
        filtration=filtration,
        diagram=diagram,
        intervals=diagram[end:-1:(end + 1 - n_success)],
        modulus=modulus,
    )
    return CircularCoordinates(
        n_success, landmarks, partition, metric, coord_data, meta, zeros(length(landmarks))
    )
end

function Base.show(io::IO, cc::CircularCoordinates)
    return print(
        io,
        "CircularCoordinates(\n",
        "  out_dim=$(cc.out_dim),\n",
        "  radius=$(tuple((d.radius for d in cc.coordinate_data)...)),\n",
        "  n_landmarks=$(length(cc.landmarks)),\n",
        "  partition=$(cc.partition),\n",
        "  metric=$(cc.metric),\n)",
    )
end

# Transform single
function _transform(cc::CircularCoordinates, point, dim)
    cd = cc.coordinate_data[dim]

    selected_ball = 0
    min_distance = Inf
    for (j, l) in enumerate(cc.landmarks)
        dist = cc.metric(point, l)
        value = cc.partition(cd.radius, dist)
        if isinf(value)
            cc.partition_buffer .= 0.0
            cc.partition_buffer[j] = 1.0
            selected_ball = j
            break
        else
            cc.partition_buffer[j] = value
        end
        if dist < min_distance
            selected_ball = j
            min_distance = dist
        end
    end
    norm = sum(cc.partition_buffer)
    if norm > 0
        cc.partition_buffer ./= norm
    else
        return missing
    end

    j = selected_ball
    coord = cd.coordinate[j]
    for k in eachindex(cc.landmarks)
        j == k && continue
        cocycle_value = j > k ? -cd.cocycle[index((j, k))] : +cd.cocycle[index((k, j))]
        coord += cc.partition_buffer[k] * cocycle_value
    end
    return _mod_z(coord)
end

function (cc::CircularCoordinates)(points, dims=1:(cc.out_dim))
    if dims isa Integer
        result = zeros(Union{Missing,Float64}, length(points))
    else
        result = zeros(Union{Missing,Float64}, length(points), cc.out_dim)
    end
    for (j, dim) in enumerate(dims)
        for (i, p) in enumerate(points)
            result[i, j] = _transform(cc, SVector(p), dim)
        end
    end
    return result
end
