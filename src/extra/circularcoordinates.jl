export Partition
module Partition

linear(r, d) = max(r - d, 0.0)
gauss(r, d) = exp(-3 * (d - r)^2)

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
    edges_at_t = filter!(e -> birth(e) < time, edges(filtration))
    is = Int[]
    js = Int[]
    vs = Float64[]
    for e in edges_at_t
        i = index(e)
        u, v = e
        append!(js, (u, v))
        append!(is, (i, i))
        append!(vs, (-1, 1))
    end
    n = nv(filtration)
    return sparse(is, js, vs, n * (n - 1), n)
end

function _to_vector(filtration::AbstractFiltration, cocycle::Chain{Int}, time)
    n = nv(filtration)
    vector = zeros(n * (n - 1))
    for (sx, c) in cocycle
        if birth(sx) < time
            vector[index(sx)] = Float64(c)
        end
    end
    return vector
end

function _mod_z(t)
    t = rem(t, 1)
    return t < 0 ? 1 + t : t
end

function _harmonic_smoothing(filtration, chain; time)
    integral_cocycle = _to_integer_coefficients(chain)
    if !is_cocycle(filtration, integral_cocycle, time)
        error(
            "the cocycle cannot be converted to `Int` coefficients. ",
            "Try using a different `modulus`.",
        )
    end
    coboundary = _zero_coboundary_matrix(filtration, time)
    real_cocycle = _to_vector(filtration, integral_cocycle, time)

    minimizer = _mod_z.(coboundary \ real_cocycle)
    harmonic_cocycle = real_cocycle - coboundary * minimizer
    return minimizer, harmonic_cocycle
end

struct CircularCoordinateData
    radius::Float64
    coordinate::Vector{Float64}
    cocycle::Vector{Float64}
end

struct CircularCoordinates{F, P<:SVector, M}
    out_dim::Int
    landmarks::Vector{P}
    partition_function::F
    metric::M

    coordinate_data::Vector{CircularCoordinateData}

    meta::NamedTuple
    #
    partition::Vector{Float64}
end

function CircularCoordinates(args...; kwargs...)
    return CircularCoordinates(Rips, args...; kwargs...)
end
function CircularCoordinates(
    ::Type{F}, points::AbstractVector, landmarks=eachindex(points);
    out_dim=1,
    modulus=23,
    metric=Euclidean(),
    threshold=nothing,
    coverage=1,
    partition=Partition.gauss,
    progress=false,
    kwargs...,
) where {F<:AbstractFiltration}
    @prog_print progress "Determining radius... "
    landmarks, min_radius = _landmarks_and_radius(points, landmarks, metric)
    @prog_println progress "done."

    # TODO: new ripserer interface, extract filtration from diagram
    flt_kwargs = metric == Euclidean() ? NamedTuple() : (metric=Euclidean())
    flt_kwargs = isnothing(threshold) ? flt_kwargs : (; threshold=threshold, flt_kwargs...)

    # compute cohomology
    @prog_print progress "Computing cohomology... "
    filtration = F(landmarks; flt_kwargs...)
    diagram = ripserer(filtration; modulus=modulus, reps=true, kwargs...)[2]
    @prog_println progress "done."
    if length(diagram) < out_dim
        error("diagram has $(length(diagram)) intervals")
    end

    @prog_print progress "Massaging cocycles... "
    coord_data = map(1:out_dim) do d
        interval = diagram[end + 1 - d]
        b, d = interval
        if !isfinite(death(interval))
            d = Float64(maximum(birth, interval.cocycle))
        elseif max(b, min_radius) ≥ d / 2
            @warn(
                "landmarks do not cover the points well enough.\n",
                "interval: $interval\n",
                "max distance to landmarks: $min_radius",
            )
            min_radius = -Inf
        end
        radius = coverage * max(b, min_radius) + (1 - coverage) * d / 2

        # Smoothen the cocycle.
        coords, cocycle = _harmonic_smoothing(
            filtration, interval.representative; time=2 * radius
        )
        CircularCoordinateData(radius, coords, cocycle)
    end
    @prog_println progress "done."

    meta = (
        filtration=filtration,
        diagram=diagram,
        intervals=diagram[end:-1:(end + 1 - out_dim)],
        modulus=modulus,
    )

    return CircularCoordinates(
        out_dim,
        landmarks,
        partition,
        metric,
        coord_data,
        meta,
        zeros(length(landmarks)),
    )
end

function Base.show(io::IO, cc::CircularCoordinates)
    return print(
        io,
        "CircularCoordinates(\n",
        "  out_dim=$(cc.out_dim),\n",
        "  radii=$(tuple((d.radius for d in cc.coordinate_data)...)),\n",
        "  n_landmarks=$(length(cc.landmarks)),\n",
        "  partition=$(cc.partition_function),\n",
        "  metric=$(cc.metric),\n)",
    )
end

function _transform(cc::CircularCoordinates, point, dim)
    cd = cc.coordinate_data[dim]

    nearest_index = 0
    nearest_dist = Inf
    for (j, l) in enumerate(cc.landmarks)
        dist = cc.metric(point, l)
        cc.partition[j] = cc.partition_function(cd.radius, dist)
        if dist < nearest_dist
            nearest_index = j
            nearest_dist = dist
        end
    end
    norm = sum(cc.partition)
    if norm > 0
        cc.partition ./= norm
    else
        @warn "partition of unity summed to 0." maxlog=1
    end

    j = nearest_index
    coord = cd.coordinate[j]
    for k in eachindex(cc.landmarks)
        cocycle_value = j > k ? cd.cocycle[index((j, k))] : -cd.cocycle[index((k, j))]
        coord += cc.partition[k] + cocycle_value
    end
    return _mod_z(coord)
end

function (cc::CircularCoordinates)(points, dims=1:cc.out_dim)
    if length(dims) > 1
        result = zeros(length(points), cc.out_dim)
    else
        result = zeros(length(points))
    end
    for (j, dim) in enumerate(dims)
        for (i, p) in enumerate(points)
            result[i, j] = _transform(cc, SVector(p), dim)
        end
    end
    return result
end
