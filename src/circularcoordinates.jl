# TODO: Perhaps this should be done elsewhere.
# TODO: Maybe field_type should be coefficient_type
function to_integer_coefficients(chain::Chain{Mod{M},S}) where {M,S}
    if M == 2
        # TODO message
        error()
    end
    result = Chain{Int,S}()
    resize!(result, length(chain))
    for (i, (sx, coef)) in enumerate(chain)
        new_coef = Int(coef) > (M - 1) ÷ 2 ? Int(coef) - M : Int(coef)
        @inbounds result[i] = (sx, Int(new_coef))
    end
    return result
end

"""
    is_cocycle(filtration, chain::Chain{F,S}, t) where {F,S}

Test whether `chain` is a cocycle in `filtration` at time `t`.
"""
function is_cocycle(filtration, chain::Chain{F,S}, time) where {F,S}
    buffer = Chain{F,simplex_type(filtration, dim(S) + 1)}()
    for (simplex, coefficient) in chain
        for cofacet in coboundary(filtration, simplex)
            if birth(cofacet) < time
                heappush!(buffer, (cofacet, coefficient), Base.Order.Forward)
            end
        end
    end
    # Handle numerical errors with Float64
    pivot = heappop!(buffer, Base.Order.Forward)
    while !isnothing(pivot) && coefficient(pivot) < √eps(Float64)
        pivot = heappop!(buffer, Base.Order.Forward)
    end
    return isnothing(pivot)
end

"""
    zero_coboundary_matrix(filtration::AbstractFiltration, time)

Return the zeroth coboundary matrix as a sparse matrix of Float64.
"""
function zero_coboundary_matrix(filtration::AbstractFiltration, time)
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

function cocycle_to_vector(filtration::AbstractFiltration, cocycle::Chain{Int}, time)
    n = nv(filtration)
    vector = zeros(n * (n - 1))
    for (sx, c) in cocycle
        if birth(sx) < time
            vector[index(sx)] = Float64(c)
        end
    end
    return vector
end

function vector_to_cocycle(filtration, vector)
    S = simplex_type(filtration, 1)
    cocycle = Chain{Float64, S}()
    for (i, c) in enumerate(vector)
        if c ≠ 0
            u, v = _vertices(i, Val(2))
            push!(cocycle, (simplex(filtration, Val(1), (u, v)), c))
        end
    end
    return cocycle
end

function mod_z(t)
    t = rem(t, 1)
    return t < 0 ? 1 + t : t
end

function harmonic_cocycle(
    filtration, interval::PersistenceInterval; time=death(interval), kwargs...
)
    return harmonic_cocycle(filtration, interval.representative; time=time, kwargs...)
end

function harmonic_cocycle(filtration, chain; time, check_cocycle=true)
    integral_cocycle = to_integer_coefficients(chain)
    if check_cocycle && !is_cocycle(filtration, integral_cocycle, time)
        error("no longer a cocycle in `Int`. Try using a different `modulus`.")
    end
    coboundary = zero_coboundary_matrix(filtration, time)
    real_cocycle = cocycle_to_vector(filtration, integral_cocycle, time)

    minimizer = mod_z.(coboundary \ real_cocycle)
    harmonic_cocycle = real_cocycle - coboundary * minimizer
    return minimizer, harmonic_cocycle
end

export circular_coordinates

# TODO: mention picking a later time a good idea in docs.
# TODO: mention that time is an open interval.
function circular_coordinates(points, landmarks=eachindex(points); kwargs...)
    return circular_coordinates(Rips, points, landmarks; kwargs...)
end

function circular_coordinates(
    F::Type, points, landmarks::AbstractVector{Int}=eachindex(points);
    coverage=1,
    n_intervals=1,
    progress=false,
    modulus=29,
)
    allunique(landmarks) || throw(ArgumentError("landmarks not unique"))

    filtration = F(points[landmarks])
    intervals = ripserer(filtration; progress=progress, modulus=modulus, reps=true)[2]
    selected = intervals[end:-1:(end - n_intervals + 1)]

    result = zeros(length(points), n_intervals)

    if landmarks == eachindex(points)
        for i in 1:n_intervals
            interval = selected[i]
            time = birth(interval) + coverage * persistence(interval)
            result[:, i] .= circular_coordinates(filtration, interval; time=time)
        end
    else
        for i in 1:n_intervals
            interval = selected[i]
            result[:, i] .= circular_coordinates(
                points, points[landmarks], filtration, interval;
                coverage=coverage,
                metric=Euclidean(),
            )
        end
    end
    return result
end

function circular_coordinates(
    filtration::AbstractFiltration, interval; time=death(interval), kwargs...
)
    return harmonic_cocycle(filtration, interval.representative; time=time, kwargs...)[1]
end

function circular_coordinates(
    points, landmarks, filtration, interval; coverage=0.999, metric=Euclidean()
)
    dists_to_landmarks = [
        metric(SVector(l), SVector(p)) for l in landmarks, p in points
    ]

    # Determine ball radius.
    rL = maximum(minimum(dists_to_landmarks; dims=1))
    birth, death = interval
    if max(birth, rL) ≥ death / 2
        error(
            "landmarks do not cover the points well enough.\n",
            "interval: $interval\n",
            "max distance to landmarks: $rL",
        )
    end
    radius = coverage * max(birth, rL) + (1 - coverage) * death / 2

    # Smoothen the cocycle.
    τ, harmonic = harmonic_cocycle(filtration, interval; time=2 * radius)

    circular_coordinates = zeros(length(points))
    partition = zeros(length(landmarks))

    # TODO: this could be hidden away in a struct to allow transforming new data.
    for (i, p) in enumerate(points)
        for j in eachindex(landmarks)
            partition[j] = max(0.0, radius - dists_to_landmarks[j, i])
        end
        partition ./= sum(partition)
        # Choose an arbitrary ball
        j = findfirst(≤(radius), view(dists_to_landmarks, :, i))
        if isnothing(j)
            @warn "Can't find a ball to cover $(points[i])"
            circular_coordinates[i] = NaN # TODO: is NaN a good idea?
        else
            θ = τ[j]
            for k in eachindex(landmarks)
                θ += partition[k] + (j > k ? harm[index((j, k))] : -harm[index((k, j))])
            end
            circular_coordinates[i] = mod_z(θ)
        end
    end
    return circular_coordinates
end
