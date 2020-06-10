"""
    ripserer(dists::AbstractMatrix; kwargs...)
    ripserer(points; metric=Euclidean(), births, kwargs...)
    ripserer(filtration::AbstractFiltration; kwargs...)

Compute the persistent homology of metric space represented by `dists`, `points` and
`metric` or an `::AbstractFiltration`.

If using points, `points` must be an array of bitstypes, such as `NTuple`s or `SVectors`.

# Keyoword Arguments

* `dim_max`: compute persistent homology up to this dimension. Defaults to `1`.
* `modulus`: compute persistent homology with coefficients in the prime field of integers
  mod `modulus`. Defaults to `2`.
* `field_type`: use this type of field of coefficients. Defaults to `Mod{modulus}`.
* `threshold`: compute persistent homology up to diameter smaller than threshold.
  For non-sparse Rips filtrations, it defaults to radius of input space.
* `cutoff`: only keep intervals with `persistence(interval) > cutoff`. Defaults to `0`.
* `reps`: if `true`, return representative cocycles along with persistence intervals.
  Defaults to `false`.
* `progress`: If `true`, show a progress bar. Defaults to `false`.
* `metric`: when calculating persistent homology from points, any metric from
  [`Distances.jl`](https://github.com/JuliaStats/Distances.jl) can be used. Defaults to
  `Euclidean()`.
* `births`: when calculating persistent homology from points, births can be used to add
  birth times to vertices. Defaults to all births equal to `0`.
"""
function ripserer(
    dists::AbstractMatrix;
    dim_max=1,
    sparse=false,
    cutoff=0,
    reps=false,
    modulus=2,
    field_type=Mod{modulus},
    progress=false,
    cohomology=true,
    kwargs..., # kwargs for filtration
)
    if issparse(dists)
        filtration = SparseRips(dists; kwargs...)
    else
        filtration = Rips(dists; kwargs...)
    end
    return ripserer(
        filtration;
        dim_max=dim_max,
        reps=reps,
        cutoff=cutoff,
        field_type=field_type,
        progress=progress,
        cohomology=cohomology,
    )
end

function ripserer(points; metric=Euclidean(), births=nothing, kwargs...)
    dists = distances(metric, points, births)
    return ripserer(dists; kwargs...)
end

function ripserer(
    filtration::AbstractFiltration;
    dim_max=1, reps=false, cutoff=0, field_type=Mod{2}, progress=false, cohomology=true
)
    if cohomology
        return _cohomology(
            filtration, cutoff, progress, field_type, Val(dim_max), Val(reps)
        )
    else
        return _homology(
            filtration, cutoff, progress, field_type, Val(dim_max), Val(reps)
        )
    end
end

function _cohomology(
    filtration, cutoff, progress, ::Type{F}, ::Val{dim_max}, ::Val{reps}
) where {F, dim_max, reps}
    result = PersistenceDiagram[]
    zeroth, to_reduce, to_skip = zeroth_intervals(
        filtration, cutoff, progress, F, Val(reps)
    )
    push!(result, zeroth)
    if dim_max == 0
        return result
    else
        matrix = ReductionMatrix{true, F}(filtration, to_reduce, to_skip)
        for dim in 1:dim_max
            push!(result, compute_intervals!(matrix, cutoff, progress, Val(reps)))

            if dim < dim_max
                matrix = next_matrix(matrix, progress)
            end
        end

        return result
    end
end

# Homology is still experimental and does not always work correctly. It does not yet output
# infinite intervals.
function _homology(
    filtration, cutoff, progress, ::Type{F}, ::Val{dim_max}, ::Val{reps}
) where {F, dim_max, reps}
    result = PersistenceDiagram[]
    zeroth, to_reduce, to_skip = zeroth_intervals(
        filtration, cutoff, progress, F, Val(reps)
    )
    push!(result, zeroth)

    # We want to start with triangles. Starting at the highest dimension and going down, as
    # dictated by the twist algorithm might not be worth it. The highest dimension may be
    # sparse and twist is not supposed to bring a huge improvement to regular persistent
    # homology.
    # Constructing a matrix and throwing it away has some overhead, but is nothing compared
    # to the slowness of homology.
    matrix = next_matrix(
        ReductionMatrix{false, F}(filtration, to_reduce, to_skip), progress
    )
    for dim in 1:dim_max
        push!(result, compute_intervals!(matrix, cutoff, progress, Val(reps)))

        if dim < dim_max
            matrix = next_matrix(matrix, progress)
        end
    end

    return result
end
