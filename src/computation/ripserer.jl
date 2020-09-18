function cohomology(
    filtration, cutoff, progress, ::Type{F}, ::Val{dim_max}, ::Val{reps}
) where {F, dim_max, reps}
    result = PersistenceDiagram[]
    zeroth, to_reduce, to_skip = zeroth_intervals(
        filtration, cutoff, progress, F, Val(reps)
    )
    push!(result, zeroth)
    if dim_max > 0
        matrix = CoboundaryMatrix(F, filtration, to_reduce, to_skip)
        for dim in 1:dim_max
            push!(result, compute_intervals!(matrix, cutoff, progress, Val(reps)))
            if dim < dim_max
                matrix = next_matrix(matrix, progress)
            end
        end
    end
    return result
end

function homology(
    filtration, cutoff, progress, ::Type{F}, ::Val{dim_max}, ::Val{reps}
) where {F, dim_max, reps}
    result = PersistenceDiagram[]
    zeroth, to_reduce, to_skip = zeroth_intervals(
        filtration, cutoff, progress, F, Val(reps)
    )
    push!(result, zeroth)
    if dim_max > 0
        simplices = columns_to_reduce(filtration, Iterators.flatten((to_reduce, to_skip)))
        for dim in 1:dim_max
            if isempty(simplices)
                # not correct
                return result
            end
            matrix = BoundaryMatrix(F, filtration, simplices)
            push!(result, compute_intervals!(matrix, cutoff, progress, Val(true)))
            if dim < dim_max
                simplices = columns_to_reduce(filtration, simplices)
            end
        end
    end
    return result
end

"""
    ripserer(dists::AbstractMatrix; kwargs...)
    ripserer(points; metric=Distances.Euclidean(1e-12), births, kwargs...)
    ripserer(filtration::AbstractFiltration; kwargs...)

Compute the persistent homology of metric space represented by `dists`, `points` and
`metric`, or a [`Ripserer.AbstractFiltration`](@ref).

If `dists` or `points` are given, the `Rips` filtration is used.

If using points, `points` must be an array of `isbits` types, such as `NTuple`s or
`SVector`s.

# Keyword Arguments

* `dim_max`: compute persistent homology up to this dimension. Defaults to `1`.
* `modulus`: compute persistent homology with coefficients in the prime field of integers
  mod `modulus`. Defaults to `2`.
* `field_type`: use this type of field of coefficients. Defaults to
  [`Ripserer.Mod`](@ref)`{modulus}`.
* `threshold`: compute persistent homology up to diameter smaller than threshold. This
  parameter is only applicable when using distance matrices or points as input. When using
  filtrations, threshold must be passed to the filtration constructor. Defaults to the
  radius of the input space. When using low thresholds with points or distance matrices,
  consider using `sparse=true`.
* `cutoff`: only keep intervals with `persistence(interval) > cutoff`. Defaults to `0`.
* `reps`: if `true`, return representative cocycles along with persistence intervals.
  Defaults to `false`.
* `progress`: If `true`, show a progress bar. Defaults to `false`.
* `metric`: when calculating persistent homology from points, any metric from
  [`Distances.jl`](https://github.com/JuliaStats/Distances.jl) can be used. Defaults to
  `Distances.Euclidean(1e-12)`.
* `cohomology`: if set to `false`, compute persistent homology instead of cohomology. This
  is much slower and gives the same result, but may give more informative representatives
  when `reps` is enabled. Currently unable to compute infinite intervals in dimensions
  higher than 0. Defaults to `false`.
"""
function ripserer(
    dists::AbstractMatrix;
    threshold=nothing, sparse=false, kwargs...,
)
    filtration = Rips(dists; threshold=threshold, sparse=sparse)
    return ripserer(filtration; kwargs...)
end

function ripserer(
    points::AbstractVector;
    metric=Euclidean(1e-12), threshold=nothing, sparse=false, kwargs...
)
    filtration = Rips(points; metric=metric, threshold=threshold, sparse=sparse)
    return ripserer(filtration; kwargs...)
end

function ripserer(
    filtration::AbstractFiltration;
    dim_max=1,
    cutoff=0,
    reps=false,
    modulus=2,
    field_type=Mod{modulus},
    progress=false,
    cohomology=true,
)
    start_time = time_ns()
    index_overflow_check(filtration, field_type, dim_max)
    if cohomology
        result = Ripserer.cohomology(
            filtration, cutoff, progress, field_type, Val(dim_max), Val(reps)
        )
    else
        result = homology(
            filtration, cutoff, progress, field_type, Val(dim_max), Val(reps)
        )
    end
    elapsed = round((time_ns() - start_time) / 1e9, digits=3)
    prog_println(progress, "Done. Took ", elapsed, " seconds.")

    foreach(sort!, result)
    return result
end
