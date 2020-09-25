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

* `reps`: if `true`, attach representative (co)cycles to persistence intervals. Can also be
  set to collection of integers to only find representatives in specified dimensions,
  e.g. `reps=1:2` will only find representatives in dimensions 1 and 2. This is useful for
  large filtrations (such as cubical) where calculating zero-dimensional representatives can
  be very slow.  Defaults to `false` for cohomology and `1:dim_max` for homology.

* `progress`: If `true`, show a progress bar. Defaults to `false`.

* `metric`: when calculating persistent homology from points, any metric from
  [`Distances.jl`](https://github.com/JuliaStats/Distances.jl) can be used. Defaults to
  `Distances.Euclidean(1e-12)`.

* `alg`: select the algorithm used in computation. The options are:

  - `:cohomology`: Default and fastest algorithm. When `reps` is set, intervals are equipped
    with representative _co_cycles.

  - `:homology`: Significantly slower than `:cohomology`, but finds representative
    cycles. Does not find infinite intervals beyond dimension 0.

  - `:involuted`: Use cohomology result to compute representative cycles. Can be extremely
    efficient compared to `:homology`, especially with `Rips` filtrations. Unlike
    `:homology`, this algorithm finds infinite intervals.

* `implicit`: If `true`, an implicit reduction algorithm is used. Defaults to `true` for
  :cohomology and `:involuted`, and `false` for `:homology`. `implicit=false` is not
  recommended for `:cohomology` because it disables the emergent pairs optimization.

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
    modulus=2,
    field_type=Mod{modulus},
    progress=false,
    alg=:cohomology,
    reps=alg == :cohomology ? false : 1:dim_max,
    implicit=alg != :homology,
)
    start_time = time_ns()
    index_overflow_check(filtration, field_type, dim_max)
    result = _ripserer(
        Val(alg), filtration, cutoff, progress, field_type, dim_max, reps, implicit
    )
    elapsed = round((time_ns() - start_time) / 1e9, digits=3)
    prog_println(progress, "Done. Time: ", ProgressMeter.durationstring(elapsed))

    return result
end

_reps(reps::Bool, _) = reps
_reps(reps, dim) = dim in reps

function _ripserer(
    ::Val{:cohomology}, filtration, cutoff, progress, field_type, dim_max, reps, implicit
)
    result = PersistenceDiagram[]
    zeroth, to_reduce, to_skip = zeroth_intervals(
        filtration, cutoff, progress, field_type, _reps(reps, 0)
    )
    push!(result, zeroth)
    if dim_max > 0
        matrix = CoboundaryMatrix{implicit}(field_type, filtration, to_reduce, to_skip)
        for dim in 1:dim_max
            push!(result, compute_intervals!(matrix, cutoff, progress, _reps(reps, dim)))
            if dim < dim_max
                matrix = next_matrix(matrix, progress)
            end
        end
    end
    return result
end

function _ripserer(
    ::Val{:homology}, filtration, cutoff, progress, field_type, dim_max, reps, implicit
)
    result = PersistenceDiagram[]
    zeroth, to_reduce, to_skip = zeroth_intervals(
        filtration, cutoff, progress, field_type, _reps(reps, 0)
    )
    push!(result, zeroth)
    if dim_max > 0
        simplices = columns_to_reduce(filtration, Iterators.flatten((to_reduce, to_skip)))
        for dim in 1:dim_max
            matrix = BoundaryMatrix{implicit}(field_type, filtration, simplices)
            push!(result, compute_intervals!(matrix, cutoff, progress, _reps(reps, dim)))
            if dim < dim_max
                simplices = columns_to_reduce(filtration, simplices)
            end
        end
    end
    return result
end

function _ripserer(
    ::Val{:involuted}, filtration, cutoff, progress, field_type, dim_max, reps, implicit
)
    result = PersistenceDiagram[]
    zeroth, to_reduce, to_skip = zeroth_intervals(
        filtration, cutoff, progress, field_type, _reps(reps, 0)
    )
    push!(result, zeroth)
    if dim_max > 0
        comatrix = CoboundaryMatrix{true}(field_type, filtration, to_reduce, to_skip)
        for dim in 1:dim_max
            columns, inf_births = compute_death_simplices!(comatrix, progress, cutoff)
            matrix = BoundaryMatrix{implicit}(field_type, filtration, columns)
            diagram = compute_intervals!(matrix, cutoff, progress, _reps(reps, dim))
            for birth_simplex in inf_births
                push!(
                    diagram.intervals,
                    interval(comatrix, birth_simplex, nothing, 0, _reps(reps, dim))
                )
            end
            push!(result, diagram)
            if dim < dim_max
                comatrix = next_matrix(comatrix, progress)
            end
        end
    end
    return result
end

function _ripserer(::Val{A}, args...) where A
    throw(ArgumentError("unsupported alg=$A"))
end
