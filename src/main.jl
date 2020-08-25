"""
    overflows(filtration::AbstractFiltration, dim_max)
    overflows(::Type{AbstractSimplex}, nv, field)

Check if `dim`-dimensional simplex with largest index on `nv(filtration)` can safely
be constructed and if it can be represented as a chain element.
"""
function overflows(flt::AbstractFiltration, dim_max, field)
    return overflows(simplex_type(flt, dim_max + 1), nv(flt), field)
end

overflows(::Type{<:AbstractSimplex}, ::Any, ::Any) = false

function overflows(S::Type{<:Simplex{<:Any, T, I}}, nv, field) where {T, I}
    length(S) > nv && throw(
        ArgumentError("$S has more than $(nv) vertices.")
    )
    # Calculate index for last possible simplex in I and BigInt and check if they are equal.
    vertices = ntuple(i -> I(nv - i + 1), length(S))
    index_I = index(vertices)
    index_big = index(BigInt.(vertices))
    if index_I ≠ index_big
        return true
    else
        # Check that packing is safe.
        Element = chain_element_type(S, field)
        return index(simplex(Element(S(index_I, oneunit(T))))) ≠ index_big
    end
end

"""
    ripserer(dists::AbstractMatrix; kwargs...)
    ripserer(points; metric=Distances.Euclidean(), births, kwargs...)
    ripserer(filtration::AbstractFiltration; kwargs...)

Compute the persistent homology of metric space represented by `dists`, `points` and
`metric` or a [`Ripserer.AbstractFiltration`](@ref).

If using points, `points` must be an array of `isbits` types, such as `NTuple`s or
`SVector`s.

# Keyoword Arguments

* `dim_max`: compute persistent homology up to this dimension. Defaults to `1`.
* `modulus`: compute persistent homology with coefficients in the prime field of integers
  mod `modulus`. Defaults to `2`.
* `field_type`: use this type of field of coefficients. Defaults to
  [`Ripserer.Mod`](@ref)`{modulus}`.
* `threshold`: compute persistent homology up to diameter smaller than threshold. This
  parameter is only applicable when using distance matrices or points as input. When using
  filtrations, threshold must be passed to the filtration constructor. Defaults to radius of
  input space. When using low thresholds with points or distance matrices, consider using
  [`SparseRips`](@ref).
* `cutoff`: only keep intervals with `persistence(interval) > cutoff`. Defaults to `0`.
* `reps`: if `true`, return representative cocycles along with persistence intervals.
  Defaults to `false`.
* `progress`: If `true`, show a progress bar. Defaults to `false`.
* `metric`: when calculating persistent homology from points, any metric from
  [`Distances.jl`](https://github.com/JuliaStats/Distances.jl) can be used. Defaults to
  `Distances.Euclidean()`.
* `cohomology`: if set to `false`, compute persistent homology instead of cohomology. This
  is much slower and gives the same result, but may give more informative representatives
  when `reps` is enabled. Currently unable to compute infinite intervals in dimensions
  higher than 0. Defaults to `false`.
"""
function ripserer(
    dists::AbstractMatrix;
    threshold=nothing,
    kwargs...,
)
    if issparse(dists)
        filtration = SparseRips(dists; threshold=threshold)
    else
        filtration = Rips(dists; threshold=threshold)
    end
    return ripserer(filtration; kwargs...)
end

function ripserer(points::AbstractVector; metric=Euclidean(), threshold=nothing, kwargs...)
    return ripserer(Rips(points; metric=metric, threshold=threshold); kwargs...)
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
    if overflows(filtration, dim_max, field_type)
        S = simplex_type(filtration, dim_max + 1)
        throw(OverflowError(
            "$S on $(nv(filtration)) vertices overflows. " *
            "Try using a larger index type or a smaller `dim_max`."
        ))
    end
    return if cohomology
        _cohomology(filtration, cutoff, progress, field_type, Val(dim_max), Val(reps))
    else
        _homology(filtration, cutoff, progress, field_type, Val(dim_max), Val(reps))
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
