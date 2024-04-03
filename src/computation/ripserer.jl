function _unpack_kwargs(;
    dim_max=1,
    cutoff=0,
    modulus=2,
    field=Mod{modulus},
    verbose=false,
    alg=:cohomology,
    reps=alg == :cohomology ? false : 1:dim_max,
    implicit=alg != :homology,
    # deprecated
    progress=nothing,
    field_type=nothing,
    #
    filtration_kwargs...,
)
    if !isnothing(progress)
        @warn "`progress` is deprecated. Use `verbose` instead" maxlog = 1
        verbose = progress
    end
    if !isnothing(field_type)
        @warn "`field_type` is deprecated. Use `field` instead" maxlog = 1
        field = field_type
    end
    ripserer_kwargs = (
        dim_max=dim_max,
        cutoff=cutoff,
        field=field,
        verbose=verbose,
        alg=alg,
        reps=reps,
        implicit=implicit,
    )

    return ripserer_kwargs, filtration_kwargs
end

"""
    ripserer(Type{<:AbstractFiltration}, args...; kwargs...)
    ripserer(filtration::AbstractFiltration; kwargs...)

Compute the persistent homology of a filtration. The filtration can be given as an
[`AbstractFiltration`](@ref) type, followed by its arguments, or as an initialized object
(see examples below). If only data is given, `Rips` is used by default.

Returns a `Vector` of [`PersistenceDiagram`](@ref)s with (`dim_max` + 1) elements. The
diagrams are sorted by dimension; the first element of the result is the 0-dimensional
diagram, and the last is the (`dim_max`)-dimensional diagram.

# Keyword Arguments

* `dim_max`: compute persistent homology up to this dimension. Defaults to `1`.

* `modulus`: compute persistent homology with coefficients in the prime field of integers
  mod `modulus`. Defaults to `2`.

* `field`: use this type of field of coefficients. Defaults to
  [`Ripserer.Mod`](@ref)`{modulus}`.

* `threshold`: compute persistent homology up to diameter smaller than threshold. Note that
  this parameter is passed to the filtration constructor. When using low thresholds with
  Rips filtrations, consider setting `sparse=true` for optimal performance.

* `cutoff`: only keep intervals with `persistence(interval) > cutoff`. Defaults to `0`. When
  `cutoff < 0`, the result will also contain zero-length intervals.

* `reps`: if `true`, attach representative (co)cycles to persistence intervals. Can also be
  set to collection of integers to only find representatives in specified dimensions,
  e.g. `reps=1:2` will only find representatives in dimensions 1 and 2. This is useful for
  large filtrations (such as cubical) where calculating zero-dimensional representatives can
  be very slow. Defaults to `false` for cohomology and `1:dim_max` for homology.
  Representatives are wrapped in a [`Chain`](@ref).

* `verbose`: If `true`, show a verbose bar. Defaults to `false`.

* `alg`: select the algorithm used in computation. The options are:

  - `:cohomology`: Default and fastest algorithm. When `reps` is set, intervals are equipped
    with representative _co_cycles.

  - `:homology`: Significantly slower than `:cohomology`, but finds representative
    cycles.

  - `:involuted`: Use cohomology result to compute representative cycles. Can be extremely
    efficient compared to `:homology`, especially with `Rips` filtrations. See [this
    paper](https://arxiv.org/abs/2105.03629) for more information.

* `implicit`: If `true`, an implicit reduction algorithm is used. Defaults to `true` for
  :cohomology and `:involuted`, and `false` for `:homology`. `implicit=false` is not
  recommended for `:cohomology` because it disables the emergent pairs optimization.

Other `kwargs...` are passed to the filtration.

# Examples

```jldoctest
julia> ts = range(0, 2π; length=20)[1:(end - 1)];

julia> X = [((2 + cos(θ)) * cos(φ), (2 + cos(θ)) * sin(φ), sin(θ)) for θ in ts for φ in ts];

julia> ripserer(X)
2-element Vector{PersistenceDiagrams.PersistenceDiagram}:
 361-element 0-dimensional PersistenceDiagram
 362-element 1-dimensional PersistenceDiagram

julia> ripserer(EdgeCollapsedRips, X; modulus=7, threshold=2)
2-element Vector{PersistenceDiagrams.PersistenceDiagram}:
 361-element 0-dimensional PersistenceDiagram
 362-element 1-dimensional PersistenceDiagram

julia> ripserer(Rips(X; threshold=1); alg=:involuted)
2-element Vector{PersistenceDiagrams.PersistenceDiagram}:
 361-element 0-dimensional PersistenceDiagram
 362-element 1-dimensional PersistenceDiagram

```

"""
function ripserer(F::Type, args...; kwargs...)
    ripserer_kwargs, filtration_kwargs = _unpack_kwargs(; kwargs...)
    start_time = time_ns()
    filtration = F(args...; verbose=ripserer_kwargs.verbose, filtration_kwargs...)
    return _ripserer(filtration, start_time; ripserer_kwargs...)
end

# Default: Rips
ripserer(dist; kwargs...) = ripserer(Rips, dist; kwargs...)

function ripserer(filtration::AbstractFiltration; kwargs...)
    ripserer_kwargs, other_kwargs = _unpack_kwargs(; kwargs...)
    start_time = time_ns()
    return _ripserer(filtration, start_time; ripserer_kwargs..., other_kwargs...)
end

_reps(reps::Bool, _) = reps
_reps(reps, dim) = dim in reps

function _ripserer(
    filtration::AbstractFiltration,
    start_time;
    dim_max,
    cutoff,
    field,
    verbose,
    alg,
    reps,
    implicit,
)
    if field <: Union{Signed,Unsigned,AbstractFloat}
        error("$field is not a field! Please try a differnet field type")
    end
    index_overflow_check(filtration, field, dim_max)
    result = _ripserer(
        Val(alg), filtration, cutoff, verbose, field, dim_max, reps, implicit
    )
    elapsed = round((time_ns() - start_time) / 1e9; digits=3)
    @prog_println verbose "Done. Time: " ProgressMeter.durationstring(elapsed)

    return result
end

function _ripserer(::Val{A}, args...) where {A}
    return throw(ArgumentError("unsupported alg=$A"))
end
