import MLJModelInterface
const MMI = MLJModelInterface

using PersistenceDiagrams.MLJPersistenceDiagrams
using PersistenceDiagrams: AbstractVectorizer
export RipsPersistentHomology, AlphaPersistentHomology, CubicalPersistentHomology

abstract type RipsererModel <: MMI.Unsupervised end

PointLike{N} = AbstractVector{NTuple{N,MMI.Continuous}}
DistanceMatrix{N} = AbstractMatrix{MMI.Continuous,N}
ImageLike{N} = Union{AbstractArray{MMI.Continuous,N},MMI.Image}

function _transform(model::RipsererModel, verbosity::Int, X)
    result = NamedTuple(Symbol(:dim_, i) => PersistenceDiagram[] for i in 0:(model.dim_max))
    if verbosity > 0
        progbar = Progress(MMI.nrows(X); desc="Computing persistent homology...")
    end
    for x in X
        r = ripserer(_filtration(model, x); _ripserer_args(model)...)
        for i in 1:(model.dim_max + 1)
            push!(result[i], r[i])
        end
        verbosity > 0 && next!(progbar)
    end
    return MMI.table(result)
end

function MMI.fit(model::RipsererModel, verbosity::Int, X)
    diagrams = _transform(model, verbosity, X)
    vectorizers, _, report = MMI.fit(model.vectorizer, verbosity, diagrams)
    return (vectorizers, nothing, report)
end

function MMI.transform(model::RipsererModel, vectorizers, X)
    diagrams = _transform(model, 0, X)
    return MMI.transform(model.vectorizer, vectorizers, diagrams)
end

MMI.output_scitype(::Type{<:RipsererModel}) = MMI.Table(MMI.Continuous)

"""
    RipsPersistentHomology

Compute Vietoris-Rips persistent homology and convert the results to a table of continuous
values.

!!! warning Warning
    Computing Vietoris-Rips persistent homology may be CPU and memory intensive even for
    modestly-sized data sets. Consider using [`AlphaPersistentHomology`](@ref) for
    low-dimensional data.

# Hyperparameters

* `vectorizer = PersistenceImageVectorizer()`: `AbstractVectorizer` used to transform
  persistence diagrams to continuous vectors. See the [PersistenceDiagrams.jl
  Documentation](https://mtsch.github.io/PersistenceDiagrams.jl/dev/mlj/) for more
  information.

* `dim_max = 1`: compute persistent homology up to this dimension.

* `modulus = 2`: compute persistent homology with coefficients in the prime field of integers
  mod `modulus`.

* `threshold = nothing`: compute persistent homology up to diameter smaller than threshold.

* `cutoff = 0`: only keep intervals with `persistence(interval) > cutoff`.

* `sparse = false`: use a sparse Rips filtration.

* `collapse = false`: use the [`EdgeCollapsedRips`](@ref) filtration. This is usually a good
  choice for computations with a high `dim_max`.

# See also

* [`ripserer`](@ref)
* [`Rips`](@ref)
* [`EdgeCollapsedRips`](@ref)

"""
MMI.@mlj_model mutable struct RipsPersistentHomology <: RipsererModel
    vectorizer::AbstractVectorizer = PersistenceImageVectorizer()
    dim_max::Int = 1::(_ ≥ 0)
    modulus::Int = 2::(is_prime(_))
    threshold::Union{Nothing,Float64} = nothing::(isnothing(_) || _ > 0)
    cutoff::Float64 = 0::(_ ≥ 0)
    sparse::Bool = false
    collapse::Bool = false
end

function _filtration(model::RipsPersistentHomology, x)
    if isnothing(model.threshold)
        kwargs = (; threshold=model.threshold, sparse=model.sparse)
    else
        kwargs = (; sparse=model.sparse)
    end
    if model.collapse
        return EdgeCollapsedRips(x; kwargs...)
    else
        return Rips(x; kwargs...)
    end
end

function _ripserer_args(model::RipsPersistentHomology)
    return (; dim_max=model.dim_max, cutoff=model.cutoff, modulus=model.modulus)
end

function MMI.input_scitype(::Type{<:RipsPersistentHomology})
    return Union{
        MMI.Table(Union{PointLike{N},DistanceMatrix}),
        AbstractVector{Union{PointLike{N},DistanceMatrix}},
    } where {N}
end

"""
    AlphaPersistentHomology

Compute alpha complex persistent homology and convert the results to a table of continuous
values.

!!! warning Warning
    Using high-dimensional data with this model may be computationaly expensive. Consider
    using [`RipsPersistentHomology`](@ref).

# Hyperparameters

* `vectorizer = PersistenceImageVectorizer()`: `AbstractVectorizer` used to transform
  persistence diagrams to continuous vectors. See the [PersistenceDiagrams.jl
  Documentation](https://mtsch.github.io/PersistenceDiagrams.jl/dev/mlj/) for more
  information.

* `dim_max = 1`: compute persistent homology up to this dimension.

* `modulus = 2`: compute persistent homology with coefficients in the prime field of integers
  mod `modulus`.

* `threshold = nothing`: compute persistent homology up to diameter smaller than threshold.

* `cutoff = 0`: only keep intervals with `persistence(interval) > cutoff`.

# See also

* [`ripserer`](@ref)
* [`Alpha`](@ref)

"""
MMI.@mlj_model mutable struct AlphaPersistentHomology <: RipsererModel
    vectorizer::AbstractVectorizer = PersistenceImageVectorizer()
    dim_max::Int = 1::(_ ≥ 0)
    modulus::Int = 2::(is_prime(_))
    threshold::Union{Nothing,Float64} = nothing::(isnothing(_) || _ > 0)
    cutoff::Float64 = 0::(_ ≥ 0)
end

function _filtration(model::AlphaPersistentHomology, data)
    if isnothing(model.threshold)
        return Alpha(data)
    else
        return Alpha(data; threshold=model.threshold)
    end
end

function _ripserer_args(model::AlphaPersistentHomology)
    return (; dim_max=model.dim_max, cutoff=model.cutoff, modulus=model.modulus)
end

function MMI.input_scitype(::Type{<:AlphaPersistentHomology})
    return Union{MMI.Table(PointLike{N}),AbstractVector{PointLike{N}}} where {N}
end

"""
    CubicalPersistentHomology

Compute cubical persistent homology and convert the results to a table of continuous
values.

# Hyperparameters

* `vectorizer = PersistenceImageVectorizer()`: `AbstractVectorizer` used to transform
  persistence diagrams to continuous vectors. See the [PersistenceDiagrams.jl
  Documentation](https://mtsch.github.io/PersistenceDiagrams.jl/dev/mlj/) for more
  information.

* `dim_max = 1`: compute persistent homology up to this dimension.

* `threshold = nothing`: compute persistent homology up to diameter smaller than threshold.

* `cutoff = 0`: only keep intervals with `persistence(interval) > cutoff`.

* `negate = false`: negate the image before computation.

# See also

* [`ripserer`](@ref)
* [`Cubical`](@ref)

"""
MMI.@mlj_model mutable struct CubicalPersistentHomology <: RipsererModel
    vectorizer::AbstractVectorizer = PersistenceImageVectorizer()
    dim_max::Int = 1::(_ ≥ 0)
    threshold::Union{Nothing,Float64} = nothing::(isnothing(_) || _ > 0)
    cutoff::Float64 = 0::(_ ≥ 0)
    negate::Bool = false
end

function _filtration(model::CubicalPersistentHomology, data)
    data = model.negate ? -data : data
    if isnothing(model.threshold)
        return Cubical(data)
    else
        return Cubical(data; threshold=model.threshold)
    end
end

function _ripserer_args(model::CubicalPersistentHomology)
    return (; dim_max=model.dim_max, cutoff=model.cutoff)
end

function MMI.input_scitype(::Type{<:CubicalPersistentHomology})
    return Union{MMI.Table(ImageLike{N}),AbstractVector{ImageLike{N}}} where {N}
end

MMI.metadata_pkg.(
    (RipsPersistentHomology, AlphaPersistentHomology, CubicalPersistentHomology),
    name="Ripserer",
    uuid="aa79e827-bd0b-42a8-9f10-2b302677a641",
    url="https://github.com/mtsch/Ripserer.jl",
    license="MIT",
    julia=true,
    is_wrapper=false,
)
