abstract type RipsererModel <: MMI.Unsupervised end

PointLike = AbstractVector{Tuple{Vararg{MMI.Continuous}}}
DistanceMatrix = AbstractMatrix{MMI.Continuous}
ImageLike = Union{AbstractArray{MMI.Continuous, N} where N, MMI.Image}

MMI.fit(model::RipsererModel, _::Int, X) = (nothing, nothing, NamedTuple())

function MMI.transform(model::RipsererModel, ::Nothing, X)
    result = NamedTuple(Symbol(:dim_, i) => PersistenceDiagram[] for i in 0:model.dim_max)
    for x in X
        r = ripserer(_filtration(model, x); dim_max=model.dim_max, modulus=model.modulus)
        for i in 1:model.dim_max+1
            push!(result[i], r[i])
        end
    end
    return MMI.table(result)
end

MMI.output_scitype(::Type{<:RipsererModel}) = MMI.Table(PersistenceDiagram)

"""
"""
MMI.@mlj_model mutable struct RipsPersistentHomology <: RipsererModel
    dim_max::Int = 1::(_ ≥ 0)
    modulus::Int = 2::(is_prime(_))
    threshold::Union{Nothing, Float64} = nothing
    cutoff::Float64 = 0::(_ ≥ 0)
    sparse::Bool = false
end

function _filtration(model::RipsPersistentHomology, data)
    if isnothing(model.threshold)
        return Rips(data; sparse=model.sparse)
    else
        return Rips(data; sparse=model.sparse, threshold=model.threshold)
    end
end

MMI.input_scitype(::Type{<:RipsPersistentHomology}) = Union{
    MMI.Table(Union{PointLike, DistanceMatrix}),
    AbstractVector{Union{PointLike, DistanceMatrix}},
}

"""
"""
MMI.@mlj_model mutable struct AlphaPersistentHomology <: RipsererModel
    dim_max::Int = 1::(_ ≥ 0)
    modulus::Int = 2::(is_prime(_))
    threshold::Union{Nothing, Float64} = nothing
    cutoff::Float64 = 0::(_ ≥ 0)
end

function _filtration(model::AlphaPersistentHomology, data)
    if isnothing(model.threshold)
        return Alpha(data)
    else
        return Alpha(data; threshold=model.threshold)
    end
end

MMI.input_scitype(::Type{<:AlphaPersistentHomology}) = Union{
    MMI.Table(PointLike),
    AbstractVector{PointLike},
}

"""
"""
MMI.@mlj_model mutable struct CubicalPersistentHomology <: RipsererModel
    dim_max::Int = 1::(_ ≥ 0)
    modulus::Int = 2::(is_prime(_))
    threshold::Union{Nothing, Float64} = nothing
    cutoff::Float64 = 0::(_ ≥ 0)
end

function _filtration(model::CubicalPersistentHomology, data)
    if isnothing(model.threshold)
        return Cubical(data)
    else
        return Cubical(data; threshold=model.threshold)
    end
end

MMI.input_scitype(::Type{<:AlphaPersistentHomology}) = Union{
    MMI.Table(ImageLike),
    AbstractVector{ImageLike},
}
