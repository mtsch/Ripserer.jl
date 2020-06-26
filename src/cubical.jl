"""
    Cubelet{D, T, I} <: IndexedSimplex{D, T, I}

A `Cubelet` is similar to a `Simplex`, but it has `2^D` vertices instead of `D+1`. Like in
`Simplex`, the vertices are encoded from an index and dimension. Because a cubelet knows
nothing about the image it came from, it returns *linear* indices from `vertices`.

The vertices should be neighboring indices, but this fact is not checked anywhere.
"""
struct Cubelet{D, T, I} <: IndexedSimplex{D, T, I}
    # This is not the most efficient way to index Cubelets, since most of the indices remain
    # unused. It forces us to use large Int types for decently sized images. It does,
    # however, make comparisons very fast, which is good.
    index::I
    diam::T

    function Cubelet{D, T, I}(index::Integer, diam) where {D, T, I<:Integer}
        D ≥ 0 || throw(DomainError(D, "dimension must be a non-negative integer"))
        return new{D, T, I}(I(index), T(diam))
    end
end

function Cubelet{D}(index::I, diam::T) where {D, T, I<:Integer}
    return Cubelet{D, T, I}(index, diam)
end
function Cubelet{D}(vertices, diam::T) where {D, T}
    return Cubelet{D, T, eltype(vertices)}(vertices, diam)
end
@generated function Cubelet{D, T, I}(vertices, diam) where {D, T, I<:Integer}
    K = 2^D
    return quote
        length(vertices) == $K ||
            throw(ArgumentError(string("a `Cubelet{",$D,"}` must have ",$K," simplices")))
        vertices_svec = sort(SVector{$K, $I}(vertices), rev=true)
        return Cubelet{$D, $T, $I}(index(vertices_svec), $T(diam))
    end
end

index(cube::Cubelet) = cube.index
diam(cube::Cubelet) = cube.diam

# @generated used because inference doesn't work well for 2^D
@generated function vertices(cube::Cubelet{D, <:Any, I}) where {D, I}
    return :(vertices(index(cube), Val($(2^D))))
end

Base.length(::Type{<:Cubelet{D}}) where D = 1 << D

"""
    Cubical{I<:Signed, T} <: AbstractFiltration

`Cubical` is used to compute sublevel persistent homology on `N`-dimensional images, which
are of type `AbstractArray{T, N}`. `I` is the integer type used to represent cubelet
indices. Set it to a large value for larger images.

# Constructor

* `Cubical{I}(::AbstractArray{T, N})`: `I` defaults to `Int64`.
"""
struct Cubical{I<:Signed, T, A<:AbstractArray{T}} <: AbstractFiltration
    data::A
    threshold::T
end

function Cubical(data)
    return Cubical{Int}(data)
end
function Cubical{I}(data::AbstractArray{T}) where {I, T}
    return Cubical{I, T, typeof(data)}(data, maximum(data))
end

n_vertices(cf::Cubical) = length(cf.data)
threshold(cf::Cubical) = cf.threshold
birth(cf::Cubical, i) = cf.data[i]
simplex_type(::Cubical{I, T}, dim) where {I, T} = Cubelet{dim, T, I}

dim(cf::Cubical) = length(size(cf.data))

function to_linear(cf::Cubical{I}, vertices) where I
    indices = LinearIndices(cf.data)
    return map(v -> I(get(indices, v, 0)), vertices)
end
function to_cartesian(cf::Cubical, vertices)
    indices = CartesianIndices(cf.data)
    return map(v -> indices[v], vertices)
end

function unsafe_simplex(
    cf::Cubical{I, T}, ::Val{D}, vertices, sign=one(I)
) where {I, T, D}
    diam = typemin(T)
    for v in vertices
        d = get(cf.data, v, missing)
        if ismissing(d)
            return nothing
        else
            _d::T = d
            diam = ifelse(_d > diam, _d, diam)
        end
    end
    return simplex_type(cf, D)(sign * index(vertices), diam)
end

# Check if linear indices u and v are adjacent in array arr.
function is_adjacent(arr, u, v)
    if isassigned(arr, u) && isassigned(arr, v)
        diff = Tuple(CartesianIndices(arr)[u] - CartesianIndices(arr)[v])
        ones = abs.(diff) .== 1
        zeros = diff .== 0
        return count(ones) == 1 && count(zeros) == length(diff) - 1
    else
        return false
    end
end

function edges(cf::Cubical{I, <:Any}) where {I}
    result = edge_type(cf)[]
    for u in eachindex(cf.data)
        for d in 1:dim(cf)
            if d == 1
                v = u + 1
            else
                v = u + prod(size(cf.data, i) for i in 1:d - 1)
            end
            if is_adjacent(cf.data, u, v)
                sx = unsafe_simplex(cf, Val(1), (v, u), 1)
                if !isnothing(sx)
                    _sx::edge_type(cf) = sx
                    push!(result, _sx)
                end
            end
        end
    end
    return result
end

# coboundary ============================================================================= #
struct CubeletCoboundary{A, D, N, I, F<:Cubical{I}, K}
    filtration::F
    vertices::NTuple{K, CartesianIndex{N}}
end

function CubeletCoboundary{A}(
    filtration::F, cubelet::C
) where {A, I, D, F<:Cubical{I}, C<:Cubelet{D}}
    K = length(C)
    N = dim(filtration)
    return CubeletCoboundary{A, D, N, I, F, K}(
        filtration, to_cartesian(filtration, Tuple(vertices(cubelet)))
    )
end

function coboundary(filtration::Cubical, cubelet::Cubelet)
    return CubeletCoboundary{true}(filtration, cubelet)
end
function coboundary(filtration::Cubical, cubelet::Cubelet, ::Val{false})
    return CubeletCoboundary{false}(filtration, cubelet)
end

function all_equal_in_dim(dim, vertices)
    i1 = vertices[1][dim]
    for v in vertices
        v[dim] ≠ i1 && return false
    end
    return true
end

# Type safe way to get half of tuple a-la TupleTools.
function second_half(tup::NTuple{N}) where N
    return _half(tup, Val(N÷2), TupleTools.unsafe_tail)
end
function first_half(tup::NTuple{N}) where N
    return _half(tup, Val(N÷2), TupleTools.unsafe_front)
end
@inline function _half(tup, ::Val{N}, f) where N
    if N == 0
        return tup
    else
        return _half(f(tup), Val(N-1), f)
    end
end

function Base.iterate(
    cc::CubeletCoboundary{A, D, N, I}, (dim, dir)=(one(I), one(I))
) where {A, D, N, I}
    while true
        # Idea: expand cube by adding vertices that have ±1 added to one of the dimensions
        # where all old vertices have the same value.
        while dim ≤ N && !all_equal_in_dim(dim, cc.vertices)
            dim += one(I)
        end
        dim > N && return nothing

        diff = CartesianIndex{N}(ntuple(isequal(dim), Val(N)) .* dir)
        new_vertices = TupleTools.sort(to_linear(
            cc.filtration, TupleTools.vcat(cc.vertices, cc.vertices .+ Ref(diff))
        ), rev=true)
        dim += I(dir == -1)
        dir *= -one(I)
        if A || second_half(new_vertices) == to_linear(cc.filtration, cc.vertices)
            sx = unsafe_simplex(cc.filtration, Val(D + 1), new_vertices, -dir)
            if !isnothing(sx)
                _sx::simplex_type(cc.filtration, D + 1) = sx
                return _sx, (dim, dir)
            end
        end
    end
end

struct CubeletBoundary{D, I, N, F<:Cubical{I}, K}
    filtration::F
    vertices::NTuple{K, CartesianIndex{N}}
end

function CubeletBoundary(
    filtration::F, cubelet::C
) where {D, I, F<:Cubical{I}, C<:Cubelet{D}}
    K = length(cubelet)
    N = dim(filtration)
    return CubeletBoundary{D, I, N, F, K}(
        filtration, to_cartesian(filtration, Tuple(vertices(cubelet)))
    )
end

boundary(filtration::Cubical, cubelet::Cubelet) = CubeletBoundary(filtration, cubelet)

# Idea: split the cube in each dimension, returning one of two halves depending on dir.
function Base.iterate(cb::CubeletBoundary{D, I}, (dim, dir)=(1, 1)) where {D, I}
    if dim > D
        return nothing
    else
        lin_vertices = to_linear(
            cb.filtration, TupleTools.sort(cb.vertices, by=v -> v[dim], rev=true)
        )
        if dir == 1
            new_vertices = first_half(lin_vertices)
            sign = one(I)
            state = (dim, -dir)
        else
            new_vertices = second_half(lin_vertices)
            sign = -one(I)
            state = (dim + 1, 1)
        end
        sx = simplex(cb.filtration, Val(D - 1), new_vertices, sign)
        if !isnothing(sx)
            _sx::simplex_type(cb.filtration, D - 1) = sx
            return _sx, state
        end
    end
end
