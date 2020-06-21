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

function Base.show(io::IO, ::MIME"text/plain", cube::Cubelet{D, T, I}) where {D, T, I}
    print(io, D, "-dim Cubelet", (index(cube), diam(cube)))
    print(io, ":\n  $(sign(cube) == 1 ? '+' : '-')$(vertices(cube))")
end

function Base.show(io::IO, cube::Cubelet{D, M}) where {D, M}
    print(io, "Cubelet{$D}($(sign(cube) == 1 ? '+' : '-')$(vertices(cube)), $(diam(cube)))")
end

index(cube::Cubelet) = cube.index
diam(cube::Cubelet) = cube.diam
coface_type(::Type{Cubelet{D, T, I}}) where {D, T, I} = Cubelet{D + 1, T, I}
face_type(::Type{Cubelet{D, T, I}}) where {D, T, I} = Cubelet{D - 1, T, I}

# @generated used because inference doesn't work well for 2^D
@generated function vertices(cube::Cubelet{D, <:Any, I}) where {D, I}
    return :(vertices(index(cube), Val($(2^D))))
end

Base.lastindex(::Cubelet{D}) where D = 2^D
Base.size(::Cubelet{D}) where D = (2^D,)

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
Base.CartesianIndices(cf::Cubical) = CartesianIndices(cf.data)
Base.LinearIndices(cf::Cubical) = LinearIndices(cf.data)

function simplex(
    cf::Cubical{I, T}, ::Val{D}, vertices, sign=one(I)
) where {I, T, D}
    diam = typemin(T)
    for v in vertices
        d = get(cf.data, v, nothing)
        if isnothing(d)
            return nothing
        else
            _d::T = d
            diam = max(diam, _d)
        end
    end
    return simplex_type(cf, D)(sign * index(vertices), diam)
end

function edges(cf::Cubical{I, <:Any}) where {I}
    E = edge_type(cf)
    result = E[]
    for u_lin in eachindex(cf.data)
        u_car = CartesianIndices(cf)[u_lin]
        for d in 1:dim(cf)
            v_car = u_car + CartesianIndex{dim(cf)}(ntuple(i -> i == d ? 1 : 0, dim(cf)))
            if v_car in CartesianIndices(cf)
                v_lin = LinearIndices(cf)[v_car]
                push!(result, E(
                    index((I(v_lin), I(u_lin))), max(cf.data[u_lin], cf.data[v_lin])
                ))
            end
        end
    end
    return sort!(result)
end

# coboundary ============================================================================= #
struct CubeletCoboundary{A, I, N, C<:Cubelet, F<:Cubical{I}, K}
    filtration::F
    cubelet::C
    vertices::SVector{K, CartesianIndex{N}}
end

function CubeletCoboundary{A}(
    filtration::F, cubelet::C
) where {A, I, D, F<:Cubical{I}, C<:Cubelet{D, <:Any, I}}
    K = length(cubelet)
    N = dim(filtration)
    vxs = map(v -> CartesianIndices(filtration)[v], vertices(cubelet))
    return CubeletCoboundary{A, I, N, C, F, K}(filtration, cubelet, vxs)
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

function Base.iterate(
    cc::CubeletCoboundary{A, I, N, C}, (dim, dir)=(one(I), one(I))
) where {A, I, N, C}
    while true
        while dim ≤ N && !all_equal_in_dim(dim, cc.vertices)
            dim += one(I)
        end
        dim > N && return nothing

        diff = CartesianIndex{N}(ntuple(isequal(dim), N) .* dir)
        new_vertices = cc.vertices .+ Ref(diff)
        #FIXME
        all_vertices = sort(map(v -> I(get(LinearIndices(cc.filtration), v, 0)),
                                vcat(cc.vertices, new_vertices)), rev=true)
        dim += I(dir == -1)
        dir *= -one(I)
        if A || all_vertices[end÷2+1:end] == LinearIndices(cc.filtration)[cc.vertices]

            sx = simplex(cc.filtration, Val(Ripserer.dim(C) + 1), all_vertices, -dir)
            if !isnothing(sx)
                _sx::simplex_type(cc.filtration, Ripserer.dim(C) + 1) = sx
                return _sx, (dim, dir)
            end
        end
    end
end

struct CubeletBoundary{D, I, C<:Cubelet, N, F<:Cubical{I}, K}
    filtration::F
    cubelet::C
    vertices::SVector{K, CartesianIndex{N}}
end

function CubeletBoundary(
    filtration::F, cubelet::C
) where {D, I, F<:Cubical{I}, C<:Cubelet{D}}
    K = length(cubelet)
    N = dim(filtration)
    vxs = map(v -> CartesianIndices(filtration)[v], vertices(cubelet))
    return CubeletBoundary{D, I, C, N, F, K}(filtration, cubelet, vxs)
end

boundary(filtration::Cubical, cubelet::Cubelet) = CubeletBoundary(filtration, cubelet)

function first_half(vec::SVector{K, T}) where {K, T}
    SVector{K÷2, T}(Tuple(vec)[1:K÷2])
end
function second_half(vec::SVector{K, T}) where {K, T}
    SVector{K÷2, T}(Tuple(vec)[K÷2+1:K])
end

# Idea: split the cube in each dimension, returning one of two halves depending on dir.
function Base.iterate(cb::CubeletBoundary{D, I, C}, (dim, dir)=(1, 1)) where {D, I, C}
    if dim > D
        return nothing
    else
        lin_vertices = map(v -> I(LinearIndices(cb.filtration)[v]),
                           sort(cb.vertices, by=v -> v[dim], rev=true))
        if dir == 1
            new_vertices = first_half(lin_vertices)
            sign = one(I)
            state = (dim, -dir)
        else
            new_vertices = second_half(lin_vertices)
            sign = -one(I)
            state = (dim + 1, 1)
        end
        sx = simplex(cb.filtration, Val(D - 1), sort(new_vertices, rev=true), sign)
        if !isnothing(sx)
            _sx::simplex_type(cb.filtration, D - 1) = sx
            return _sx, state
        end
    end
end
