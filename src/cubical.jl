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

function Base.show(io::IO, ::MIME"text/plain", csx::Cubelet{D, T, I}) where {D, T, I}
    print(io, D, "-dim Cubelet", (index(csx), diam(csx)))
    print(io, ":\n  $(sign(csx) == 1 ? '+' : '-')$(vertices(csx))")
end

function Base.show(io::IO, csx::Cubelet{D, M}) where {D, M}
    print(io, "Cubelet{$D}($(sign(csx) == 1 ? '+' : '-')$(vertices(csx)), $(diam(csx)))")
end

index(csx::Cubelet) = csx.index
diam(csx::Cubelet) = csx.diam
coface_type(::Type{Cubelet{D, T, I}}) where {D, T, I} = Cubelet{D + 1, T, I}
face_type(::Type{Cubelet{D, T, I}}) where {D, T, I} = Cubelet{D - 1, T, I}

# @generated used because inference doesn't work well for 2^D
@generated function vertices(csx::Cubelet{D, <:Any, I}) where {D, I}
    return :(vertices(index(csx), Val($(2^D))))
end

Base.lastindex(sx::Cubelet{D}) where D = 2^D
Base.size(sx::Cubelet{D}) where D = (2^D,)


"""
    Cubical{T, N} <: AbstractFiltration{T, <:Cubelet{0, T}}

`Cubical` is used to compute sublevel persistent homology on `N`-dimensional
images, which are of type `AbstractArray{T, N}`.

# Constructor

    Cubical(::AbstractArray{T, N})
"""
struct Cubical{I, T, N, A<:AbstractArray{T, N}} <: AbstractFiltration{T, Cubelet}
    data::A
    threshold::T
end

function Cubical(data)
    return Cubical{Int}(data)
end
function Cubical{I}(data::AbstractArray{T, N}) where {I, T, N}
    return Cubical{I, T, N, typeof(data)}(data, maximum(data))
end

n_vertices(cf::Cubical) = length(cf.data)
threshold(cf::Cubical) = cf.threshold
birth(cf::Cubical, i) = cf.data[i]
simplex_type(::Cubical{I, T}, dim) where {I, T} = Cubelet{dim, T, I}

Base.CartesianIndices(cf::Cubical) = CartesianIndices(cf.data)
Base.LinearIndices(cf::Cubical) = LinearIndices(cf.data)

function diam(cf::Cubical{<:Any, T}, vertices) where {T}
    res = typemin(T)
    for v in vertices
        if isnothing(get(cf.data, v, nothing))
            return missing
        else
            res = max(res, cf.data[v])
        end
    end
    return res
end

function edges(cf::Cubical{I, <:Any, N}) where {I, N}
    E = edge_type(cf)
    result = E[]
    for u_lin in eachindex(cf.data)
        u_car = CartesianIndices(cf)[u_lin]
        for dim in 1:N
            v_car = u_car + CartesianIndex{N}(ntuple(i -> i == dim ? 1 : 0, N))
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
struct CubeletCoboundary{A, N, I, C<:Cubelet, F<:Cubical{I, <:Any, N}, K}
    filtration::F
    cubelet::C
    vertices::SVector{K, CartesianIndex{N}}
end

function CubeletCoboundary{A}(
    filtration::F, cubelet::C
) where {A, N, I, D, F<:Cubical{I, <:Any, N}, C<:Cubelet{D, <:Any, I}}
    K = 2^D
    vxs = map(v -> CartesianIndices(filtration)[v], vertices(cubelet))
    return CubeletCoboundary{A, N, I, C, F, K}(filtration, cubelet, vxs)
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
    cc::CubeletCoboundary{A, N, I, C}, (dim, dir)=(one(I), one(I))
) where {A, N, I, C}
    # If not all indices in a given dimension are equal, we can't create a coface by
    # expanding in that direction.
    diameter = missing
    new_vertices = cc.vertices
    while ismissing(diameter)
        while dim ≤ N && !all_equal_in_dim(dim, cc.vertices)
            dim += one(I)
        end
        if dim > N
            break
        end

        diff = CartesianIndex{N}(ntuple(isequal(dim), N) .* dir)
        new_vertices = cc.vertices .+ Ref(diff)
        diameter = max(diam(cc.cubelet), diam(cc.filtration, new_vertices))

        dim += I(dir == -1)
        dir *= -one(I)
    end

    if ismissing(diameter)
        return nothing
    else
        all_vertices = sort(map(v -> I(LinearIndices(cc.filtration)[v]),
                                vcat(cc.vertices, new_vertices)), rev=true)
        return coface_type(C)(-dir * index(all_vertices), diameter), (dim, dir)
    end
end

struct CubeletBoundary{D, I, C<:Cubelet, N, F<:Cubical{I, <:Any, N}, K}
    filtration::F
    cubelet::C
    vertices::SVector{K, CartesianIndex{N}}
end

function CubeletBoundary(
    filtration::F, cubelet::C
) where {N, D, I, F<:Cubical{I, <:Any, N}, C<:Cubelet{D}}

    K = 2^D
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

# Idea: split the cube in each dimension, returning two halves depending on dir.
function Base.iterate(cb::CubeletBoundary{D, I, C}, (dim, dir)=(1, 1)) where {D, I, C}
    if dim > D
        return nothing
    else
        lin_vertices = map(v -> I(LinearIndices(cb.filtration)[v]),
                           sort(cb.vertices, by=v -> v[dim], rev=true))
        if dir == 1
            new_vertices = first_half(lin_vertices)
            diameter = diam(cb.filtration, new_vertices)
            return face_type(C)(new_vertices, diameter), (dim, -dir)
        else
            new_vertices = second_half(lin_vertices)
            diameter = diam(cb.filtration, new_vertices)
            return -face_type(C)(new_vertices, diameter), (dim + 1, 1)
        end
    end
end
