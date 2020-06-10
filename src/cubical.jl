"""
    Cubelet{D, T, I} <: IndexedSimplex{D, T, I}

A `Cubelet` is similar to a `Simplex`, but it has `2^D` vertices instead of `D+1`. Like in
`Simplex`, the vertices are encoded from an index and dimension. Because a cubelet knows
nothing about the image it came from, it returns *linear* indices from `vertices`.

The vertices should be neighboring indices, but this fact is not checked anywhere.
"""
struct Cubelet{D, T, I} <: IndexedSimplex{D, T, I}
    index ::I
    diam  ::T

    function Cubelet{D, T, I}(index, diam) where {D, T, I}
        D ≥ 0 || throw(DomainError(D, "dimension must be a non-negative integer"))
        return new{D, T, I}(I(index), T(diam))
    end
end
Cubelet{D}(index::I, diam::T) where {D, T, I} = Cubelet{D, T, I}(index, diam)

function Cubelet{D}(vertices::NTuple{K, I}, diam::T) where {D, T, K, I}
    K == 2^D || throw(ArgumentError("a `Cubelet` must have 2^D simplices"))
    return Cubelet{D, T, I}(index(vertices), diam)
end
function Cubelet{D, T, I}(vertices::NTuple, diam) where {D, T, I}
    return Cubelet{D, T, I}(index(vertices), diam)
end

function Base.show(io::IO, ::MIME"text/plain", csx::Cubelet{D, T, I}) where {D, T, I}
    print(io, D, "-dim Cubelet", (index(csx), diam(csx)))
    if I ≢ Int64
        print(io, " with ", I, " index")
    end
    print(io, ":\n  $(sign(csx) == 1 ? '+' : '-')$(vertices(csx))")
end

function Base.show(io::IO, csx::Cubelet{D, M}) where {D, M}
    print(io, "Cubelet{$D}($(sign(csx) == 1 ? '+' : '-')$(vertices(csx)), $(diam(csx)))")
end

diam(csx::Cubelet) = csx.diam
index(csx::Cubelet) = csx.index

# @generated used because inference doesn't work well for 2^D
@generated function vertices(csx::Cubelet{D, <:Any, I}) where {D, I}
    return :(vertices(index(csx), Val($(2^D))))
end

coface_type(::Type{Cubelet{D, T, I}}) where {D, T, I} = Cubelet{D + 1, T, I}
face_type(::Type{Cubelet{D, T, I}}) where {D, T, I} = Cubelet{D - 1, T, I}

"""
    Cubical{T, N} <: AbstractFiltration{T, <:Cubelet{0, T}}

`Cubical` is used to compute sublevel persistent homology on `N`-dimensional
images, which are of type `AbstractArray{T, N}`.

# Constructor

    Cubical(::AbstractArray{T, N})
"""
struct Cubical{
    T, N, V<:Cubelet{0, T}, A<:AbstractArray{T, N}
} <: AbstractFiltration{T, V}
    data      ::A
    threshold ::T
end

function Cubical(
    data::AbstractArray{T, N};
    vertex_type::DataType=Cubelet{0, T, Int64}
) where {T, N}
    return Cubical{T, N, vertex_type, typeof(data)}(data, maximum(data))
end

n_vertices(cf::Cubical) = length(cf.data)
threshold(cf::Cubical) = cf.threshold
birth(cf::Cubical, i) = cf.data[i]

Base.CartesianIndices(cf::Cubical) = CartesianIndices(cf.data)
Base.LinearIndices(cf::Cubical) = LinearIndices(cf.data)

# doesn't quite follow interface.
function diam(cf::Cubical{T}, vertices) where {T}
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

function edges(cf::Cubical{<:Any, N}) where N
    E = edge_type(cf)
    result = E[]
    for u_lin in eachindex(cf.data)
        u_car = CartesianIndices(cf)[u_lin]
        for dim in 1:N
            v_car = u_car + CartesianIndex(ntuple(i -> i == dim ? 1 : 0, N))
            if v_car in CartesianIndices(cf)
                v_lin = LinearIndices(cf)[v_car]
                push!(result, E(index((v_lin, u_lin)), max(cf.data[u_lin], cf.data[v_lin])))
            end
        end
    end
    return sort!(result)
end

# coboundary ============================================================================= #
struct CubeletCoboundary{A, N, C<:Cubelet, F<:Cubical{<:Any, N}, K}
    filtration ::F
    cubelet    ::C
    vertices   ::NTuple{K, CartesianIndex{N}}
end

function CubeletCoboundary{A}(
    filtration::F, cubelet::C
) where {A, N, D, F<:Cubical{<:Any, N}, C<:Cubelet{D}}

    K = 2^D
    vxs = map(v -> CartesianIndices(filtration)[v], vertices(cubelet))
    return CubeletCoboundary{A, N, C, F, K}(filtration, cubelet, vxs)
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

function Base.iterate(cc::CubeletCoboundary{A, N, C}, (dim, dir)=(1, 1)) where {A, N, C}
    # If not all indices in a given dimension are equal, we can't create a coface by
    # expanding in that direction.
    diameter = missing
    new_vertices = cc.vertices
    while ismissing(diameter)
        while dim ≤ N && !all_equal_in_dim(dim, cc.vertices)
            dim += 1
        end
        if dim > N
            break
        end

        diff = CartesianIndex{N}(ntuple(isequal(dim), N) .* dir)
        new_vertices = cc.vertices .+ Ref(diff)
        diameter = max(diam(cc.cubelet), diam(cc.filtration, new_vertices))

        dim += Int(dir == -1)
        dir *= -1
    end

    if ismissing(diameter)
        return nothing
    # We swapped the direction of dir at the end of the loop so we use -dir everywhere.
    elseif dir == -1
        all_vertices = map(v -> LinearIndices(cc.filtration)[v],
                           TupleTools.vcat(new_vertices, cc.vertices))
    else
        all_vertices = map(v -> LinearIndices(cc.filtration)[v],
                           TupleTools.vcat(cc.vertices, new_vertices))
    end
    all_vertices = TupleTools.sort(all_vertices, rev=true)
    return coface_type(C)(-dir * index(all_vertices), diameter), (dim, dir)
end

struct CubeletBoundary{D, C<:Cubelet, N, F<:Cubical{<:Any, N}, K}
    filtration ::F
    cubelet    ::C
    vertices   ::NTuple{K, CartesianIndex{N}}
end

function CubeletBoundary(
    filtration::F, cubelet::C
) where {N, D, F<:Cubical{<:Any, N}, C<:Cubelet{D}}

    K = 2^D
    vxs = map(v -> CartesianIndices(filtration)[v], vertices(cubelet))
    return CubeletBoundary{D, C, N, F, K}(filtration, cubelet, vxs)
end

boundary(filtration::Cubical, cubelet::Cubelet) = CubeletBoundary(filtration, cubelet)

# Idea: split the cube in each dimension, returning two halves depending on dir.
function Base.iterate(cb::CubeletBoundary{D, C}, (dim, dir)=(1, 1)) where {D, C}
    if dim > D
        return nothing
    else
        lin_vertices = map(v -> LinearIndices(cb.filtration)[v],
                           TupleTools.sort(cb.vertices, by=v -> v[dim], rev=true))
        if dir == 1
            new_vertices = lin_vertices[1:end÷2]
            diameter = diam(cb.filtration, new_vertices)
            return face_type(C)(index(new_vertices), diameter), (dim, -dir)
        else
            new_vertices = lin_vertices[end÷2+1:end]
            diameter = diam(cb.filtration, new_vertices)
            return face_type(C)(-index(new_vertices), diameter), (dim + 1, 1)
        end
    end
end
