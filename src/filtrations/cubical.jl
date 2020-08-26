# We get cube diameters from a cubemap. The cubemap of a 3 × 3 image is an array with the
# following structure,
#
#      1 2 3 4 5
#    1 ⋅ — ⋅ — ⋅
#    2 | □ | □ |
#    3 ⋅ — ⋅ — ⋅
#    4 | □ | □ |
#    5 ⋅ — ⋅ — ⋅
#
# where `⋅` are the vertices, `—` and `|` are 1-cubes and `□` are 2-cubes. Each value in the
# resulting array is equal to the birth time of a cube, which is equal to the `maximum` of
# the birth times of cofaces.
#
# See https://link.springer.com/chapter/10.1007%2F978-3-642-23175-9_7 for more info.
@generated function cubemap(input::Array{T, N}) where {T, N}
    quote
        result = similar(input, size(input) .* 2 .- 1)
        result .= typemin(T)
        @inbounds @nloops $N i input begin
            index = CartesianIndex{$N}(@ntuple $N i)
            value = input[index]
            @nexprs $N j -> (range_j = max(2i_j - 2, 1):min(2i_j, size(result, j)))
            (@nref $N result range) .= max.((@nref $N result range), value)
        end
        return result
    end
end

one_hot(i, ::Val{N}) where N = CartesianIndex{N}(ntuple(isequal(i), Val(N)))

# Convenience functions from converting cubemap index to vertices and back.
@inline function from_cubemap(root::CartesianIndex{K}, ::Val{N}) where {K, N}
    2^count(iseven, Tuple(root)) == N || throw(ArgumentError("invalid N"))
    result = ntuple(_ -> root, Val(N))
    for (i, j) in enumerate(Tuple(root))
        if iseven(j)
            result = result .+ ntuple(Val(N)) do k
                ifelse(isodd(k), 1, -1) * one_hot(i, Val(K))
            end
        end
        result = TupleTools.sort(result, by=Base.Fix2(getindex, i))
    end
    result = SVector(TupleTools.sort(result))
    return map(result) do c
        CartesianIndex{K}((Tuple(c) .+ 1) .÷ 2)
    end
end

function to_cubemap(vertices::NTuple{N}) where N
    floor(log2(N)) == log2(N) || throw(ArgumentError("number vertices is not 2^d"))
    allunique(vertices) || throw(ArgumentError("duplicate vertices in $vertices"))
    K = length(first(vertices))
    for i in 1:K
        l, h = extrema(getindex.(vertices, i))
        h - l ≤ 1 || throw(ArgumentError("vertices not adjacent"))
    end
    result = ntuple(Val(K)) do i
        sum(2 .* getindex.(vertices, i) .- 1) ÷ N
    end
    return CartesianIndex{K}(result)
end

"""
    Cube{D, T, K} <: AbstractSimplex{D, T, CartesianIndex{K}}

A `Cube` is similar to a `Simplex`, but it has `2^D` vertices instead of `D+1`. The vertices
are encoded as the position in the CubeMap (see reference in [`Cubical`](@ref)). A `Cube`'s
vertices are of type `CartesianIndex{K}`.

# Example

```jldoctest
julia> Cube{1}(CartesianIndex(1, 2), 1.0)
1-dimensional Cube(index=CartesianIndex(1, 2), birth=1.0):
  +CartesianIndex{2}[CartesianIndex(1, 1), CartesianIndex(1, 2)]

```
"""
struct Cube{D, T, K} <: AbstractSimplex{D, T, CartesianIndex{K}}
    root::NTuple{K, Int32}
    birth::T

    function Cube{D, T, K}(root::CartesianIndex{K}, birth) where {D, T, K}
        D < 0 && throw(DomainError("dimension must be a non-negative integer"))
        return new{D, T, K}(Int32.(Tuple(root)), convert(T, birth))
    end
end
function Cube{D}(root::CartesianIndex{K}, birth::T) where {D, T, K}
    return Cube{D, T, K}(root, birth)
end
function Cube{D}(vertices, birth) where D
    root = to_cubemap(tuple(CartesianIndex.(vertices)...))
    return Cube{D}(root, birth)
end

birth(cube::Cube) = cube.birth
index(cube::Cube) = CartesianIndex(cube.root)
Base.sign(cube::Cube) = 1
Base.:-(cube::Cube) = cube
Base.abs(cube::Cube) = cube

@generated Base.length(::Type{<:Cube{D}}) where D = :($(2^D))
vertices(cube::Cube{D}) where D = from_cubemap(index(cube), Val(length(cube)))

"""
    Cubical{T, K} <: AbstractFiltration{CartesianIndex{K}, T}

`Cubical` is used to compute sublevel persistent homology on `N`-dimensional images, which
are of type `AbstractArray{T, N}`.

This type uses the CubeMap structure to find birth times of cubes (see reference).

# Constructor

* `Cubical(image::AbstractArray{T, N}, threshold=maximum(image))`

# Reference

Wagner, H., Chen, C., & Vuçini, E. (2012). [Efficient computation of persistent homology for
cubical data.](https://link.springer.com/chapter/10.1007/978-3-642-23175-9_7) In Topological
methods in data analysis and visualization II (pp. 91-106). Springer, Berlin, Heidelberg.

# Example

```jldoctest
julia> image = [1 0 0; 0 2 0; 0 0 0];

julia> ripserer(Cubical(image))[1]
1-element 0-dimensional PersistenceDiagram:
 [0.0, ∞)

julia> ripserer(Cubical(image))[2]
1-element 1-dimensional PersistenceDiagram:
 [1.0, 2.0)

```
"""
struct Cubical{K, T, A<:AbstractArray{T, K}} <: AbstractFiltration{CartesianIndex{K}, T}
    data::A
    cubemap::A
    threshold::T
end

function Cubical(data::AbstractArray{T, K}; threshold=maximum(data)) where {T, K}
    Cubical{K, T, typeof(data)}(data, cubemap(data), T(threshold))
end

nv(cf::Cubical) = length(cf.data)
threshold(cf::Cubical) = cf.threshold
simplex_type(::Type{<:Cubical{K, T}}, D) where {K, T} = Cube{D, T, K}

# when working with vertices, we don't cubemap them.
birth(cf::Cubical, i) = cf.data[i]
birth(cf::Cubical) = cf.data
vertices(cf::Cubical) = CartesianIndices(cf.data)

function edges(cf::Cubical{K}) where K
    result = edge_type(cf)[]
    for i in CartesianIndices(cf.cubemap)
        if count(iseven, Tuple(i)) == 1
            edge = unsafe_simplex(cf, Val(1), i, 1)
            if !isnothing(edge)
                push!(result, edge)
            end
        end
    end
    return result
end

function simplex(cf::Cubical{N, T}, ::Val{D}, vertices, sign=1) where {D, T, N}
    root = to_cubemap(vertices)
    return unsafe_simplex(cf, Val(D), root, sign)
end

function unsafe_simplex(cf::Cubical{N, T}, ::Val{D}, new_root, _) where {D, T, N}
    birth = get(cf.cubemap, new_root, missing)
    if ismissing(birth) || birth > cf.threshold
        return nothing
    else
        return Cube{D, T, N}(new_root, birth)
    end
end

struct CubeCoboundary{K, D, F<:Cubical{K}, C<:Cube{D}}
    filtration::F
    cube::C
end

function coboundary(filtration::Cubical{K}, cube::Cube{D}) where {K, D}
    return CubeCoboundary{K, D, typeof(filtration), typeof(cube)}(filtration, cube)
end

function Base.iterate(cob::CubeCoboundary{K, D}, (i,sign)=(K,1)) where {K, D}
    orig_root = index(cob.cube)
    while i ≤ K
        # To get coboundary in decreasing root order, we have to traverse i from K to 1 and
        # then back to K with sign reversed.
        next_i = i - sign
        next_sign = sign
        if next_i ≤ 0
            next_sign = -1
            next_i = 1
        end

        if isodd(orig_root[i])
            new_root = orig_root + sign * one_hot(i, Val(K))
            cofacet = unsafe_simplex(
                cob.filtration, Val(D + 1), new_root, sign
            )
            if !isnothing(cofacet)
                _cofacet::simplex_type(cob.filtration, D + 1) = cofacet
                return _cofacet, (next_i, next_sign)
            end
        end

        i, sign = next_i, next_sign
    end
    return nothing
end

struct CubeBoundary{K, D, F<:Cubical{K}, C<:Cube{D}}
    filtration::F
    cube::C
end

function boundary(filtration::Cubical{K}, cube::Cube{D}) where {K, D}
    return CubeBoundary{K, D, typeof(filtration), typeof(cube)}(filtration, cube)
end

function Base.iterate(bnd::CubeBoundary{K, D}, (i,sign)=(K,-1)) where {K, D}
    orig_root = index(bnd.cube)
    while i ≤ K
        next_i = i + sign
        next_sign = sign
        if next_i ≤ 0
            next_sign = 1
            next_i = 1
        end

        if iseven(orig_root[i])
            new_root = orig_root + sign * one_hot(i, Val(K))
            facet = unsafe_simplex(
                bnd.filtration, Val(D - 1), new_root, sign
            )
            if !isnothing(facet)
                _facet::simplex_type(bnd.filtration, D - 1) = facet
                return _facet, (next_i, next_sign)
            end
        end

        i, sign = next_i, next_sign
    end
    return nothing
end

struct CubicalColumnsToReduce{D, K, F<:Cubical{K}}
    filtration::F
end

function columns_to_reduce(filtration::Cubical{K}, itr) where K
    D = dim(eltype(itr)) + 1
    return CubicalColumnsToReduce{D, K, typeof(filtration)}(filtration)
end

function Base.iterate(cols::CubicalColumnsToReduce{D, K}, i=1) where {D, K}
    last_i = lastindex(cols.filtration.cubemap)
    cart = CartesianIndices(cols.filtration.cubemap)
    while i < last_i
        root = cart[i]
        if count(iseven, Tuple(root)) == D
            simplex = unsafe_simplex(cols.filtration, Val(D), root, 1)
            if !isnothing(simplex)
                _simplex::simplex_type(cols.filtration, D) = simplex
                return simplex, i + 1
            end
        end
        i += 1
    end
    return nothing
end

# Having a custom chain element type may save a few bytes per element. Allowing fields other
# than Mod{2} makes no sense and may present extra work with having to keep track of cube
# signs.
struct CubicalChainElement{C<:Cube} <: AbstractChainElement{C, Mod{2}}
    cube::C
    coefficient::Bool

    function CubicalChainElement{C}(cube, coef=Mod{2}(1)) where C<:Cube
        new{C}(cube, Mod{2}(coef) == Mod{2}(1))
    end
end

function coefficient(cce::CubicalChainElement)
    return ifelse(cce.coefficient, Mod{2}(1), Mod{2}(0))
end
function simplex(cce::CubicalChainElement)
    return cce.cube
end
function chain_element_type(::Type{C}, ::Type{Mod{2}}) where {C<:Cube}
    return CubicalChainElement{C}
end
function chain_element_type(::Type{C}, ::Type{F}) where {C<:Cube, F}
    error("only Mod{2} allowed as field type for cubical homology")
end
