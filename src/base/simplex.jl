"""
    Simplex{D, T, I<:Integer} <: AbstractSimplex{D, T, I}

The vanilla simplex type represented by dimension `D`, an index of type `I`, and a birth
time of type `T`.

# Constructors

* `Simplex{D[, T, I]}(index, birth)`
* `Simplex{D}(vertices, birth)`: vertices must be sorted descending. This constructor mainly
  exists for debugging purposes. Using [`simplex`](@ref) is usually the better option.

# Examples

```jldoctest
julia> Simplex{2}(2, 1)
2-dimensional Simplex(index=2, birth=1):
  +(4, 2, 1)

julia> Simplex{10}(Int128(-10), 1.0)
10-dimensional Simplex(index=10, birth=1.0):
  -(12, 11, 10, 9, 8, 7, 6, 5, 4, 2, 1)

julia> Simplex{2}((5, 2, 1), 1)
2-dimensional Simplex(index=5, birth=1):
  +(5, 2, 1)

```
"""
struct Simplex{D,T,I<:Integer} <: AbstractSimplex{D,T,I}
    index::I
    birth::T

    function Simplex{D,T,I}(index::Integer, birth) where {D,T,I<:Integer}
        D ≥ 0 || throw(DomainError(D, "dimension must be a non-negative integer"))
        return new{D,T,I}(I(index), convert(T, birth))
    end
end

function Simplex{D}(index::I, birth::T) where {D,T,I<:Integer}
    return Simplex{D,T,I}(index, birth)
end
function Simplex{D}(vertices, birth::T) where {D,T}
    length(vertices) == D + 1 ||
        throw(ArgumentError("invalid number of vertices $(length(vertices))"))
    vertices_svec = sort(SVector{D + 1}(vertices); rev=true)
    return Simplex{D,T,eltype(vertices)}(index(vertices_svec), T(birth))
end

index(sx::Simplex) = abs(sx.index)
birth(sx::Simplex) = sx.birth
Base.sign(sx::Simplex) = sign(sx.index)
Base.:-(sx::S) where {S<:Simplex} = S(-sx.index, birth(sx))
