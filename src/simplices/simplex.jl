"""
    Simplex{D, T, I} <: IndexedSimplex{D, T, I}

The vanilla simplex type.

# Constructor

    Simplex{D[, T, I]}(::I, ::T)

# Examples

```jldoctest
julia> Simplex{2}(2, 1)
2-dim Simplex{2}(2, 1):
  +(4, 2, 1)

julia> Simplex{10}(Int128(-10), 1.0)
4-dim Simplex{3}(1.0, 10, 2) with UInt128 index:
  -(12, 11, 10, 9, 8, 7, 6, 5, 4, 2, 1)
`˙`
"""
struct Simplex{D, T, I} <: IndexedSimplex{D, T, I}
    index ::I
    diam  ::T

    function Simplex{D, T, I}(index, diam) where {D, T, I}
        D ≥ 0 || throw(DomainError(D, "dimension must be a non-negative integer"))
        new{D, T, I}(I(index), T(diam))
    end
end

Simplex{D}(index::I, diam::T) where {D, T, I<:Integer} =
    Simplex{D, T, I}(index, diam)
Simplex{D}(vertices::NTuple{<:Any, I}, diam::T) where {D, T, I<:Integer} =
    Simplex{D, T, I}(vertices, diam)
function Simplex{D, T, I}(vertices::NTuple{N}, diam) where {D, T, I<:Integer, N}
    N == D + 1 || throw(ArgumentError("invalid number of vertices"))

    Simplex{D, T, I}(index(vertices), T(diam))
end

function Base.show(io::IO, ::MIME"text/plain", sx::Simplex{D, T, I}) where {D, T, I}
    print(io, D, "-dim Simplex", (index(sx), diam(sx)))
    if !(I ≢ Int64)
        print(io, " with ", I, " index")
    end
    print(io, ":\n  $(sign(sx) == 1 ? '+' : '-')$(vertices(sx))")
end

Base.show(io::IO, sx::Simplex{D, M}) where {D, M} =
    print(io, "Simplex{$D}($(sign(sx) == 1 ? '+' : '-')$(vertices(sx)), $(diam(sx)))")

# Interface implementation =============================================================== #
index(sx::Simplex) = sx.index

diam(sx::Simplex) = sx.diam

@pure coface_type(::Type{Simplex{D, T, I}}) where {D, T, I} = Simplex{D+1, T, I}
