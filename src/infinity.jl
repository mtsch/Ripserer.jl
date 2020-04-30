# infinity =============================================================================== #
"""
    Infinity

`Infinity()` is bigger than _anything_ else, except `missing` and `Inf`. It is used to:

* Avoiding using `typemax(T)` in persistence intervals. Getting death times of
  `9223372036854775807` doesn't look good.
* Returned by `diam(::AbstractFiltration, args...)` to signal that a simplex should be
  skipped.
"""
struct Infinity end

Base.show(io::IO, ::Infinity) =
    print(io, "∞")

(::Type{T})(::Infinity) where T<:AbstractFloat =
    typemax(T)
for op in (:<, :>, :isless, :isequal, :(==))
    @eval (Base.$op)(x::Real, ::Infinity) =
        $op(x, Inf)
    @eval (Base.$op)(::Infinity, x::Real) =
        $op(Inf, x)
end
Base.isapprox(::Infinity, x::Real; args...) =
    isapprox(Inf, x; args...)
Base.isapprox(x::Real, ::Infinity; args...) =
    isapprox(x, Inf; args...)

Base.isless(::Infinity, ::Missing) =
    false
Base.isless(::Missing, ::Infinity) =
    true
Base.isless(::Infinity, ::Infinity) =
    false
Base.:>(::Infinity, ::Infinity) =
    false
Base.isless(a, ::Infinity) =
    true
Base.isless(::Infinity, a) =
    false

Base.isfinite(::Infinity) =
    false

const ∞ = Infinity()
