# infinity =============================================================================== #
"""
    Infinity

`Infinity()` is a "strong" infinity and is bigger than _anything_ else, except `missing`. It
is used to:

* Avoiding using `typemax(T)` in persistence intervals. Getting death times of
  `9223372036854775807` doesn't look good. This also fixes the problem of calculating
  distances between infinite and finite intervals.
* Returned by `diam(::AbstractFiltration, args...)` to signal that a simplex should be
  skipped.
"""
  struct Infinity end

Base.show(io::IO, ::Infinity) = print(io, "∞")

Base.isless(::Infinity, a) = false
Base.isless(a, ::Infinity) = true
Base.:>(a, ::Infinity) = false
Base.:>(::Infinity, a) = true

Base.isless(::Infinity, ::Infinity) = false
Base.:>(::Infinity, ::Infinity) = false

Base.isless(::Infinity, ::Missing) = false
Base.isless(::Missing, ::Infinity) = true
Base.:>(::Missing, ::Infinity) = missing
Base.:>(::Infinity, ::Missing) = missing

Base.isapprox(::Infinity, x::Real; args...) = x == Inf
Base.isapprox(x::Real, ::Infinity; args...) = x == Inf

Base.isfinite(::Infinity) = false

const ∞ = Infinity()
