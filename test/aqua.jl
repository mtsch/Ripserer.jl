using Ripserer
using Aqua

# Ignore ambiguities in StaticArrays.
Aqua.test_all(
    Ripserer,
    ambiguities=(exclude=[Base.convert, Base.unsafe_convert, Base.all],)

)
