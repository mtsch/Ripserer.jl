using Documenter
using Distances
using Ripserer
using Test

DocMeta.setdocmeta!(
    Ripserer,
    :DocTestSetup,
    :(using Ripserer; using Distances);
    recursive=true,
)
doctest(Ripserer)
