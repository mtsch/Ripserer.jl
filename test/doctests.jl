using Documenter
using Distances
using Ripserer
using Test

if v"1.6-DEV" > VERSION ≥ v"1.5-DEV"
    DocMeta.setdocmeta!(
        Ripserer,
        :DocTestSetup,
        :(using Ripserer; using Distances);
        recursive=true,
    )
    doctest(Ripserer)
else
    @warn "Doctests were set up on Julia v1.5. Skipping."
end
