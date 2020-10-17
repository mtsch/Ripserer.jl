using Documenter
using Distances
using Ripserer
using Test

if Sys.iswindows()
    @warn "Skipping doctests due to problems with `Alpha` on Windows"
elseif VERSION ≥ v"1.6-DEV" || VERSION < v"1.5-DEV"
    @warn "Doctests were set up on Julia v1.5. Skipping."
else
    DocMeta.setdocmeta!(
        Ripserer, :DocTestSetup, :(using Ripserer; using Distances); recursive=true
    )
    doctest(Ripserer)
end
