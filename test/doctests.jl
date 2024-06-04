using Documenter
using Distances
using Ripserer

if VERSION ≥ v"1.8-DEV" || VERSION < v"1.7-DEV" || !Sys.islinux()
    @warn "Doctests only run on Linux and Julia 1.7"
else
    DocMeta.setdocmeta!(
        Ripserer, :DocTestSetup, :(using Ripserer; using Distances); recursive=true
    )
    doctest(Ripserer)
end
