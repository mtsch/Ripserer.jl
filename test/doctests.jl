using Documenter
using Distances
using Ripserer

if VERSION ≥ v"1.6-DEV" || VERSION < v"1.5-DEV" || !Sys.islinux()
    @warn "Doctests only run on Linux and on stable Julia"
else
    DocMeta.setdocmeta!(
        Ripserer, :DocTestSetup, :(using Ripserer; using Distances); recursive=true
    )
    doctest(Ripserer)
end
