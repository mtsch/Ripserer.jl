using Documenter
using Distances
using Ripserer

if VERSION ≥ v"1.11-DEV" || VERSION < v"1.10-DEV" || !Sys.islinux()
    @warn "Doctests only run on Linux and Julia 1.10"
else
    DocMeta.setdocmeta!(
        Ripserer, :DocTestSetup, :(using Ripserer; using Distances); recursive=true
    )
    doctest(Ripserer; fix=true)
end
