using Documenter
using Ripserer
using SparseArrays

makedocs(sitename="Ripserer",
         pages=[
             "index.md",
             "api.md",
         ])

deploydocs(
    repo = "github.com/mtsch/Ripserer.jl.git",
)
