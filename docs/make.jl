using Documenter
using Ripserer
using SparseArrays

makedocs(sitename="Ripserer.jl",
         format = Documenter.HTML(
             # Use clean URLs, unless built as a "local" build
             prettyurls = !("local" in ARGS),
             assets = ["assets/favicon.ico"],
         ),
         pages=[
             "Home" => "index.md",
             "Quick Start" => "quickstart.md",
             "API" => "api.md",
         ])

deploydocs(
    repo = "github.com/mtsch/Ripserer.jl.git",
)
