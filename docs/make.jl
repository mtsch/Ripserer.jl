@info "build started."
using Documenter
using Literate
using Ripserer
using Plots

for example in readdir(joinpath(@__DIR__, "src/literate"), join=true)
    endswith(example, "jl") || continue
    Literate.markdown(example, "src/generated", documenter=true)
end

makedocs(sitename="Ripserer.jl",
         format = Documenter.HTML(
             # Use clean URLs, unless built as a "local" build
             prettyurls = !("local" in ARGS),
             assets = ["assets/favicon.ico"],
         ),
         pages=[
             "Home" => "index.md",
             "Examples" => [
                 "generated/basics.md",
                 "generated/stability.md",
                 "generated/cocycles.md",
                 "generated/time_series_sublevel.md",
             ],
             "API" => "api.md",
         ])

deploydocs(
    repo = "github.com/mtsch/Ripserer.jl.git",
)
