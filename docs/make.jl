@info "build started."
using Documenter
using Literate
using Ripserer
using PersistenceDiagrams
using Plots
gr()
ENV["GKSwstype"] = "100"

EXAMPLES_INPUT = joinpath(@__DIR__, "src/examples")
EXAMPLES_OUTPUT = joinpath(@__DIR__, "src/generated")

for example in readdir(EXAMPLES_INPUT; join=true)
    endswith(example, ".jl") || continue
    Literate.markdown(example, EXAMPLES_OUTPUT; documenter=true)
end

makedocs(;
    sitename="Ripserer.jl",
    format=Documenter.HTML(
        # Use clean URLs, unless built as a "local" build;
        ;
        prettyurls=!("local" in ARGS),
        assets=["assets/favicon.ico"],
    ),
    pages=[
        "Home" => "index.md",
        "Examples" => [
            "generated/basics.md",
            "generated/stability.md",
            "generated/cocycles.md",
            "generated/cubical.md",
            "generated/malaria.md",
        ],
        "Benchmarks" => "benchmarks.md",
        "API" => ["Public" => "api/ripserer.md", "Interfaces" => "api/extensions.md"],
        "Related Julia Packages" => "related-work.md",
        "Acknowledgements and References" => "references.md",
    ],
    doctest=false, # Doctests are run as part of testing -- no need to run them twice.
)

deploydocs(; repo="github.com/mtsch/Ripserer.jl.git")
