using SafeTestsets
using Test

include("testhelpers.jl")

@safetestset "Aqua" begin
    include("aqua.jl")
end
@safetestset "simplex" begin
    include("simplex.jl")
end
@safetestset "ripsfiltration" begin
    include("ripsfiltration.jl")
end
@safetestset "cubical" begin
    include("cubical.jl")
end
@safetestset "primefield" begin
    include("primefield.jl")
end
@safetestset "chainelement" begin
    include("chainelement.jl")
end
@safetestset "zerodimensional" begin
    include("zerodimensional.jl")
end
@safetestset "reductionmatrix" begin
    include("reductionmatrix.jl")
end
@safetestset "ripserer" begin
    include("ripserer.jl")
end
@safetestset "plotting" begin
    include("plotting.jl")
end
