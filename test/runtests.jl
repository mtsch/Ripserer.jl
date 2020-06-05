using SafeTestsets
using Test

include("test_helpers.jl")

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
@safetestset "reductionmatrix" begin
    include("reductionmatrix.jl")
end
@safetestset "reduction" begin
    include("reduction.jl")
end
@safetestset "plotting" begin
    include("plotting.jl")
end
