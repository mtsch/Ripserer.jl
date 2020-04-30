using SafeTestsets
using Test

@testset "Ripserer" begin
    @safetestset "Aqua" begin
        include("aqua.jl")
    end
    @safetestset "infinity" begin
        include("infinity.jl")
    end
    @safetestset "simplex" begin
        include("simplex.jl")
    end
    @safetestset "ripsfiltration" begin
        include("ripsfiltration.jl")
    end
    @safetestset "diagram" begin
        include("diagram.jl")
    end
    @safetestset "reduction" begin
        include("reduction.jl")
    end
end
