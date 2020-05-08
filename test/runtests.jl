using SafeTestsets
using Test

@testset "Ripserer" begin
    @safetestset "Aqua" begin
        include("aqua.jl")
    end
    @safetestset "infinity" begin
        include("infinity.jl")
    end
    @safetestset "simplices" begin
        include("simplex.jl")
    end
    @safetestset "ripsfiltration" begin
        include("ripsfiltration.jl")
    end
    @safetestset "primefield" begin
        include("primefield.jl")
    end
    @safetestset "diagram" begin
        include("diagram.jl")
    end
    @safetestset "chainelement" begin
        include("chainelement.jl")
    end
    @safetestset "reduction_structs" begin
        include("reduction_structs.jl")
    end
    @safetestset "reduction" begin
        include("reduction.jl")
    end
    @safetestset "plotting" begin
        include("plotting.jl")
    end
end
