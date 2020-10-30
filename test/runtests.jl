using SafeTestsets
using Test

@safetestset "primefield" begin
    include("primefield.jl")
end
@safetestset "simplex" begin
    include("simplex.jl")
end
@safetestset "chain" begin
    include("chain.jl")
end
@safetestset "zerodimensional" begin
    include("zerodimensional.jl")
end
@safetestset "reductionmatrix" begin
    include("reductionmatrix.jl")
end

@safetestset "rips" begin
    include("filtrations/rips.jl")
end
@safetestset "cubical" begin
    include("filtrations/cubical.jl")
end
@safetestset "custom" begin
    include("filtrations/custom.jl")
end
@safetestset "alpha" begin
    include("filtrations/alpha.jl")
end
@safetestset "edgecollapse" begin
    include("filtrations/edgecollapse.jl")
end
@safetestset "new filtrations" begin
    include("filtrations/new-filtrations.jl")
end

@safetestset "cycles" begin
    include("cycles.jl")
end
@safetestset "plotting" begin
    include("plotting.jl")
end

@safetestset "aqua" begin
    include("aqua.jl")
end
@safetestset "doctests" begin
    include("doctests.jl")
end
