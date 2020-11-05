using SafeTestsets
using Test

###
### Basic stuff
###
@safetestset "primefield" begin
    include("primefield.jl")
end
@safetestset "simplex" begin
    include("simplex.jl")
end
@safetestset "chainelement" begin
    include("chainelement.jl")
end
@safetestset "chain" begin
    include("chain.jl")
end

###
### Computation structs
###
@safetestset "reducedmatrix" begin
    include("reducedmatrix.jl")
end
@safetestset "zerodimensional" begin
    include("zerodimensional.jl")
end

###
### Main tests
###
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

###
### Extra features
###
@safetestset "cycles" begin
    include("cycles.jl")
end
@safetestset "circularcoordinates" begin
    include("circularcoordinates.jl")
end
@safetestset "plotting" begin
    include("plotting.jl")
end

###
### Sanity checks
###
@safetestset "aqua" begin
    include("aqua.jl")
end
@safetestset "doctests" begin
    include("doctests.jl")
end
