using Compat
using Distances
using PersistenceDiagrams
using Ripserer
using StaticArrays
using Test

# Julia 1.0 does not allow these to be defined inside @testset.
struct NewRips <: Ripserer.AbstractRipsFiltration{Int, Float64} end
Ripserer.dist(::NewRips) = Float64[i ≠ j for i in 1:10, j in 1:10]
Ripserer.dist(::NewRips, i, j) = 1.0

@testset "New rips filtration" begin
    d0, d1 = ripserer(NewRips())
    @test d0 == [fill((0.0, 1.0), 9); (0.0, Inf)]
    @test d1 == []
end

struct NewFiltration <: Ripserer.AbstractFiltration{Int, Int} end

function Ripserer.unsafe_simplex(
    ::Type{Simplex{0, Int, Int}},
    ::NewFiltration,
    (v,),
    sign,
)
    return Simplex{0}(sign * v, 0)
end
function Ripserer.unsafe_simplex(
    ::Type{Simplex{D, Int, Int}},
    ::NewFiltration,
    vertices,
    sign,
) where D
    return Simplex{D}(sign * index(vertices), 1)
end
Ripserer.n_vertices(::NewFiltration) = 10
Ripserer.simplex_type(::Type{NewFiltration}, D) = Simplex{D, Int, Int}
Ripserer.edges(::NewFiltration) = Simplex{1}.(10:-1:1, 1)
function Ripserer.postprocess_diagram(::NewFiltration, diagram)
    result = PersistenceInterval[]
    for int in diagram
        if !hasproperty(int, :representative)
            push!(result, PersistenceInterval(birth(int) + 1, death(int) + 1))
        end
    end
    PersistenceDiagram(sort!(result); diagram.meta...)
end

@testset "New filtration" begin
    d0, d1 = ripserer(NewFiltration())
    @test d0 == [fill((1.0, 2.0), 4); fill((1.0, Inf), 6)]
    @test d1 == []

    d0, d1 = ripserer(NewFiltration(), reps=true)
    @test d0 == []
    @test d1 == []
end
