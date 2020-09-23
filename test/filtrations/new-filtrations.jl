using Compat
using Distances
using PersistenceDiagrams
using Ripserer
using SparseArrays
using StaticArrays
using Test

include(joinpath(@__DIR__, "test-datasets.jl"))
using Ripserer: nv, edges

# Julia 1.0 does not allow these to be defined inside @testset.

# This test is mostly supposed to test that pre-finding apparent pairs with
# `find_apparent_pairs` first works correctly. On the side, it also tests that rips
# filtrations that only overload `adjacency_matrix` and `threshold` work fine.
struct ApparentPairsRips{I, T, R<:Rips{I, T}} <: Ripserer.AbstractRipsFiltration{I, T}
    rips::R
    apparent::Ref{Int}
    normal::Ref{Int}
end

function ApparentPairsRips(data; kwargs...)
    return ApparentPairsRips(Rips(data; kwargs...), Ref(0), Ref(0))
end

Ripserer.adjacency_matrix(rw::ApparentPairsRips) = Ripserer.adjacency_matrix(rw.rips)
Ripserer.threshold(rw::ApparentPairsRips) = Ripserer.threshold(rw.rips)

# This works fine but is slower than the original algorithm (even if the pair finding code
# was optimized) as apparent pairs are found in `initialize_coboundary!` anyway. Doing this
# in parallel or on the GPU is an option, however.  See https://arxiv.org/abs/2003.07989
function Ripserer.find_apparent_pairs(rw::ApparentPairsRips, columns, _)
    S = eltype(columns)
    C = Ripserer.simplex_type(rw, dim(S)+1)
    cols_left = S[]
    apparent = Tuple{S, C}[]
    for σ in columns
        τ = minimum(Ripserer.coboundary(rw, σ)) # This is broken if coboundary is empty.
        σ′ = maximum(Ripserer.boundary(rw, τ))
        if σ′ == σ
            rw.apparent[] += 1
            push!(apparent, (σ, τ))
        else
            rw.normal[] += 1
            push!(cols_left, σ)
        end
    end
    return cols_left, apparent
end

@testset "Rips with `find_apparent_pairs` overloaded." begin
    for m in (2, 3, 5), t in (2, 1)
        @testset "projective plane with modulus=$m, threshold=$t" begin
            appa = ApparentPairsRips(projective_plane; threshold=t)
            res_rips = ripserer(Rips(projective_plane; threshold=t); modulus=m)
            res_appa = ripserer(appa; modulus=m)
            @test res_rips == res_appa
            # The nv - 1 edges were cleared by the twist algorithm.
            @test appa.apparent[] + appa.normal[] + nv(appa) - 1 == length(edges(appa))
        end
    end
    for m in (2, 3, 5), t in (8, 4)
        @testset "cycle with modulus=$m, threshold=$t" begin
            appa = ApparentPairsRips(cycle; threshold=t)
            res_rips = ripserer(Rips(cycle; threshold=t); modulus=m, dim_max=3)
            res_appa = ripserer(appa; modulus=m, dim_max=3)
            @test res_rips == res_appa
            @test appa.apparent[] > 0
        end
    end

    for m in (2, 3)
        @testset "Results are the same for all algorithms with modulus=$m" begin
            appa = ApparentPairsRips(projective_plane)
            _, hom_imp = ripserer(appa; alg=:homology, implicit=true, modulus=m)
            _, hom_exp = ripserer(appa; alg=:homology, implicit=false, modulus=m)
            _, coh_imp = ripserer(appa; alg=:cohomology, implicit=true, modulus=m)
            _, coh_exp = ripserer(appa; alg=:cohomology, implicit=false, modulus=m)

            @test hom_imp == hom_exp == coh_imp == coh_exp
            @test all(
                i.representative == j.representative for (i, j) in zip(hom_imp, hom_exp)
            )
        end
    end
end

# This test checks that apparent pairs can produce correct intervals. On the side, it checks
# `AbstractCustomFiltration`s.
struct ApparentPairsCustom <: Ripserer.AbstractCustomFiltration{Int, Int}
    found::Ref{Bool}
end
ApparentPairsCustom() = ApparentPairsCustom(Ref(false))

function Ripserer.simplex_dicts(::ApparentPairsCustom)
    return [
        Dict(1 => 0, 2 => 0, 3 => 0),
        Dict(1 => 1, 2 => 1, 3 => 2),
        Dict(1 => 3)
    ]
end
Ripserer.adjacency_matrix(::ApparentPairsCustom) = sparse(Bool[0 1 1; 1 0 1; 1 1 0])
function Ripserer.find_apparent_pairs(a::ApparentPairsCustom, columns, _)
    σ = columns[1]
    a.found[] = true
    return typeof(σ)[], [(σ, minimum(Ripserer.coboundary(a, σ)))]
end

@testset "Apparent pairs can produce valid intervals." begin
    @testset "cohomology" begin
        appa = ApparentPairsCustom()
        @test !appa.found[]
        d0, d1 = ripserer(appa, reps=true)
        @test appa.found[]
        @test d0 == [(0, 1), (0, 1), (0, Inf)]
        @test d1 == [(2, 3)]
        @test d1[1].birth_simplex == Simplex{1}(3, 2)
        @test d1[1].death_simplex == Simplex{2}(1, 3)
        @test simplex(only(d1[1].representative)) == Simplex{1}(3, 2)
    end
    @testset "Results are the same for all algorithms" begin
        appa = ApparentPairsCustom()
        _, hom_imp = ripserer(appa; alg=:homology, implicit=true)
        _, hom_exp = ripserer(appa; alg=:homology, implicit=false)
        _, coh_imp = ripserer(appa; alg=:cohomology, implicit=true)
        _, coh_exp = ripserer(appa; alg=:cohomology, implicit=false)

        @test hom_imp == hom_exp == coh_imp == coh_exp
        @test all(i.representative == j.representative for (i, j) in zip(hom_imp, hom_exp))
    end
end

# The main idea of this test is to test postprocess_diagrams
struct NoRepsFiltration <: Ripserer.AbstractFiltration{Int, Int} end

function Ripserer.unsafe_simplex(
    ::Type{Simplex{0, Int, Int}},
    ::NoRepsFiltration,
    (v,),
    sign,
)
    return Simplex{0}(sign * v, 0)
end
function Ripserer.unsafe_simplex(
    ::Type{Simplex{D, Int, Int}},
    ::NoRepsFiltration,
    vertices,
    sign,
) where D
    return Simplex{D}(sign * index(vertices), 1)
end
Ripserer.nv(::NoRepsFiltration) = 10
Ripserer.simplex_type(::Type{NoRepsFiltration}, D) = Simplex{D, Int, Int}
Ripserer.edges(::NoRepsFiltration) = Simplex{1}.(10:-1:1, 1)
function Ripserer.postprocess_diagram(::NoRepsFiltration, diagram)
    result = PersistenceInterval[]
    for int in diagram
        if !hasproperty(int, :representative)
            push!(result, PersistenceInterval(birth(int) + 1, death(int) + 1))
        end
    end
    PersistenceDiagram(sort!(result); diagram.meta...)
end

@testset "New filtration" begin
    d0, d1 = ripserer(NoRepsFiltration())
    @test d0 == [fill((1.0, 2.0), 4); fill((1.0, Inf), 6)]
    @test d1 == []

    d0, d1 = ripserer(NoRepsFiltration(), reps=true)
    @test d0 == []
    @test d1 == []
end
