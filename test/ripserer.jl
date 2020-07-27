using Compat
using Distances
using PersistenceDiagrams
using Ripserer
using StaticArrays
using Suppressor
using Test

using Ripserer: zeroth_intervals, ChainElement, PackedElement


include("data.jl")

@testset "Rips" begin
    @testset "No threshold" begin
        @testset "Icosahedron" begin
            d0, d1, d2 = ripserer(icosahedron; dim_max=2)
            @test d0 == [fill((0.0, 1.0), 11); (0.0, Inf)]
            @test d1 == []
            @test d2 == [(1.0, 2.0)]
        end
        @testset "Cycle with various fields" begin
            d0_2, d1_2, d2_2, d3_2 = ripserer(Rips{Int32}(cycle); dim_max=3)
            d0_7, d1_7, d2_7, d3_7 = ripserer(Rips(cycle); dim_max=3, modulus=7)
            d0_r, d1_r, d2_r, d3_r = ripserer(cycle; dim_max=3, field_type=Rational{Int})

            @test d0_2 == d0_7 == d0_r == [fill((0, 1), size(cycle, 1) - 1); (0, Inf)]
            @test d1_2 == d1_7 == d1_r == [(1, 6)]
            @test d2_2 == d2_7 == d2_r == fill((6, 7), 5)
            @test d3_2 == d3_7 == d3_r == [(7, 8)]
        end
        @testset "RP2 with various fields" begin
            _, d1_2, d2_2 = ripserer(projective_plane; dim_max=2)
            _, d1_3, d2_3 = ripserer(projective_plane; dim_max=2, modulus=3)
            _, d1_331, d2_331 = ripserer(projective_plane; dim_max=2, field_type=Mod{5})
            _, d1_r, d2_r = ripserer(projective_plane; dim_max=2, field_type=Rational{Int})
            @test d1_2 == [(1, 2)]
            @test d2_2 == [(1, 2)]
            @test d1_3 == d1_331 == d1_r == []
            @test d2_3 == d2_331 == d1_r == []
        end
    end
    @testset "Threshold" begin
        @testset "Icosahedron, threshold=1" begin
            d0, d1, d2 = ripserer(Rips(icosahedron; threshold=1); dim_max=2)
            @test d0 == [fill((0.0, 1.0), 11); (0.0, Inf)]
            @test d1 == []
            @test d2 == [(1.0, Inf)]
        end
        @testset "Icosahedron, threshold=0.5" begin
            d0, d1, d2 = ripserer(icosahedron; dim_max=2, threshold=0.5)
            @test d0 == fill((0.0, Inf), 12)
            @test d1 == []
            @test d2 == []
        end
        @testset "RP2 with various fields, threshold=1" begin
            _, d1_2, d2_2 = ripserer(projective_plane;
                                     dim_max=2, threshold=1)
            _, d1_3, d2_3 = ripserer(projective_plane;
                                     dim_max=2, modulus=3, threshold=1)
            _, d1_331, d2_331 = ripserer(projective_plane;
                                         dim_max=2, field_type=Mod{5}, threshold=1)
            _, d1_r, d2_r = ripserer(projective_plane;
                                     dim_max=2, field_type=Rational{Int}, threshold=1)
            @test d1_2 == [(1, Inf)]
            @test d2_2 == [(1, Inf)]
            @test d1_3 == d1_331 == d1_r == []
            @test d2_3 == d2_331 == d1_r == []
        end
    end
    @testset "Points as input" begin
        for metric in (Euclidean(), Cityblock())
            pts = torus_points(9)
            @test ripserer(pts; metric=metric) == ripserer(Rips(pts; metric=metric))
        end
    end
    @testset "Cutoff" begin
        d0, d1 = ripserer(rand_dist_matrix(20), cutoff=0.5)
        @test all(persistence.(d0) .> 0.5)
        @test all(persistence.(d1) .> 0.5)
    end
end

@testset "SparseRips" begin
    @testset "Icosahedron" begin
        d0, d1, d2 = ripserer(sparse(icosahedron); dim_max=2)
        @test d0 == [fill((0.0, 1.0), 11); (0.0, Inf)]
        @test d1 == []
        @test d2 == [(1.0, 2.0)]
    end
    @testset "RP2 with various fields, threshold=1" begin
        _, d1_2, d2_2 = ripserer(sparse(projective_plane);
                                 dim_max=2, threshold=1)
        _, d1_3, d2_3 = ripserer(SparseRips(projective_plane; threshold=1);
                                 dim_max=2, modulus=3)
        _, d1_331, d2_331 = ripserer(sparse(projective_plane);
                                     dim_max=2, field_type=Mod{5}, threshold=1)
        _, d1_r, d2_r = ripserer(sparse(projective_plane);
                                 dim_max=2, field_type=Rational{Int}, threshold=1)
        @test d1_2 == [(1, Inf)]
        @test d2_2 == [(1, Inf)]
        @test d1_3 == d1_331 == d1_r == []
        @test d2_3 == d2_331 == d1_r == []
    end
    @testset "Equal to Rips" begin
        for thresh in (nothing, 1, 0.5, 0.126)
            data = rand_torus(100)
            r_res = ripserer(data, threshold=thresh, dim_max=2)
            s_res_1 = ripserer(SparseRips(data, threshold=thresh), dim_max=2)

            # Add zeros to diagonal. Adding ones first actually changes the structure of
            # the matrix.
            data2 = sparse(data)
            for i in axes(data2, 1)
                data2[i, i] = 1
                data2[i, i] = 0
            end
            s_res_2 = ripserer(data2, threshold=thresh, dim_max=2)

            @test r_res == s_res_1 == s_res_2
        end
    end
end

@testset "Representatives" begin
    @testset "Known example" begin
        # This example was generated by getting representatives from ripser.
        _, d1, d2 = ripserer(projective_plane; dim_max=2, reps=true)

        @test simplex.(representative(only(d1))) == [
            Simplex{1}((11, 10), 1),
            Simplex{1}((10, 7), 1),
            Simplex{1}((10, 6), 1),
            Simplex{1}((8, 1), 1),
            Simplex{1}((7, 3), 1),
            Simplex{1}((7, 1), 1),
            Simplex{1}((6, 2), 1),
            Simplex{1}((5, 1), 1),
            Simplex{1}((2, 1), 1),
        ]
        @test coefficient.(representative(only(d1))) == fill(Mod{2}(1), 9)
        @test simplex.(representative(only(d2))) == [Simplex{2}((6, 2, 1), 1)]
        @test coefficient.(representative(only(d2))) == [Mod{2}(1)]
    end
    @testset "Types" begin
        d0, d1, d2, d3 = ripserer(cycle; dim_max=3, reps=true)
        @test eltype(d0) <: RepresentativeInterval{
            PersistenceInterval,
            Simplex{0, Int, Int},
            Union{Nothing, Simplex{1, Int, Int}},
            <:Vector{<:PackedElement{Simplex{0, Int, Int}, Mod{2}}}}
        @test eltype(d1) <: RepresentativeInterval{
            PersistenceInterval,
            Simplex{1, Int, Int},
            Union{Nothing, Simplex{2, Int, Int}},
            <:Vector{<:PackedElement{Simplex{1, Int, Int}, Mod{2}}}}
        @test eltype(d2) <: RepresentativeInterval{
            PersistenceInterval,
            Simplex{2, Int, Int},
            Union{Nothing, Simplex{3, Int, Int}},
            <:Vector{<:PackedElement{Simplex{2, Int, Int}, Mod{2}}}}
        @test eltype(d3) <: RepresentativeInterval{
            PersistenceInterval,
            Simplex{3, Int, Int},
            Union{Nothing, Simplex{4, Int, Int}},
            <:Vector{<:PackedElement{Simplex{3, Int, Int}, Mod{2}}}}

        d0, d1, d2, d3 = ripserer(cycle; dim_max=3, reps=true,
                                  field_type=Rational{Int})
        @test eltype(d0) ≡ RepresentativeInterval{
            PersistenceInterval,
            Simplex{0, Int, Int},
            Union{Nothing, Simplex{1, Int, Int}},
            Vector{ChainElement{Simplex{0, Int, Int}, Rational{Int}}}}
        @test eltype(d1) ≡ RepresentativeInterval{
            PersistenceInterval,
            Simplex{1, Int, Int},
            Union{Nothing, Simplex{2, Int, Int}},
            Vector{ChainElement{Simplex{1, Int, Int}, Rational{Int}}}}
        @test eltype(d2) ≡ RepresentativeInterval{
            PersistenceInterval,
            Simplex{2, Int, Int},
            Union{Nothing, Simplex{3, Int, Int}},
            Vector{ChainElement{Simplex{2, Int, Int}, Rational{Int}}}}
        @test eltype(d3) ≡ RepresentativeInterval{
            PersistenceInterval,
            Simplex{3, Int, Int},
            Union{Nothing, Simplex{4, Int, Int}},
            Vector{ChainElement{Simplex{3, Int, Int}, Rational{Int}}}}
    end
    @testset "Infinite interval" begin
        _, d1 = ripserer(cycle; dim_max=1, reps=true, threshold=1, field_type=Rational{Int})
        rep = representative(only(d1))
        @test rep == []
        @test eltype(rep) == ChainElement{Simplex{1, Int, Int}, Rational{Int}}
    end
    @testset "Critical simplices" begin
        result = ripserer(torus(100); reps=true, threshold=0.5)
        for diag in result
            @test birth.(diag) == birth.(birth_simplex.(diag))
            finite = filter(isfinite, diag)
            @test death.(finite) == birth.(death_simplex.(finite))
            infinite = filter(!isfinite, diag)
            @test all(isnothing, death_simplex.(infinite))
        end
    end
end

@testset "Zero-dimensional sublevel set persistence" begin
    @testset "with SparseRips" begin
        data = [1, 0, 1, 2, 3, 4, 3, 2, 3, 2, 1, 2]

        # Create distance matrix from data, where neighboring points are connected by edges
        # and the edge weights are equal to the max of both vertex births.
        n = length(data)
        dists = spzeros(n, n)
        for i in 1:n
            dists[i, i] = data[i]
        end
        for i in 1:n-1
            j = i + 1
            dists[i, j] = dists[j, i] = max(dists[i, i], dists[j, j])
        end
        # 0-dimensional persistence should find values of minima and maxima of our data.
        d0 = ripserer(dists; dim_max=0)[1]
        @test d0 == [(0, Inf), (1, 4), (2, 3)]
    end

    @testset "with Rips" begin
        data = [1, 0, 1, 2, 3, 4, 3, 2, 3, 2, 1, 2]

        n = length(data)
        dists = zeros(Int, (n, n))
        for i in 1:n, j in 1:n
            if abs(i - j) ≤ 1
                dists[i, j] = max(data[i], data[j])
            else
                dists[i, j] = 5
            end
        end
        d0 = ripserer(dists; dim_max=0)[1]
        @test d0 == [(0, Inf), (1, 4), (2, 3)]
    end
end

@testset "Cubical" begin
    @testset "1D curve" begin
        data = [1, 0, 1, 2, 3, 4, 3, 2, 3, 2, 1, 2]
        d0, d1 = ripserer(Cubical(data); dim_max=2)

        @test d0 == [(0, Inf), (1, 4), (2, 3)]
        @test d1 == []
    end
    @testset "2D image" begin
        data = [0 0 0 0 0;
                0 2 2 2 0;
                0 2 1 2 0;
                0 2 2 2 0;
                0 0 0 0 0]

        d0, d1, d2 = ripserer(Cubical(data); reps=true, dim_max=2)

        @test d0 == [(0, Inf), (1, 2)]
        @test d1 == [(0, 2)]
        @test d2 == []

        @test sort(vertices.(representative(d0[1]))) ==
            sort(SVector.(vec(CartesianIndices(data))))
        @test vertices(only(representative(d0[2]))) ==
            SVector(CartesianIndex(3, 3))
    end
    @testset "3D image" begin
        # Cube with hole in the middle.
        data = zeros(5, 5, 5)
        data[2, 2:4, 2:4] .= 1
        data[3, :, :] .= [0 0 0 0 0; 0 1 1 1 0; 0 1 0 1 0; 0 1 1 1 0; 0 0 0 0 0]
        data[4, 2:4, 2:4] .= 1

        d0, d1, d2 = ripserer(Cubical(data); dim_max=2)

        @test d0 == [(0, 1.0), (0, Inf)]
        @test isempty(d1)
        @test d2 == [(0, 1)]
    end
end

@testset "Persistent homology" begin
    @testset "Produces the same diagram as cohomology" begin
        res_hom = ripserer(cycle, cohomology=false, dim_max=3)
        res_coh = ripserer(cycle, cohomology=true, dim_max=3)

        @test res_hom == res_coh
    end
    @testset "Has the same cirical simplices as cohomology" begin
        # Add some noise because critical simplices might be different if values are exactly
        # the same.
        cyc = cycle .+ 0.01 .* rand_dist_matrix(18)
        res_hom = ripserer(cyc, cohomology=false, dim_max=3, reps=true)
        res_coh = ripserer(cyc, cohomology=true, dim_max=3, reps=true)

        for i in 1:4
            @test birth_simplex.(res_hom[i]) == birth_simplex.(res_coh[i])
            @test death_simplex.(res_hom[i]) == death_simplex.(res_coh[i])
        end
    end
    @testset "Representative cycle" begin
        res_hom = ripserer(cycle, cohomology=false, reps=true, dim_max=3)
        @test vertices.(simplex.(representative(res_hom[2][1]))) == sort!(vcat(
            [SVector(i+1, i) for i in 1:17], [SVector(18, 1)]
        ))
    end
    @testset "Infinite intervals" begin
        @test_broken ripserer(cycle, cohomology=false, threshold=2)[2][1] == (1.0, Inf)
    end
    @testset "Cubical" begin
        data = zeros(5, 5, 5)
        data[2, 2:4, 2:4] .= 1
        data[3, :, :] .= [0 0 0 0 0; 0 1 1 1 0; 0 1 0 1 0; 0 1 1 1 0; 0 0 0 0 0]
        data[4, 2:4, 2:4] .= 1

        d0, d1, d2 = ripserer(Cubical(data); dim_max=2, cohomology=false)

        @test d0 == [(0, 1.0), (0, Inf)]
        @test isempty(d1)
        @test d2 == [(0, 1)]
    end
end

@testset "Only print to stderr and only when progress is enabled" begin
    @suppress begin
        @test (@capture_out ripserer(torus(16); dim_max=5)) == ""
        @test (@capture_out ripserer(torus(16); dim_max=5, progress=true)) == ""

        @test (@capture_err ripserer(torus(16); dim_max=5)) == ""
        @test (@capture_err ripserer(torus(16); dim_max=5, progress=true)) != ""
    end
end

@testset "Overflow checking" begin
    @test_throws OverflowError ripserer(Rips{Int16}(zeros(1000, 1000)))
end

# Julia 1.0 does not allow these to be defined inside @testset.
struct CustomRips <: Ripserer.AbstractRipsFiltration{Int, Float64} end
Ripserer.dist(::CustomRips) = Float64[i ≠ j for i in 1:10, j in 1:10]
Ripserer.dist(::CustomRips, i, j) = 1.0

@testset "Custom rips filtration" begin
    d0, d1 = ripserer(CustomRips())
    @test d0 == [fill((0.0, 1.0), 9); (0.0, Inf)]
    @test d1 == []
end

struct CustomFiltration <: Ripserer.AbstractFiltration{Int, Int} end

function Ripserer.unsafe_simplex(
    ::Type{Simplex{0, Int, Int}},
    ::CustomFiltration,
    (v,),
    sign,
)
    return Simplex{0}(sign * v, 0)
end
function Ripserer.unsafe_simplex(
    ::Type{Simplex{D, Int, Int}},
    ::CustomFiltration,
    vertices,
    sign,
) where D
    return Simplex{D}(sign * Ripserer._index(vertices), 1)
end
Ripserer.n_vertices(::CustomFiltration) = 10
Ripserer.simplex_type(::Type{CustomFiltration}, D) = Simplex{D, Int, Int}
Ripserer.edges(::CustomFiltration) = Simplex{1}.(10:-1:1, 1)
function Ripserer.postprocess_interval(::CustomFiltration, int::PersistenceInterval)
    return PersistenceInterval(birth(int) + 1, death(int) + 1)
end
function Ripserer.postprocess_interval(::CustomFiltration, ::RepresentativeInterval)
    return nothing
end

@testset "Custom filtration" begin
    d0, d1 = ripserer(CustomFiltration())
    @test d0 == [fill((1.0, 2.0), 4); fill((1.0, Inf), 6)]
    @test d1 == []

    d0, d1 = ripserer(CustomFiltration(), reps=true)
    @test d0 == []
    @test d1 == []
end
