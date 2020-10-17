using Compat
using Distances
using Ripserer
using SparseArrays
using Random
using StaticArrays
using Suppressor
using Test

using Ripserer:
    distances,
    births,
    adjacency_matrix,
    edges,
    nv,
    unsafe_simplex,
    ChainElement,
    PackedElement

include("test-datasets.jl")

@testset "distances" begin
    for points in (
        [(0, 0), (0, 1), (1, 1), (1, 0)],
        [SVector(0, 0), SVector(0, 1), SVector(1, 1), SVector(1, 0)],
    )
        @test distances(points) ≈ [0 1 √2 1; 1 0 1 √2; √2 1 0 1; 1 √2 1 0]
        @test distances(points, Cityblock()) == [0 1 2 1; 1 0 1 2; 2 1 0 1; 1 2 1 0]
    end
end

for sparse in (true, false)
    @testset "With no threshold, sparse=$sparse" begin
        filtration = Rips(Float64[0 1 2 9; 1 0 3 9; 2 3 0 4; 9 9 4 0]; sparse=sparse)

        @test sprint(show, filtration) == "Rips{Int64, Float64}(nv=4, sparse=$sparse)"

        @test nv(filtration) == 4
        if issparse(adjacency_matrix(filtration))
            @test adjacency_matrix(filtration) == [0 1 2 0; 1 0 3 0; 2 3 0 4; 0 0 4 0]
        else
            @test adjacency_matrix(filtration) == [0 1 2 9; 1 0 3 9; 2 3 0 4; 9 9 4 0]
        end

        @test threshold(filtration) == 4.0
        @test births(filtration) == zeros(4)

        @test adjacency_matrix(filtration) == filtration.adj

        @test unsafe_simplex(filtration, Val(0), (1,)) === Simplex{0}(1, 0.0)
        @test unsafe_simplex(filtration, Val(1), (2, 1), -1) === -Simplex{1}(1, 1.0)
        @test simplex(filtration, Val(2), (1, 3, 2)) === Simplex{2}(1, 3.0)
        @test_throws ArgumentError simplex(filtration, Val(2), (1, 1, 2))
    end
    @testset "With threshold and index type, sparse=$sparse" begin
        filtration = Rips{Int128}([1 2 2; 2 2 3; 2 3 2]; threshold=2, sparse=sparse)

        @test sprint(show, filtration) == "Rips{Int128, Int64}(nv=3, sparse=$sparse)"

        @test nv(filtration) == 3
        @test threshold(filtration) == 2
        if sparse
            @test adjacency_matrix(filtration) == [1 2 2; 2 2 0; 2 0 2]
        else
            @test adjacency_matrix(filtration) == [1 2 2; 2 2 3; 2 3 2]
        end
        @test threshold(filtration) == 2

        @test births(filtration) == [1, 2, 2]

        @test adjacency_matrix(filtration) == filtration.adj

        @test unsafe_simplex(filtration, Val(0), (1,)) === Simplex{0,Int,Int128}(1, 1)
        @test unsafe_simplex(filtration, Val(2), (3, 2, 1)) === nothing
        @test simplex(filtration, Val(1), (1, 2), -1) === -Simplex{1,Int,Int128}(1, 2)
        @test_throws ArgumentError simplex(filtration, Val(2), (0, 1, 2))
    end
end

@testset "Rips points constructor, sparse=false" begin
    filtration = Rips([(sin(x), cos(x)) for x in range(0, 2π; length=101)[1:(end - 1)]])
    adj = adjacency_matrix(filtration)
    @test all(x -> x > 0, adj[i, j] for i in 1:100 for j in (i + 1):100)
    @test eltype(edges(filtration)) === Simplex{1,Float64,Int}

    filtration = Rips{Int32}([
        (sin(x), cos(x)) for x in range(0.0f0, 2.0f0π; length=101)[1:(end - 1)]
    ])
    adj = adjacency_matrix(filtration)
    @test !issparse(adj)
    @test all(x -> x > 0, adj[i, j] for i in 1:100 for j in (i + 1):100)
    @test eltype(edges(filtration)) === Simplex{1,Float32,Int32}
end

@testset "Rips points constructor, sparse=true" begin
    filtration = Rips(
        [(sin(x), cos(x)) for x in range(0, 2π; length=101)[1:(end - 1)]];
        threshold=0.1,
        sparse=true,
    )
    adj = adjacency_matrix(filtration)
    @test issparse(adj)
    @test maximum(adj) ≤ 0.1
    @test maximum(adj) ≤ threshold(filtration)
end

@testset "Construction does not alter input" begin
    dist = [
        0 9 1 2
        9 0 3 4
        1 3 0 4
        2 4 4 0
    ]
    orig_dist = copy(dist)
    Rips(dist; threshold=1)
    @test dist == orig_dist
    Rips(dist; threshold=1, sparse=true)
    @test dist == orig_dist

    dist = sparse(dist)
    orig_dist = copy(dist)
    Rips(dist; threshold=1)
    @test dist == orig_dist
end

@testset "Warns with duplicate points" begin
    pts = [(1, 0), (1, 1), (0, 1), (0, 0), (0, 0)]
    @test @capture_err(Rips(pts)) ≠ ""
    @test @capture_err(Rips(pts[1:4])) == ""
    @test nv(@suppress Rips(pts)) == 4
end

@testset "Errors" begin
    @testset "Non-square matrices throw an error" begin
        @test_throws DimensionMismatch Rips(zeros(3, 2))
        @test_throws DimensionMismatch Rips(zeros(3, 2); sparse=true)
    end
    @testset "Asymmetric matrices throw an error" begin
        @test_throws ArgumentError Rips([1 1 1; 1 1 1; 1 2 1])
        @test_throws ArgumentError Rips([1 1 1; 2 1 1; 1 2 1]; sparse=true)
    end
    @testset "Edge births must be larger than vertex births" begin
        @test_throws ArgumentError Rips([1 1 1; 1 1 1; 1 1 2])
        @test_throws ArgumentError Rips([1 1 1; 1 2 1; 1 1 1]; sparse=true)
    end
end

@testset "ripserer" begin
    @testset "Rips" begin
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
            _, d1_2, d2_2 = ripserer(projective_plane; dim_max=2, threshold=1)
            _, d1_3, d2_3 = ripserer(projective_plane; dim_max=2, modulus=3, threshold=1)
            _, d1_331, d2_331 = ripserer(
                projective_plane; dim_max=2, field_type=Mod{5}, threshold=1
            )
            _, d1_r, d2_r = ripserer(
                projective_plane; dim_max=2, field_type=Rational{Int}, threshold=1
            )
            @test d1_2 == [(1, Inf)]
            @test d2_2 == [(1, Inf)]
            @test d1_3 == d1_331 == d1_r == []
            @test d2_3 == d2_331 == d1_r == []
        end
        @testset "Points as input" begin
            for metric in (Euclidean(), Cityblock())
                pts = torus_points(9)
                @test ripserer(pts; metric=metric) == ripserer(Rips(pts; metric=metric))
            end
        end
        @testset "Cutoff" begin
            d0, d1 = ripserer(rand_dist_matrix(20); cutoff=0.5)
            @test all(persistence.(d0) .> 0.5)
            @test all(persistence.(d1) .> 0.5)
        end
    end

    @testset "sparse" begin
        @testset "Icosahedron" begin
            d0, d1, d2 = ripserer(sparse(icosahedron); dim_max=2)
            @test d0 == [fill((0.0, 1.0), 11); (0.0, Inf)]
            @test d1 == []
            @test d2 == [(1.0, 2.0)]
        end
        @testset "RP2 with various fields, threshold=1" begin
            _, d1_2, d2_2 = ripserer(sparse(projective_plane); dim_max=2, threshold=1)
            _, d1_3, d2_3 = ripserer(
                Rips(projective_plane; threshold=1, sparse=true); dim_max=2, modulus=3
            )
            _, d1_331, d2_331 = ripserer(
                sparse(projective_plane); dim_max=2, field_type=Mod{5}, threshold=1
            )
            _, d1_r, d2_r = ripserer(
                sparse(projective_plane); dim_max=2, field_type=Rational{Int}, threshold=1
            )
            @test d1_2 == [(1, Inf)]
            @test d2_2 == [(1, Inf)]
            @test d1_3 == d1_331 == d1_r == []
            @test d2_3 == d2_331 == d1_r == []
        end
        @testset "Equal to Rips" begin
            for thresh in (nothing, 1, 0.5, 0.126)
                data = torus_points(100)
                r_res = ripserer(data; threshold=thresh, dim_max=2)
                s_res_1 = ripserer(data; threshold=thresh, sparse=true, dim_max=2)

                # Add zeros to diagonal. Adding ones first actually changes the structure of
                # the matrix.
                data2 = sparse(Ripserer.distances(data))
                for i in axes(data2, 1)
                    data2[i, i] = 1
                    data2[i, i] = 0
                end
                s_res_2 = ripserer(data2; threshold=thresh, dim_max=2)

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
            @test eltype(d0) isa DataType
            @test eltype(d1) isa DataType
            @test eltype(d2) isa DataType
            @test eltype(d3) isa DataType
            @test eltype(d1[1].representative) <: PackedElement

            d0, d1, d2, d3 = ripserer(cycle; dim_max=3, reps=true, field_type=Rational{Int})
            @test eltype(d0) isa DataType
            @test eltype(d1) isa DataType
            @test eltype(d2) isa DataType
            @test eltype(d3) isa DataType
            @test eltype(d3[1].representative) <: ChainElement
        end
        @testset "Infinite interval very low threshold" begin
            _, d1 = ripserer(
                cycle; dim_max=1, reps=true, threshold=1, field_type=Rational{Int}
            )
            rep = representative(only(d1))
            @test isempty(rep)
            @test eltype(rep) == ChainElement{Simplex{1,Int,Int},Rational{Int}}
        end
        @testset "Infinite interval higher threshold" begin
            _, d1 = ripserer(
                cycle; dim_max=1, reps=true, threshold=3, field_type=Rational{Int}
            )
            rep = representative(only(d1))
            @test !isempty(rep)
            @test eltype(rep) == ChainElement{Simplex{1,Int,Int},Rational{Int}}
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

    @testset "Diagram metadata" begin
        filtration = Rips(cycle)
        d0, d1, d2, d3 = ripserer(
            filtration; dim_max=3, reps=true, field_type=Rational{Int}
        )
        @test d0.dim == 0
        @test d1.dim == 1
        @test d2.dim == 2
        @test d3.dim == 3
        thresh = Float64(threshold(Rips(cycle)))
        @test d0.threshold ≡ thresh
        @test d1.threshold ≡ thresh
        @test d2.threshold ≡ thresh
        @test d3.threshold ≡ thresh
        field_type = Rational{Int}
        @test d0.field_type ≡ field_type
        @test d1.field_type ≡ field_type
        @test d2.field_type ≡ field_type
        @test d3.field_type ≡ field_type
        @test d0.filtration == filtration
        @test d1.filtration == filtration
        @test d2.filtration == filtration
        @test d3.filtration == filtration
    end

    @testset "Zero-dimensional sublevel set persistence" begin
        @testset "with sparse matrix" begin
            data = [1, 0, 1, 2, 3, 4, 3, 2, 3, 2, 1, 2]

            # Create distance matrix from data, where neighboring points are connected by edges
            # and the edge weights are equal to the max of both vertex births.
            n = length(data)
            dists = spzeros(n, n)
            for i in 1:n
                dists[i, i] = data[i]
            end
            for i in 1:(n - 1)
                j = i + 1
                dists[i, j] = dists[j, i] = max(dists[i, i], dists[j, j])
            end
            # 0-dimensional persistence should find values of minima and maxima of our data.
            d0 = ripserer(dists; dim_max=0)[1]
            @test d0 == [(2, 3), (1, 4), (0, Inf)]
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
            @test d0 == [(2, 3), (1, 4), (0, Inf)]
        end
    end

    @testset "Homology and explicit cohomology" begin
        for m in (2, 17)
            @testset "Same for projective plane with modulus=$m" begin
                r = Rips(projective_plane)
                _, hom_imp1, hom_imp2 = ripserer(
                    r; alg=:homology, implicit=true, dim_max=2, modulus=m
                )
                _, hom_exp1, hom_exp2 = ripserer(
                    r; alg=:homology, implicit=false, dim_max=2, modulus=m
                )
                _, hom_ass1, hom_ass2 = ripserer(r; alg=:involuted, dim_max=2, modulus=m)
                _, coh_imp1, coh_imp2 = ripserer(
                    r; implicit=true, reps=true, dim_max=2, modulus=m
                )
                _, coh_exp1, coh_exp2 = ripserer(
                    r; implicit=true, reps=true, dim_max=2, modulus=m
                )
                @test hom_imp1 == hom_exp1 == hom_ass1 == coh_imp1 == coh_exp1
                @test hom_imp2 == hom_exp2 == hom_ass2 == coh_imp2 == coh_exp2
                @test representative.(hom_imp1) == representative.(hom_exp1)
                @test representative.(hom_imp1) == representative.(hom_ass1)
                @test representative.(hom_imp2) == representative.(hom_exp2)
                @test representative.(hom_imp2) == representative.(hom_ass2)
                @test representative.(coh_imp1) == representative.(coh_exp1)
                @test representative.(coh_imp2) == representative.(coh_exp2)
            end
        end
        @testset "Representative cycle" begin
            res_hom = ripserer(cycle; alg=:homology, reps=true, dim_max=3)
            @test vertices.(simplex.(representative(res_hom[2][1]))) ==
                  sort!(vcat([SVector(i + 1, i) for i in 1:17], [SVector(18, 1)]))
        end
        @testset "Infinite intervals" begin
            @test_broken ripserer(cycle; alg=:homology, threshold=2, implicit=true)[2][1] ==
                         (1.0, Inf)
            @test_broken ripserer(cycle; alg=:homology, threshold=2, implicit=false)[2][1] ==
                         (1.0, Inf)
            @test ripserer(cycle; alg=:involuted, threshold=2)[2][1] == (1.0, Inf)
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
        @test_throws OverflowError ripserer(Rips{Int16}(ones(1000, 1000)))
    end

    @testset "Unsupported algirithms" begin
        @test_throws ArgumentError ripserer(Rips(ones(5, 5)); alg=:something)
    end
end
