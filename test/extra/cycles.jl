using LightGraphs
using Ripserer
using Test

using Ripserer: distance_matrix, OneSkeleton

@testset "OneSkeleton respects the LightGraphs interface" begin
    @testset "no threshold or removed simplices" begin
        flt = Rips(
            [
                0 1 2 3 2 1
                1 0 1 2 3 2
                2 1 0 1 2 3
                3 2 1 0 1 2
                2 3 2 1 0 1
                1 2 3 2 1 0
            ]
        )

        g = OneSkeleton(flt)
        @test eltype(g) ≡ Int
        @test nv(g) == 6
        @test ne(g) == binomial(6, 2)
        @test has_vertex(g, 1)
        @test has_vertex(g, 6)
        @test !has_vertex(g, 7)
        @test !has_vertex(g, 0)
        @test edges(g) == [Edge(j, i) for i in 1:5 for j in (i + 1):6]
        @test vertices(g) == 1:6
        @test is_directed(g) == is_directed(typeof(g)) == false
        @test weights(g) == distance_matrix(flt)
        @test has_edge(g, 1, 2)
        @test has_edge(g, 3, 6)
        @test !has_edge(g, 3, 7)
        @test neighbors(g, 1) == [6, 5, 4, 3, 2]
    end

    @testset "threshold and removed vertices" begin
        dists = [
            0 1 2 3 2 1
            1 0 1 2 3 2
            2 1 0 1 2 3
            3 2 1 0 1 2
            2 3 2 1 0 1
            1 2 3 2 1 0
        ]
        flt = Rips(dists)

        g = OneSkeleton(flt, 2, [Simplex{1}((2, 1), 1), Simplex{1}((3, 1), 2)])
        @test eltype(g) ≡ Int
        @test nv(g) == 6
        @test ne(g) == binomial(6, 2) - 2 - 3
        @test has_vertex(g, 1)
        @test has_vertex(g, 6)
        @test !has_vertex(g, 7)
        @test !has_vertex(g, 0)
        @test edges(g) == [
            Edge(j, i) for i in 1:5 for
            j in (i + 1):6 if dists[i, j] ≠ 3 && (j, i) ∉ [(2, 1), (3, 1)]
        ]
        @test vertices(g) == 1:6
        @test is_directed(g) == is_directed(typeof(g)) == false
        @test weights(g) == distance_matrix(flt)
        @test !has_edge(g, 1, 2)
        @test !has_edge(g, 3, 6)
        @test has_edge(g, 3, 5)
        @test !has_edge(g, 3, 7)
        @test neighbors(g, 1) == [6, 5]
    end

    @testset "simplex as threshold" begin
        dists = [
            0 1 1 1
            1 0 1 1
            1 1 0 1
            1 1 1 0
        ]
        flt = Rips(dists)
        g = OneSkeleton(flt, Simplex{1}((3, 2), 1))

        @test eltype(g) ≡ Int
        @test nv(g) == 4
        @test ne(g) == 4
        @test has_vertex(g, 1)
        @test has_vertex(g, 4)
        @test !has_vertex(g, 7)
        @test !has_vertex(g, 0)
        @test edges(g) == [Edge(j, i) for i in 1:3 for j in (i + 1):4 if (j, i) ≥ (3, 2)]
        @test vertices(g) == 1:4
        @test is_directed(g) == is_directed(typeof(g)) == false
        @test weights(g) == distance_matrix(flt)
        @test !has_edge(g, 1, 2)
        @test !has_edge(g, 3, 6)
        @test has_edge(g, 3, 2)
        @test has_edge(g, 4, 2)
        @test neighbors(g, 1) == [4]
        @test neighbors(g, 4) == [3, 2, 1]
    end

    @testset "Cubical" begin
        image = [
            0 0 0 0 0 0 0
            1 1 1 1 1 1 1
            1 1 2 3 2 1 2
            1 1 2 5 2 1 1
            1 1 1 1 1 1 1
        ]
        flt = Cubical(image)
        g = OneSkeleton(flt, 2)

        @test eltype(g) ≡ Int
        @test nv(g) == *(size(image)...)
        @test ne(g) == 4 * 6 * 2 + 4 + 6 - 7
        @test has_vertex(g, 1)
        @test has_vertex(g, 4)
        @test !has_vertex(g, 36)
        @test !has_vertex(g, 0)
        @test vertices(g) == 1:35
        @test is_directed(g) == is_directed(typeof(g)) == false
        @test weights(g) == distance_matrix(flt)
        @test has_edge(g, 1, 2)
        @test !has_edge(g, 3, 6)
        @test neighbors(g, 1) == [6, 2]
        @test neighbors(g, 14) == [15, 13, 9]
    end
end

@testset "Cycle reconstruction" begin
    @testset "small cycle" begin
        flt = Rips(
            [
                0 1 2 3 2 1
                1 0 1 2 3 2
                2 1 0 1 2 3
                3 2 1 0 1 2
                2 3 2 1 0 1
                1 2 3 2 1 0
            ]
        )

        interval = ripserer(flt; reps=true)[2][1]
        cyc1 = reconstruct_cycle(flt, interval)
        @test eltype(cyc1) == Simplex{1,Int,Int}
        @test all(birth.(cyc1) .== 1)
        @test allunique(cyc1)
        @test length(cyc1) == 6

        cyc2 = reconstruct_cycle(flt, interval, 2)
        @test isempty(cyc2)
    end

    @testset "circle points" begin
        pts = [(sin(t), cos(t)) for t in range(0, 2π; length=22)[1:(end - 1)]]
        flt = Rips(pts)
        interval = ripserer(flt; reps=true)[2][1]

        prev_length = nv(flt)
        for r in range(birth(interval), death(interval); length=10)[1:(end - 1)]
            cyc = reconstruct_cycle(flt, interval)
            @test maximum(birth, cyc) ≤ r

            # Sum of lengths approximates length around circle
            @test π ≤ sum(birth, cyc) ≤ 2π

            # Length decreases for larger r
            @test length(cyc) ≤ prev_length
            prev_length = length(cyc)

            # Each vertex appears exactly twice
            vxs = sort!(collect(Iterators.flatten(vertices.(cyc))))
            @test vxs[1:2:end] == vxs[2:2:end]
            @test allunique(vxs[1:2:end])
        end
        @test isempty(reconstruct_cycle(flt, interval, death(interval)))
        @test isempty(reconstruct_cycle(flt, interval, birth(interval) - eps(Float64)))
    end

    @testset "circle points over a grid with a hole" begin
        # The idea here is that at birth time, the cycle is a circle and at time 1, the
        # cycle is a square surrounding the hole.
        pts = unique!(
            vcat(
                [(3sin(t), 3cos(t)) for t in range(0, 2π; length=22)[1:(end - 1)]],
                vec([
                    (i - 5, j - 5) for
                    (i, j) in Iterators.product(0.0:10.0, 0.0:10.0) if !(i ∈ 4:6 && j ∈ 4:6)
                ]),
            ),
        )

        flt = Rips(pts)
        interval = sort!(ripserer(flt; reps=true)[2]; by=persistence)[end]

        cyc1 = reconstruct_cycle(flt, interval)
        cyc2 = reconstruct_cycle(flt, interval, 1)
        cyc3 = reconstruct_cycle(flt, interval, 2)
        cyc4 = reconstruct_cycle(flt, interval, 3)

        @test sum(birth, cyc1) ≈ 6π atol = 0.1
        @test sum(birth, cyc2) == 16
        @test sum(birth, cyc3) ≈ 8 + 4 * √2
        @test sum(birth, cyc4) ≈ 8 * √2
    end

    @testset "Cubical" begin
        image = [
            0 0 0 0 0 0 0
            1 1 1 1 1 1 1
            1 1 2 3 2 1 2
            1 1 2 5 2 1 1
            1 1 1 1 1 1 1
        ]
        cf = Cubical(image)

        interval = ripserer(cf; reps=true)[2][1]
        c(t1, t2) = simplex(cf, Val(1), (CartesianIndex(t1), CartesianIndex(t2)))
        @test sort(reconstruct_cycle(cf, interval, 1)) == sort([
            c((5, 2), (4, 2)),
            c((4, 2), (3, 2)),
            c((3, 2), (2, 2)),
            c((2, 2), (2, 3)),
            c((2, 3), (2, 4)),
            c((2, 4), (2, 5)),
            c((2, 5), (2, 6)),
            c((2, 6), (3, 6)),
            c((3, 6), (4, 6)),
            c((4, 6), (5, 6)),
            c((5, 6), (5, 5)),
            c((5, 5), (5, 4)),
            c((5, 4), (5, 3)),
            c((5, 3), (5, 2)),
        ])
        @test sort(reconstruct_cycle(cf, interval, 2)) == sort([
            c((5, 3), (4, 3)),
            c((4, 3), (3, 3)),
            c((3, 3), (2, 3)),
            c((2, 3), (2, 4)),
            c((2, 4), (2, 5)),
            c((2, 5), (3, 5)),
            c((3, 5), (4, 5)),
            c((4, 5), (5, 5)),
            c((5, 5), (5, 4)),
            c((5, 4), (5, 3)),
        ])
        @test sort(reconstruct_cycle(cf, interval, 3)) == sort([
            c((5, 3), (4, 3)),
            c((4, 3), (3, 3)),
            c((3, 3), (3, 4)),
            c((3, 4), (3, 5)),
            c((3, 5), (4, 5)),
            c((4, 5), (5, 5)),
            c((5, 5), (5, 4)),
            c((5, 4), (5, 3)),
        ])
        @test reconstruct_cycle(cf, interval, 4) == reconstruct_cycle(cf, interval, 3)
        @test isempty(reconstruct_cycle(cf, interval, 0))
        @test isempty(reconstruct_cycle(cf, interval, 5))
    end

    @testset "Overriding distance matrix" begin
        flt = Rips([
            0 2 4 2 2
            2 0 2 4 2
            4 2 0 2 6
            2 4 2 0 6
            2 2 6 6 0
        ])

        interval = ripserer(flt; reps=true)[2][1]

        dists = [
            0 2 4 2 0
            2 0 2 4 0
            4 2 0 2 6
            2 4 2 0 6
            0 0 6 6 0
        ]
        cyc1 = reconstruct_cycle(flt, interval, 2)
        cyc2 = reconstruct_cycle(flt, interval, 2; distances=dists)
        @test length(cyc1) == 4
        @test length(cyc2) == 5
    end

    @testset "Errors" begin
        flt = Rips(
            [
                0 1 2 3 2 1
                1 0 1 2 3 2
                2 1 0 1 2 3
                3 2 1 0 1 2
                2 3 2 1 0 1
                1 2 3 2 1 0
            ]
        )
        _, d1, d2 = ripserer(flt; reps=2, dim_max=2)

        @test_throws ArgumentError reconstruct_cycle(flt, d1[1])
        @test_throws ArgumentError reconstruct_cycle(flt, d2[1])
    end
end
