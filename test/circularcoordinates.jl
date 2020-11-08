using Distances
using LinearAlgebra
using Ripserer
using StaticArrays
using Test

using Ripserer:
    _landmarks_and_radius,
    _to_integer_coefficients,
    _zero_coboundary_matrix,
    _to_vector,
    _harmonic_smoothing,
    is_cocycle,
    ChainElement

include("filtrations/test-datasets.jl")

@testset "Helpers" begin
    @testset "Partition functions" begin
        for fun in (Partition.linear, Partition.quadratic, Partition.exponential)
            for r in range(1, 2; length=5), d in range(0, 2; length=5)
                if r > d
                    @test fun(r, d) > 0
                else
                    @test fun(r, d) == 0
                end
            end
        end
    end

    @testset "Landmark selection" begin
        points = [(sin(t), cos(t)) for t in range(0, 2π; length=41)[1:40]]
        point_matrix = reshape(reinterpret(Float64, points), (2, 40))

        for landmark_selector in (10, 1:4:40, points[1:4:40])
            landmarks, radius = _landmarks_and_radius(
                points, landmark_selector, Euclidean()
            )
            @test length(landmarks) == 10
            @test radius ≈ 0.3128689 atol = 1e-5
            @test eltype(landmarks) <: SVector{2,Float64}
        end
    end

    @testset "Cocycles" begin
        rips = Rips([
            0.0 1.1 2.0 2.0 1.2
            1.1 0.0 1.0 2.0 2.0
            2.0 1.0 0.0 1.0 2.0
            2.0 2.0 1.0 0.0 1.0
            1.2 2.0 2.0 1.0 0.0
        ])

        cocycle = ripserer(rips; modulus=17, reps=true)[2][1].representative
        integral = _to_integer_coefficients(cocycle)
        @test simplex.(integral) == simplex.(cocycle)
        @test coefficient.(integral) == [-1, -1, -1]
        @test _to_vector(rips, integral, Inf) == [0, -1, 0, -1, 0, 0, -1, 0, 0, 0]
        @test _to_vector(rips, integral, 2) == [0, 0, 0, 0, 0, 0, -1, 0, 0, 0]

        @test _zero_coboundary_matrix(rips, Inf) == [
            1 -1 0 0 0
            1 0 -1 0 0
            0 1 -1 0 0
            1 0 0 -1 0
            0 1 0 -1 0
            0 0 1 -1 0
            1 0 0 0 -1
            0 1 0 0 -1
            0 0 1 0 -1
            0 0 0 1 -1
        ]
        @test _zero_coboundary_matrix(rips, 2) == [
            1 -1 0 0 0
            0 0 0 0 0
            0 1 -1 0 0
            0 0 0 0 0
            0 0 0 0 0
            0 0 1 -1 0
            1 0 0 0 -1
            0 0 0 0 0
            0 0 0 0 0
            0 0 0 1 -1
        ]
        coords, harmonic = _harmonic_smoothing(rips, cocycle; time=2)
        @test all(diff(coords) .≈ 0.2)

        # The harmonic interval should have a minimal norm.
        @test norm(harmonic) < norm(_to_vector(rips, integral, 2))
    end
end

@testset "CircularCoordinates" begin
    @testset "landmark selection" begin
        points = torus_points(1024)

        cc1 = CircularCoordinates(points, 100)
        @test length(cc1.landmarks) == 100
        @test allunique(cc1.landmarks)

        cc2 = CircularCoordinates(points, 1:100)
        @test cc2.landmarks == SVector.(points[1:100])

        cc3 = CircularCoordinates(points, points[1:100])
        @test cc3.landmarks == SVector.(points[1:100])

        @test cc2.coordinate_data[1].radius == cc3.coordinate_data[1].radius
    end
    @testset "selecting out_dim" begin
        points = torus_points(1024)

        cc = CircularCoordinates(points, 1:10:1024)
        @test cc.out_dim == 1

        cc = CircularCoordinates(points, 1:10:1024; out_dim=3, warn=false)
        @test cc.out_dim == 1
        @test size(cc(points[1:10])) == (10, 1)
        @test size(cc(points[11:20], 1)) == (10,)

        cc = CircularCoordinates(points, 1:5:1024; out_dim=3, warn=false)
        @test cc.out_dim == 2
        out_both = cc(points[1:10])
        out_first = cc(points[1:10], 1)
        out_second = cc(points[1:10], 2)
        @test out_both[:, 1] == out_first
        @test out_both[:, 2] == out_second

        @test_throws ErrorException CircularCoordinates(points, 1:10)
    end
    @testset "transformation" begin
        ts1 = range(0, 1; length=100)[1:(end - 1)]
        circle = [(sinpi(2t), cospi(2t)) for t in ts1]

        cc = CircularCoordinates(circle, 1:10:100)
        transformed1 = cc(circle, 1)
        transformed1 .-= transformed1[1]
        transformed1 .= Ripserer._mod_z.(transformed1)

        @test transformed1 ≈ ts1 atol = 0.05

        @test cc([(100, 100)])[1] === missing

        ts2 = [0; rand(20)]
        new_points = [(sinpi(2t), cospi(2t)) for t in ts2]
        transformed2 = cc(new_points, 1)
        transformed2 .-= transformed2[1]
        transformed2 .= Ripserer._mod_z.(transformed2)

        @test transformed2 ≈ ts2 atol = 0.05
    end
end
