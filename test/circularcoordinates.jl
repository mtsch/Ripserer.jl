using Distances
using LinearAlgebra
using Ripserer
using StaticArrays
using Test

using Ripserer: _landmarks_and_radius, _to_integer_coefficients, _zero_coboundary_matrix,
    _to_vector, _harmonic_smoothing, is_cocycle, ChainElement

@testset "Helpers" begin
    @testset "Landmark selection" begin
        points = [(sin(t), cos(t)) for t in range(0, 2π, length=41)[1:40]]
        point_matrix = reshape(reinterpret(Float64, points), (2, 40))

        for landmark_selector in (10, 1:4:40, points[1:4:40])
            landmarks, radius = _landmarks_and_radius(
                points, landmark_selector, Euclidean()
            )
            @test length(landmarks) == 10
            @test radius ≈ 0.3128689 atol=1e-5
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
            1 -1  0  0  0
            1  0 -1  0  0
            0  1 -1  0  0
            1  0  0 -1  0
            0  1  0 -1  0
            0  0  1 -1  0
            1  0  0  0 -1
            0  1  0  0 -1
            0  0  1  0 -1
            0  0  0  1 -1
        ]
        @test _zero_coboundary_matrix(rips, 2) == [
            1 -1  0  0  0
            0  0  0  0  0
            0  1 -1  0  0
            0  0  0  0  0
            0  0  0  0  0
            0  0  1 -1  0
            1  0  0  0 -1
            0  0  0  0  0
            0  0  0  0  0
            0  0  0  1 -1
        ]
        coords, harmonic = _harmonic_smoothing(rips, cocycle; time=2)
        @test all(diff(coords) .≈ 0.2)

        # The harmonic interval should have a minimal norm.
        @test norm(harmonic) < norm(_to_vector(rips, integral, 2))
    end
end
