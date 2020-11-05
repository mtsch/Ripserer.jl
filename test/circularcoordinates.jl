using Distances
using LinearAlgebra
using Ripserer
using StaticArrays
using Test

using Ripserer:
    _landmarks_and_radius, _to_integer_coefficients, _to_vector, _harmonic_smoothing,
    is_cocycle, ChainElement

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
        points = [(sin(t), cos(t)) for t in range(0, 2π, length=11)[1:10]]
        rips = Rips(points)
        interval = ripserer(rips; modulus=17, reps=true)[end][end]
        cocycle = interval.representative

        int_cocycle = _to_integer_coefficients(cocycle)
        @test is_cocycle(Rips(points), int_cocycle, death(interval))
        @test eltype(int_cocycle) <: ChainElement{Simplex{1,Float64,Int},Int}
        vec = _to_vector(rips, int_cocycle, death(interval))
        @test length(vec) == 90
        @test count(vec .≠ 0) == count(birth.(cocycle) .< death(interval))

        @test_throws ErrorException _harmonic_smoothing(rips, cocycle; time=Inf)
        mini, harm = _harmonic_smoothing(rips, cocycle; time=death(interval))
        @test length(mini) == 10
        @test length(harm) == 90
    end
end
