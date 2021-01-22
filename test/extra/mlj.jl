using Ripserer
using Test
using MLJBase
using PersistenceDiagrams
using PersistenceDiagrams.MLJPersistenceDiagrams

include("../testdatasets.jl")

function test_mlj_model(model, X)
    mach = machine(model, X)
    fit!(mach; verbosity=0)
    Xt = transform(mach, X)

    @test nrows(Xt) == length(X)
    @test length(first(Xt)) ≥ 30

    model.vectorizer = PersistenceCurveVectorizer()
    fit!(mach; verbosity=0)
    Xt = transform(mach, X)

    @test nrows(Xt) == length(X)
    @test length(first(Xt)) == 20

    model.vectorizer.length = 100
    fit!(mach; verbosity=0)
    Xt = transform(mach, X)
    @test nrows(Xt) == length(X)
    @test length(first(Xt)) == 200

    model.dim_max = 2
    fit!(mach; verbosity=0)
    Xt = transform(mach, X)
    @test nrows(Xt) == length(X)
    @test length(first(Xt)) == 300
end

@testset "Rips" begin
    @testset "point-like input" begin
        X = [torus_points(100), torus_points(50), torus_points(10)]
        test_mlj_model(RipsPersistentHomology(), X)
    end
    @testset "distance matrix input, modulus" begin
        X = [icosahedron, cycle, projective_plane]

        model = RipsPersistentHomology()
        mach = machine(model, X)
        fit!(mach; verbosity=0)
        Xt1 = transform(mach, X)

        model.modulus = 3
        mach = machine(model, X)
        fit!(mach; verbosity=0)
        Xt2 = transform(mach, X)

        @test Xt1 ≠ Xt2
    end
end

@testset "Alpha" begin
    X = [torus_points(500), torus_points(400), torus_points(300)]
    test_mlj_model(AlphaPersistentHomology(), X)
end

@testset "Cubical" begin
    X = [rand(5, 10, 20), rand(10, 10, 10), rand(15, 5, 11)]
    test_mlj_model(CubicalPersistentHomology(), X)
end
