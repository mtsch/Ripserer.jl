using Ripserer

using DataStructures
using Ripserer: DisjointSetsWithBirth, find_leaves!

s = DisjointSetsWithBirth(collect(1:10))

@testset "basic tests" begin
    @test length(s) == 10
    for i = 1:10
        @test find_root!(s, i) == i
        @test find_leaves!(s, i) == [i]
        @test birth(s, i) == i
    end
    @test_throws BoundsError find_root!(s, 11)
end

@testset "union!, find_root!, find_leaves!" begin
    union!(s, 2, 3)
    @test find_root!(s, 3) == 2
    @test find_leaves!(s, 2) == [2, 3]

    @test union!(s, 8, 7) == 8
    @test union!(s, 5, 6) == 5
    @test find_leaves!(s, 5) == [5, 6]
    @test union!(s, 8, 5) == 8
    @test find_root!(s, 6) == 8
    @test find_leaves!(s, 8) == [5, 6, 7, 8]
    union!(s, 2, 6)
    @test find_root!(s, 2) == 8
    root1 = find_root!(s, 6)
    root2 = find_root!(s, 2)
    @test root_union!(s, Int(root1), Int(root2)) == 8
    @test union!(s, 5, 6) == 8
end
