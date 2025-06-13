using Ripserer
using Test

using DataStructures
using Ripserer: DisjointSetsWithBirth, find_leaves!

for arr in (1:10, CartesianIndices([1 2 3 4 5; 1 2 3 4 5]))
    @testset "with $(typeof(collect(arr)))" begin
        s = DisjointSetsWithBirth(arr, reshape(1:10, size(arr)))

        @testset "basic tests" begin
            @test length(s) == 10
            for i in 1:10
                @test find_root!(s, arr[i]) == arr[i]
                @test find_leaves!(s, arr[i]) == [arr[i]]
                @test birth(s, arr[i]) == (i, arr[i])
            end
            @test_throws BoundsError find_root!(s, arr[11])
        end

        @testset "union!, find_root!, find_leaves!" begin
            union!(s, arr[2], arr[3])
            @test find_root!(s, arr[3]) == arr[2]
            @test find_leaves!(s, arr[2]) == [arr[2], arr[3]]

            @test union!(s, arr[8], arr[7]) == arr[8]
            @test union!(s, arr[5], arr[6]) == arr[5]
            @test find_leaves!(s, arr[5]) == [arr[5], arr[6]]
            @test union!(s, arr[8], arr[5]) == arr[8]
            @test find_root!(s, arr[6]) == arr[8]
            @test find_leaves!(s, arr[8]) == [arr[5], arr[6], arr[7], arr[8]]
            union!(s, arr[2], arr[6])
            @test find_root!(s, arr[2]) == arr[8]
            root1 = find_root!(s, arr[6])
            root2 = find_root!(s, arr[2])
            @test root_union!(s, root1, root2) == arr[8]
            @test union!(s, arr[5], arr[6]) == arr[8]
        end
    end
end

@testset "merge tree" begin
    for data in
    end
end
