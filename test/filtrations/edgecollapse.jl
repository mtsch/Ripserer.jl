using Ripserer
using Test

include("test-datasets.jl")

@testset "Result equal to Rips" begin
    for m in (2, 3)
        for t in (2, 3)
            for data in (projective_plane, [tuple(3 .* rand(5)...) for _ in 1:20])
                rips = Rips(data; threshold=t)
                col1 = EdgeCollapsedRips(rips)
                col2 = EdgeCollapsedRips(data; threshold=t)
                col3 = EdgeCollapsedRips{Int32}(EdgeCollapsedRips(data); threshold=t)
                res_rips = ripserer(rips; modulus=m, dim_max=2)
                res_col1 = ripserer(col1; modulus=m, dim_max=2)
                res_col2 = ripserer(col2; modulus=m, dim_max=2)
                res_col3 = ripserer(col3; modulus=m, dim_max=2)

                ne_rips = length(Ripserer.edges(rips))
                ne_col1 = length(Ripserer.edges(col1))
                ne_col2 = length(Ripserer.edges(col2))
                ne_col3 = length(Ripserer.edges(col3))

                @test Ripserer.simplex_type(col3, 2) == Simplex{2,eltype(rips.adj),Int32}

                @test ne_rips > ne_col1
                @test ne_col1 == ne_col2 == ne_col3
                @test res_rips == res_col1 == res_col2 == res_col3
            end
        end
    end
end
