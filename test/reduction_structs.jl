using Ripserer
using Ripserer: ChainElement, coef,
    ReductionMatrix, insert_column!, has_column,
    Column, pop_pivot!, pivot,
    DisjointSetsWithBirth

using Compat

@testset "ReductionMatrix" begin
    for T in (Mod{3}, Rational{Int})
        CE = ChainElement{Simplex{3, Int, Int}, T}
        SE = ChainElement{Simplex{2, Int, Int}, T}
        rm = @inferred ReductionMatrix{Simplex{3, Int, Int}, SE}()
        @test length(rm) == 0
        @test eltype(rm) ≡ SE

        insert_column!(rm, -Simplex{3}(3, 1))
        push!(rm, SE(Simplex{2}(1, 1), 1))
        push!(rm, SE(Simplex{2}(2, 1), 1))
        push!(rm, SE(Simplex{2}(3, 1), 1))
        push!(rm, SE(Simplex{2}(4, 1), 1))

        insert_column!(rm, CE(Simplex{3}(10, 1), 1))
        push!(rm, SE(Simplex{2}(5, 2), 2))
        push!(rm, SE(Simplex{2}(6, 2), 2))
        push!(rm, SE(Simplex{2}(7, 2), 2))

        insert_column!(rm, Simplex{3}(1, 1))

        insert_column!(rm, Simplex{3}(15, 1))
        push!(rm, SE(Simplex{2}(8, 3), -1))

        @test has_column(rm, Simplex{3}(3, 1))
        @test collect(rm[Simplex{3}(3, 1)]) == SE.(Simplex{2}.([1, 2, 3, 4], 1), 1)

        @test has_column(rm, Simplex{3}(10, 1))
        @test length(rm[Simplex{3}(10, 1)]) == 3
        @test all(ce -> coef(ce) == 2, rm[Simplex{3}(10, 1)])

        @test has_column(rm, -Simplex{3}(1, 1))
        @test length(rm[Simplex{3}(1, 1)]) == 0
        @test eltype(collect(rm[Simplex{3}(1, 1)])) ≡ SE

        @test has_column(rm, CE(Simplex{3}(15, 1), 0))
        @test length(rm[CE(Simplex{3}(15, 1))]) == 1
        @test first(rm[CE(Simplex{3}(15, 1))]) == SE(Simplex{2}(8, 3), -1)

        @test !has_column(rm, Simplex{3}(2, 1))
        @test !has_column(rm, CE(Simplex{3}(100, 1)))
    end
end

@testset "Column" begin
    @testset "single element" begin
        CE = ChainElement{Simplex{1, Float64, Int}, Mod{2}}
        col = Column{CE}()
        push!(col, Simplex{1}(1, 2.0))
        push!(col, Simplex{1}(1, 2.0))
        push!(col, Simplex{1}(1, 2.0))
        push!(col, Simplex{1}(1, 2.0))
        push!(col, Simplex{1}(1, 2.0))

        @test pop_pivot!(col) == CE(Simplex{1}(1, 2.0), 1)
        @test isnothing(pivot(col))
        @test isempty(col)

        CE = ChainElement{Simplex{2, Float64, Int}, Mod{3}}
        col = Column{CE}()
        push!(col, Simplex{2}(1, 2.0))
        push!(col, Simplex{2}(1, 2.0))
        push!(col, Simplex{2}(1, 2.0))

        @test isnothing(pivot(col))
        @test isnothing(pop_pivot!(col))
        @test isempty(col)

        CE = ChainElement{Simplex{2, Int, Int}, Rational{Int}}
        col = Column{CE}()
        push!(col, Simplex{2}(1, 2))
        push!(col, Simplex{2}(2, 3))
        push!(col, Simplex{2}(-1, 2))
        push!(col, Simplex{2}(-2, 3))
        push!(col, -Simplex{2}(1, 2))

        @test coef(pivot(col)) == -1
        @test pop_pivot!(col) == CE(Simplex{2}(1, 2), -1)
        @test isnothing(pivot(col))
        @test isnothing(pop_pivot!(col))
        @test isempty(col)
    end
    @testset "multiple" begin
        CE = ChainElement{Simplex{3, Float64, Int}, Mod{5}}
        col = Column{CE}()
        push!(col, CE(Simplex{3}(2, 1.0), 3))
        push!(col, CE(Simplex{3}(3, 2.0), 4))
        push!(col, CE(Simplex{3}(2, 1.0), 2))
        push!(col, CE(Simplex{3}(1, 3.0), 2))
        push!(col, CE(Simplex{3}(3, 2.0), 1))
        push!(col, CE(Simplex{3}(4, 4.0), 4))
        push!(col, CE(Simplex{3}(4, 4.0), 4))
        push!(col, CE(Simplex{3}(4, 4.0), 4))
        push!(col, CE(Simplex{3}(4, 5.0), 4))
        push!(col, CE(Simplex{3}(4, 5.0), 1))

        @test pop_pivot!(col) == CE(Simplex{3}(1, 3.0), 2)
        @test pivot(col) == CE(Simplex{3}(4, 4.0), 2)
        @test pop_pivot!(col) == CE(Simplex{3}(4, 4.0), 2)
        @test isnothing(pop_pivot!(col))
        @test isnothing(pop_pivot!(col))
    end
end
