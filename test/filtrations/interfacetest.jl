using Test
using Ripserer
using Suppressor

const R = Ripserer

function test_filtration(F, args...; test_verbose=true, flt_kwargs=(), kwargs...)
    if !isempty(kwargs)
        desc = "$F with $(summary.(args)) and $(kwargs...)"
    else
        desc = "$F with $(summary.(args))"
    end
    all_kw = (; flt_kwargs..., kwargs...)
    @testset "$desc" begin
        flt = F(args...; flt_kwargs...)

        @testset "Interface" begin
            @test length(R.vertices(flt)) == R.nv(flt)
            @test length(R.births(flt)) == R.nv(flt)
            edges = R.edges(flt)
            @test first(edges) isa R.simplex_type(flt, 1)
            @test R.simplex(flt, Val(1), R.vertices(first(edges))) == first(edges)

            triangles = R.columns_to_reduce(flt, edges)
            @test first(triangles) isa R.simplex_type(flt, 2)
            @test R.simplex(flt, Val(2), R.vertices(first(triangles))) == first(triangles)
            @test R.distance_matrix(flt) isa AbstractMatrix

            @test R.threshold(flt) isa eltype(R.births(flt))
            @test R.emergent_pairs(flt) isa Bool
        end

        @testset "Basic ripserer" begin
            @test ripserer(F, args...; all_kw...) == ripserer(flt; kwargs...)
        end
        if test_verbose
            @testset "Verbose ripserer" begin
                @suppress begin
                    @test (@capture_out ripserer(F, args...; all_kw...)) == ""
                    @test (@capture_out ripserer(F, args...; verbose=true, all_kw...)) == ""

                    @test (@capture_err ripserer(F, args...; all_kw...)) == ""
                    @test (@capture_err ripserer(F, args...; verbose=true, all_kw...)) != ""
                end
            end
        end
        @testset "Algorithms and representatives" begin
            coh_imp = ripserer(flt; implicit=true, reps=true, kwargs...)
            coh_exp = ripserer(flt; implicit=false, kwargs...)
            hom_imp = ripserer(flt; alg=:homology, implicit=true, reps=true, kwargs...)
            hom_exp = ripserer(flt; alg=:homology, implicit=false, reps=true, kwargs...)
            hom_inv = ripserer(flt; alg=:involuted, reps=true, kwargs...)

            n_diagrams = length(coh_imp)
            for i in 2:n_diagrams
                @test coh_imp[i] == coh_exp[i]
                @test hom_imp[i] == hom_exp[i] == hom_inv[i]
                @test hom_inv[i] == coh_imp[i]
                @test representative.(hom_imp[i]) == representative.(hom_exp[i])
                @test representative.(hom_inv[i]) == representative.(hom_exp[i])
            end
        end
    end
end
