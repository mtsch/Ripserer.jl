using Ripserer

@test ∞ > 0.0
@test typemax(Int) < ∞
@test ∞ == Inf
@test Inf == ∞
@test !(∞ > Inf)
@test !(∞ < Inf)
@test !(∞ > NaN)
@test !(∞ < NaN)
@test !(∞ > ∞)
@test !(∞ < ∞)
@test !isless(∞, ∞)
@test ∞ == ∞
@test ∞ != NaN
@test ∞ != 1.0
@test ∞ != "infinity"
@test ∞ > "infinity plus one"
@test !isless(∞, "infinity plus one")
@test "infinity plus two" < ∞
@test isless("infinity plus two", ∞)
@test ∞ > ["infinity", +, 3]
@test ismissing(∞ > missing)
@test ismissing(∞ < missing)
@test !isless(∞, missing)
@test isless(missing, ∞)
@test ∞ ≈ Inf
@test Inf ≈ ∞
@test 1 ≉ Inf
@test Inf ≉ 1

@test sprint(print, ∞) == "∞"
