"""
    compute_death_simplices!(matrix::CoboundaryMatrix{true}, verbose, cutoff)

Fully reduce `matrix`, but only compute (homological) death simplices. Return all death
simplices up to the last that produces an interval with persistence greater than `cutoff`.

Used for involuted homology.
"""
function compute_death_simplices!(matrix::CoboundaryMatrix{true}, verbose, cutoff)
    columns, apparent = find_apparent_pairs(
        matrix.filtration, matrix.columns_to_reduce, verbose
    )
    bulk_insert!(matrix.reduced, apparent)
    deaths = simplex_type(matrix.filtration, dim(matrix) + 1)[]
    inf_births = simplex_type(matrix.filtration, dim(matrix))[]
    if isempty(columns)
        return deaths, inf_births
    else
        dim(matrix) > 1 && sort!(columns; rev=true)
        thresh = typemin(birth(first(columns)))
        for pair in apparent
            if birth(pair[2]) - birth(pair[1]) > cutoff
                thresh = max(thresh, birth(pair[2]))
            end
            push!(deaths, pair[2])
        end
        if verbose
            progbar = Progress(length(columns); desc="Precomputing columns...   ")
        end
        for column in columns
            pivot = reduce_column!(matrix, column)
            if !isnothing(pivot)
                if birth(pivot) - birth(column) > cutoff
                    thresh = max(thresh, birth(pivot))
                end
                push!(deaths, simplex(pivot))
            else
                push!(inf_births, column)
            end
            verbose && next!(progbar; showvalues=((:simplices, length(deaths)),))
        end
        return filter!(x -> birth(x) ≤ thresh, deaths), inf_births
    end
end

function _ripserer(
    ::Val{:involuted}, filtration, cutoff, verbose, field, dim_max, reps, implicit
)
    result = PersistenceDiagram[]
    zeroth, to_reduce, to_skip = zeroth_intervals(
        filtration, cutoff, verbose, field, _reps(reps, 0)
    )
    push!(result, zeroth)
    if dim_max > 0
        comatrix = CoboundaryMatrix{true}(field, filtration, to_reduce, to_skip)
        for dim in 1:dim_max
            columns, inf_births = compute_death_simplices!(comatrix, verbose, cutoff)
            birth_candidates = comatrix.columns_to_reduce
            matrix = BoundaryMatrix{implicit}(
                field, filtration, birth_candidates, columns;
                infinite_intervals=false,
            )
            diagram = compute_intervals!(matrix, cutoff, verbose, _reps(reps, dim))
            for birth_simplex in inf_births
                push!(
                    diagram.intervals,
                    interval(comatrix, birth_simplex, nothing, 0, _reps(reps, dim)),
                )
            end
            push!(result, diagram)
            if dim < dim_max
                comatrix = next_matrix(comatrix, verbose)
            end
        end
    end
    return result
end
