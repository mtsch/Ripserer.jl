function prog_print(progress, args...)
    return progress && printstyled(stderr, args...; color=:green)
end
function prog_println(progress, args...)
    return prog_print(progress, args..., '\n')
end

simplex_name(::Type{<:Simplex{1}}) = "edges"
simplex_name(::Type{<:Simplex{2}}) = "triangles"
simplex_name(::Type{<:Simplex{3}}) = "tetrahedra"
simplex_name(::Type{<:AbstractSimplex{D}}) where {D} = "$D-simplices"

function fmt_number(i)
    # Stolen from Humanize.jl
    value = string(i)
    group_ends = reverse(collect(length(value):-3:1))
    groups = [value[max(end_index - 2, 1):end_index] for end_index in group_ends]
    return join(groups, ",")
end
