macro prog_print(verbose, args...)
    quote
        if $(esc(verbose))
            printstyled(stderr, $(esc.(args)...); color=:green)
        end
    end
end
macro prog_println(verbose, args...)
    quote
        if $(esc(verbose))
            printstyled(stderr, $(esc.(args)...), '\n'; color=:green)
        end
    end
end

simplex_name(::Type{<:Simplex{1}}) = "edges"
simplex_name(::Type{<:Simplex{2}}) = "triangles"
simplex_name(::Type{<:Simplex{3}}) = "tetrahedra"
simplex_name(::Type{<:AbstractSimplex{D}}) where {D} = "$D-simplices"
simplex_name(::Type{<:AbstractCell{D}}) where {D} = "$D-cells"

function fmt_number(i)
    # Stolen from Humanize.jl
    value = string(i)
    group_ends = reverse(collect(length(value):-3:1))
    groups = [value[max(end_index - 2, 1):end_index] for end_index in group_ends]
    return join(groups, ",")
end
