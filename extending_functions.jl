import Base.maximum
function maximum(mat::Matrix)
    maximum([maximum(vec) for vec in eachrow(mat)])
end

import Base.*
function *(a::T, b::S; type::String="Symbol") where {T, S}
    if type == "Symbol"
        return Symbol(String(a) * String(b)) 
    elseif type == "String"
        return String(a) * String(b)
    end 
end
function *(a::Vector{T}; delim=:none) where T
    if delim == :none
        return join(a) |> T
    else
        a = reshape(a, 1, :)
        b = [delim for _ in a]
        b = b .|> T
        c = vcat(a, b)[:][1:end-1]
        return *(c)
    end
end

# adds functionality to groupby, by also giving an array of the keys, used for for loops
function grouping(df, variables=[:TOK, :DATE, :SHOT, :TIME])
	gbdf = groupby(df, variables)
	return gbdf, keys(gbdf) |> collect |> sort
end