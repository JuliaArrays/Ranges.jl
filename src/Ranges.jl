__precompile__()

module Ranges

import Base: first, last, start, next, done, step, convert, promote_rule,
             show, size, isempty, length, minimum, maximum,
             ctranspose, transpose, copy, getindex, intersect, findin, vcat,
             reverse, issorted, sort, sort!, sortperm, sum, mean, median,
             in, map, float, big
import Base: ==, -, +, .+, .-, .*, ./

include("floatrange.jl")
include("linspace.jl")

end # module
