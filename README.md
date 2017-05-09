# Ranges

[![Build Status](https://travis-ci.org/JuliaArrays/Ranges.jl.svg?branch=master)](https://travis-ci.org/JuliaArrays/Ranges.jl)

[![codecov.io](http://codecov.io/github/JuliaArrays/Ranges.jl/coverage.svg?branch=master)](http://codecov.io/github/JuliaArrays/Ranges.jl?branch=master)

This package exists to support more general Ranges for Julia 0.5. This package is not necessary on later Julia versions due to https://github.com/JuliaLang/julia/pull/18777.

## LinSpace

Currently this package contains just one range, `LinSpace`. This is a different implementation than the one in Julia proper; its main advantage is that it supports any type supporting "linear interpolation,"

```jl
x = (1-t)*a + t*b
```

Here, `t` is a number between 0 and 1, and `a` and `b` are the
endpoints of the range.  While `a` and `b` might be numbers, this
supports any type for which `*` and `+` return an object of the same
type as `a` and `b`. While this formula can be evaluated for any `t`,
a `LinSpace` also encodes a length `n` and can be evaluated only for
indexes `i = 1, 2, ..., n` (for which `t = (i-1)/(n-1)`).

The interface is the same as for all other Ranges. Aside from indexing
and iteration, the main functions are `Ranges.linspace` (to construct
the object) and `step` (to extract the step between adjacent values).

Note that the LinSpace type is not equivalent to the version of LinSpace present in Julia 0.6.
