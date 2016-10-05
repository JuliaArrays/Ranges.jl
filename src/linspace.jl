# Standalone https://github.com/JuliaLang/julia/pull/18777 for use
# with earlier versions of Julia

import Base: first, last, start, next, done, step, convert, promote_rule,
             show, size, isempty, length, minimum, maximum,
             ctranspose, transpose, copy, getindex, intersect, findin, vcat,
             reverse, issorted, sort, sort!, sortperm, sum, mean, median,
             in, map, float, big
import Base: ==, -, +, .+, .-, .*, ./

immutable LinSpace{T} <: Range{T}
    start::T
    stop::T
    len::Int
    lendiv::Int

    function LinSpace(start,stop,len)
        len >= 0 || error("linspace($start, $stop, $len): negative length")
        if len == 1
            start == stop || error("linspace($start, $stop, $len): endpoints differ")
            return new(start, stop, 1, 1)
        end
        new(start,stop,len,max(len-1,1))
    end
end

function LinSpace(start, stop, len::Integer)
    T = typeof((stop-start)/len)
    LinSpace{T}(start, stop, len)
end

linspace(start, stop, len::Real=50) = LinSpace(start, stop, Int(len))

range(start, step, len) = linspace(start, fma(len-1, step, start), len)

function show(io::IO, r::LinSpace)
    print(io, "Ranges.linspace(")
    show(io, first(r))
    print(io, ',')
    show(io, last(r))
    print(io, ',')
    show(io, length(r))
    print(io, ')')
end

logspace(start::Real, stop::Real, n::Integer=50) = 10.^linspace(start, stop, n)

## interface implementations

isempty(r::LinSpace) = length(r) == 0

step(r::LinSpace) = (last(r)-first(r))/r.lendiv

length(r::LinSpace) = r.len

first(r::LinSpace) = r.start

last(r::LinSpace) = r.stop

start(r::LinSpace) = 1
done(r::LinSpace, i::Int) = length(r) < i
@inline function next(r::LinSpace, i::Int)
    unsafe_getindex(r, i), i+1
end

@inline function getindex(r::LinSpace, i::Integer)
    @boundscheck checkbounds(r, i)
    unsafe_getindex(r, i)
end

# This is separate to make it useful even when running with --check-bounds=yes
function unsafe_getindex(r::LinSpace, i::Integer)
    d = r.lendiv
    j, a, b = ifelse(2i >= length(r), (i-1, r.start, r.stop), (length(r)-i, r.stop, r.start))
    lerpi(j, d, a, b)
end

@inline function lerpi{T}(j::Integer, d::Integer, a::T, b::T)
    t = j/d
    # computes (1-t)*a + t*b
    # T(fma(t, b, fma(-t, a, a)))
    # Until issues of precision, etc, are worked out in
    # https://github.com/JuliaLang/julia/pull/18777, for safety use
    # the simplest approach
    convert(T, (1-t)*a + t*b)
end

@inline function getindex{T}(r::LinSpace{T}, s::OrdinalRange)
    @boundscheck checkbounds(r, s)
    sl = length(s)
    ifirst = first(s)
    ilast = last(s)
    vfirst = unsafe_getindex(r, ifirst)
    vlast  = unsafe_getindex(r, ilast)
    return linspace(vfirst, vlast, sl)
end

=={T<:LinSpace}(r::T, s::T) = (first(r) == first(s)) & (length(r) == length(s)) & (last(r) == last(s))

-(r::LinSpace) = LinSpace(-r.start, -r.stop, r.len)

function .+(x::Real, r::LinSpace)
    LinSpace(x + r.start, x + r.stop, r.len)
end
function .+(x::Number, r::LinSpace)
    LinSpace(x + r.start, x + r.stop, r.len)
end
function .+{T}(x::Ref{T}, r::LinSpace{T})
    LinSpace(x + r.start, x + r.stop, r.len)
end
function .-(x::Real, r::LinSpace)
    LinSpace(x - r.start, x - r.stop, r.len)
end
function .-(x::Number, r::LinSpace)
    LinSpace(x - r.start, x - r.stop, r.len)
end
function .-{T}(x::Ref{T}, r::LinSpace{T})
    LinSpace(x - r.start, x - r.stop, r.len)
end
function .-(r::LinSpace, x::Real)
    LinSpace(r.start - x, r.stop - x, r.len)
end
function .-(r::LinSpace, x::Number)
    LinSpace(r.start - x, r.stop - x, r.len)
end
function .-{T}(r::LinSpace{T}, x::Ref{T})
    LinSpace(r.start - x, r.stop - x, r.len)
end

.*(x::Real, r::LinSpace)     = LinSpace(x * r.start, x * r.stop, r.len)
.*(r::LinSpace, x::Real)     = x .* r

./(r::LinSpace, x::Real)     = LinSpace(r.start / x, r.stop / x, r.len)

promote_rule{T1,T2}(::Type{LinSpace{T1}},::Type{LinSpace{T2}}) =
    LinSpace{promote_type(T1,T2)}
convert{T<:AbstractFloat}(::Type{LinSpace{T}}, r::LinSpace{T}) = r
convert{T<:AbstractFloat}(::Type{LinSpace{T}}, r::LinSpace) =
    LinSpace{T}(r.start, r.stop, r.len)

promote_rule{F,OR<:OrdinalRange}(::Type{LinSpace{F}}, ::Type{OR}) =
    LinSpace{promote_type(F,eltype(OR))}
convert{T<:AbstractFloat}(::Type{LinSpace{T}}, r::OrdinalRange) =
    linspace(convert(T, first(r)), convert(T, last(r)), length(r))
convert{T}(::Type{LinSpace}, r::OrdinalRange{T}) =
    convert(LinSpace{typeof(float(first(r)))}, r)

# Promote FloatRange to LinSpace
promote_rule{F,OR<:FloatRange}(::Type{LinSpace{F}}, ::Type{OR}) =
    LinSpace{promote_type(F,eltype(OR))}
convert{T<:AbstractFloat}(::Type{LinSpace{T}}, r::FloatRange) =
    linspace(convert(T, first(r)), convert(T, last(r)), length(r))
convert{T<:AbstractFloat}(::Type{LinSpace}, r::FloatRange{T}) =
    convert(LinSpace{T}, r)

reverse(r::LinSpace)     = LinSpace(r.stop, r.start, length(r))

# Extra methods
function map{T<:AbstractFloat}(::Type{T}, r::LinSpace)
    LinSpace(T(r.start), T(r.stop), length(r))
end

float(r::LinSpace) = LinSpace(float(r.start), float(r.stop), length(r))
big(r::LinSpace) = LinSpace(big(r.start), big(r.stop), length(r))

function lerpi(j::Integer, d::Integer, A::AbstractArray, B::AbstractArray)
    broadcast((a,b) -> lerpi(j, d, a, b), A, B)
end
function lerpi(j::Integer, d::Integer, a::BigFloat, b::BigFloat)
    t = BigFloat(j)/d
    fma(t, b, fma(-t, a, a))
end
function lerpi(j::Integer, d::Integer, a::Rational, b::Rational)
    ((d-j)*a)/d + (j*b)/d
end

function show(io::IO, ::MIME"text/plain", r::LinSpace)
    # show for linspace, e.g.
    # linspace(1,3,7)
    # 7-element LinSpace{Float64}:
    #   1.0,1.33333,1.66667,2.0,2.33333,2.66667,3.0
    print(io, summary(r))
    if !isempty(r)
        println(io, ":")
        Base.print_range(io, r)
    end
end
