using Core.Intrinsics: box, unbox
using Base: tail

typealias RoundingFloat Union{Float64, Float32, Float16}

immutable FloatRange{T} <: Range{T}
    ref_hi::T
    ref_lo::T
    step_hi::T
    step_lo::T
    offset::Int
    length::Int

    function FloatRange(ref_hi, ref_lo, step_hi, step_lo, offset, length)
        1 <= offset <= max(length,1) || error("FloatRange: offset $offset must be in [1,$length]")
        new(ref_hi, ref_lo, step_hi, step_lo, offset, length)
    end
end

FloatRange{T<:AbstractFloat}(ref_hi::T, ref_lo::T, step_hi::T, step_lo::T, offset::Integer, length::Integer) = FloatRange{T}(ref_hi, ref_lo, step_hi, step_lo, Int(offset), Int(length))

length(r::FloatRange) = r.length
size(r::FloatRange) = (r.length,)

function linspace{T<:RoundingFloat}(start::T, stop::T, len::Integer)
    len < 2 && return _linspace1(T, start, stop, len)
    # Attempt to find exact rational approximations
    start_n, start_d = rat(start)
    stop_n, stop_d = rat(stop)
    den = lcm(start_d, stop_d)
    start_n = round(Int, den*start)
    stop_n = round(Int, den*stop)
    if T(start_n/den) == start && T(stop_n/den) == stop && T((stop_n-start_n)/den) ≈ stop-start
        return _linspace(T, start_n, stop_n, den, len)
    end
    # Find the index that returns the smallest-magnitude element
    Δ, Δfac = stop-start, 1
    if ~isfinite(Δ)   # handle overflow
        Δ, Δfac = stop/len - start/len, Int(len)
    end
    tmin = -(start/Δ)/Δfac            # interpolation t such that return value is 0
    imin = round(Int, tmin*(len-1)+1)
    if 1 < imin < len
        # The smallest-magnitude element is in the interior
        t = (imin-1)/(len-1)
        ref = T((1-t)*start + t*stop)
        # if abs(ref) < eps(max(abs(start), abs(stop)))
        #     ref = zero(T)
        # end
        step = imin-1 < len-imin ? (ref-start)/(imin-1) : (stop-ref)/(len-imin)
    elseif imin <= 1
        imin = 1
        ref = start
        step = (Δ/(len-1))*Δfac
    else
        imin = Int(len)
        ref = stop
        step = (Δ/(len-1))*Δfac
    end
    if ~isfinite(step) && len == 2
        # For very large endpoints where step overflows, use the splitting to handle the overflow
        ref_hi, step_hi = start, -start
        ref_lo, step_lo = zero(T), stop
    else
        # Double-double calculations to get high precision endpoint matching
        nb = ceil(UInt, log2(len))  # number of trailing zeros needed for exact multiplication
        step_hi = truncbits(step, nb)
        x1_hi, x1_lo = add2((1-imin)*step_hi, ref)
        x2_hi, x2_lo = add2((len-imin)*step_hi, ref)
        a, b = (start - x1_hi) - x1_lo, (stop - x2_hi) - x2_lo
        step_lo = (b - a)/(len - 1)
        ref_lo = abs(ref) < eps(max(abs(start), abs(stop))) ? zero(T) : a - (1 - imin)*step_lo
    end
    FloatRange(ref, ref_lo, step_hi, step_lo, imin, Int(len))
end

# linspace for rational numbers, start = start_n/den, stop = stop_n/den
function _linspace{T}(::Type{T}, start_n::Integer, stop_n::Integer, den::Integer, len::Integer)
    len < 2 && return _linspace1(T, start_n/den, stop_n/den, len)
    tmin = -start_n/(Float64(stop_n) - Float64(start_n))
    imin = round(Int, tmin*(len-1)+1)
    imin = clamp(imin, 1, Int(len))
    # These compute (1-t)*a and t*b separately (itp = interpolant)...
    start_itp_hi, start_itp_lo = proddiv(T, (len-imin, start_n), (den, len-1))
    stop_itp_hi, stop_itp_lo = proddiv(T, (imin-1, stop_n), (den, len-1))
    # ...and then combine them to make ref
    ref_hi, ref_lo = add2(start_itp_hi, stop_itp_hi)
    ref_hi, ref_lo = add2(ref_hi, ref_lo + (start_itp_lo + stop_itp_lo))
    # This computes step to double-double precision...
    step_hi_pre, step_lo = proddiv(T, (stop_n-start_n,), (len-1, den))
    # ...and then truncates enough low-bits of step_hi to ensure that
    # multiplication by 1:len is exact
    nb = ceil(UInt, log2(len))
    step_hi = truncbits(step_hi_pre, nb)
    step_lo += step_hi_pre - step_hi
    FloatRange(ref_hi, ref_lo, step_hi, step_lo, imin, Int(len))
end

function _linspace1{T}(::Type{T}, start, stop, len)
    len >= 0 || error("linspace($start, $stop, $len): negative length")
    if len <= 1
        len == 1 && (start == stop || error("linspace($start, $stop, $len): endpoints differ"))
        step_hi = stop-start
        return FloatRange(start, zero(T), -step_hi, (start+step_hi) - stop, 1, len)
    end
end

# promote types without risking a StackOverflowError
linspace(start::Real, stop::Real, len::Integer) = _fr(promote(start, stop)..., len)
_fr{T<:AbstractFloat}(start::T, stop::T, len) = linspace(start, stop, len)
_fr(start::BigFloat, stop::BigFloat, len) = LinSpace{BigFloat}(start, stop, len)
_fr{T<:Integer}(start::T, stop::T, len) = _linspace(Float64, start, stop, 1, len)
_fr{T<:Rational}(start::T, stop::T, len) = LinSpace{T}(start, stop, len)
_fr(start, stop, len) = error("$start::$(typeof(start)) and $stop::$(typeof(stop)) cannot be promoted to a common type")

@inline function getindex(r::FloatRange, i::Integer)
    @boundscheck checkbounds(r, i)
    unsafe_getindex(r, i)
end

# This is separate to make it useful even when running with --check-bounds=yes
function unsafe_getindex(r::FloatRange, i::Integer)
    u = i - r.offset
    shift_hi = u*r.step_hi
    shift_lo = u*r.step_lo
    x_hi, x_lo = add2(r.ref_hi, shift_hi)
    x_hi + (x_lo + (shift_lo + r.ref_lo))
end

@inline function getindex(r::FloatRange, s::OrdinalRange)
    @boundscheck checkbounds(r, s)
    sl = length(s)
    ifirst = first(s)
    ilast = last(s)
    vfirst = unsafe_getindex(r, ifirst)
    vlast  = unsafe_getindex(r, ilast)
    return linspace(vfirst, vlast, sl)
end

step(r::FloatRange) = r.step_hi + r.step_lo
first(r::FloatRange) = unsafe_getindex(r, 1)
last(r::FloatRange) = unsafe_getindex(r, length(r))

=={T<:FloatRange}(r::T, s::T) = (first(r) == first(s)) & (length(r) == length(s)) & (last(r) == last(s))

-(r::FloatRange) = FloatRange(-r.ref_hi, -r.ref_lo, -r.step_hi, -r.step_lo, r.offset, length(r))
reverse(r::FloatRange) = FloatRange(r.ref_hi, r.ref_lo, -r.step_hi, -r.step_lo, length(r)-r.offset+1, length(r))

function map{T<:AbstractFloat}(::Type{T}, r::FloatRange)
    FloatRange(T(r.ref_hi), T(r.ref_lo), T(r.step_hi), T(r.step_lo), r.offset, length(r))
end

float{T<:AbstractFloat}(r::FloatRange{T}) = r
big(r::FloatRange{BigFloat}) = r
big{T<:AbstractFloat}(r::FloatRange{T}) = linspace(big(first(r)), big(last(r)), length(r))

function .+(x::Real, r::FloatRange)
    linspace(x + first(r), x + last(r), length(r))
end
function .+(x::Number, r::FloatRange)
    linspace(x + first(r), x + last(r), length(r))
end
function .-(x::Real, r::FloatRange)
    linspace(x - first(r), x - last(r), length(r))
end
function .-(x::Number, r::FloatRange)
    linspace(x - first(r), x - last(r), length(r))
end
function .-(r::FloatRange, x::Real)
    linspace(first(r) - x, last(r) - x, length(r))
end
function .-(r::FloatRange, x::Number)
    linspace(first(r) - x, last(r) - x, length(r))
end

.*(x::Real, r::FloatRange)     = linspace(x * first(r), x * last(r), length(r))
.*(r::FloatRange, x::Real)     = x .* r

./(r::FloatRange, x::Real)     = linspace(first(r) / x, last(r) / x, length(r))

promote_rule{S,T}(::Type{FloatRange{S}}, ::Type{FloatRange{T}}) =
    FloatRange{promote_type(S,T)}
promote_rule{S,T}(::Type{FloatRange{S}}, ::Type{Base.FloatRange{T}}) =
    FloatRange{promote_type(S,T)}
promote_rule{S,O<:OrdinalRange}(::Type{FloatRange{S}}, ::Type{O}) =
    FloatRange{promote_type(S,eltype(O))}

convert{T}(::Type{FloatRange{T}}, r::Range) = linspace(T(first(r)), T(last(r)), length(r))
convert(::Type{FloatRange}, r::Range) = convert(FloatRange{eltype(r)}, r)


### Numeric utilities

function rat(x)
    y = x
    a = d = 1
    b = c = 0
    m = maxintfloat(narrow(typeof(x)))
    while abs(y) <= m
        f = trunc(Int,y)
        y -= f
        a, c = f*a + c, a
        b, d = f*b + d, b
        max(abs(a), abs(b)) <= convert(Int,m) || return c, d
        oftype(x,a)/oftype(x,b) == x && break
        y = inv(y)
    end
    return a, b
end

narrow(::Type{BigFloat}) = Float64
narrow(::Type{Float64}) = Float32
narrow(::Type{Float32}) = Float16
narrow(::Type{Float16}) = Float16

truncbits(x::Float16, nb) = box(Float16, unbox(UInt16, _truncbits(box(UInt16, unbox(Float16, x)), nb)))
truncbits(x::Float32, nb) = box(Float32, unbox(UInt32, _truncbits(box(UInt32, unbox(Float32, x)), nb)))
truncbits(x::Float64, nb) = box(Float64, unbox(UInt64, _truncbits(box(UInt64, unbox(Float64, x)), nb)))
function _truncbits{U<:Unsigned}(xi::U, nb::Integer)
    mask = typemax(U)
    xi & (mask << nb)
end

function add2{T}(u::T, v::T)
    if abs(v) > abs(u)
        return add2(v, u)
    end
    w = u + v
    w, (u-w) + v
end

function mul2(u_hi, u_lo, v::Integer)
    v == 0 && return zero(u_hi), zero(u_lo)
    nb = ceil(Int, log2(abs(v)))
    ut = truncbits(u_hi, nb)
    add2(ut*v, ((u_hi-ut) + u_lo)*v)
end

function div2(u_hi, u_lo, v::Integer)
    hi = u_hi/v
    w_hi, w_lo = mul2(hi, zero(hi), v)
    lo = (((u_hi - w_hi) - w_lo) + u_lo)/v
    add2(hi, lo)
end

@inline function proddiv(T, num, den)
    v_hi = T(num[1])
    v_hi, v_lo = _prod(v_hi, zero(T), tail(num)...)
    _div(v_hi, v_lo, den...)
end
@inline _prod(v_hi, v_lo, x, y...) = _prod(mul2(v_hi, v_lo, x)..., y...)
_prod(v_hi, v_lo) = v_hi, v_lo
@inline _div(v_hi, v_lo, x, y...) = _div(div2(v_hi, v_lo, x)..., y...)
_div(v_hi, v_lo) = v_hi, v_lo
