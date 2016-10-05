# This file is a part of Julia. License is MIT: http://julialang.org/license

# Test whether r[i] has full precision
function test_linspace{T<:AbstractFloat}(r::Ranges.LinSpace{T})
    isempty(r) && return nothing
    n = length(r)
    f, l = first(r), last(r)
    a, b, d = BigFloat(f), BigFloat(l), max(n-1,1)
    Δ = 2.1*max(eps(f), eps(l))  # FIXME: tighten these bounds
    for i = 1:n
        c = b*(BigFloat(i-1)/d) + a*(BigFloat(d-i+1)/d)
        @test abs(r[i]-T(c)) <= Δ
    end
    nothing
end

L32 = Ranges.linspace(Int32(1), Int32(4), 4)
L64 = Ranges.linspace(Int64(1), Int64(4), 4)
@test L32[1] == 1 && L64[1] == 1
@test L32[2] == 2 && L64[2] == 2
@test L32[3] == 3 && L64[3] == 3
@test L32[4] == 4 && L64[4] == 4

test_linspace(Ranges.linspace(0.1,0.3,3))
test_linspace(Ranges.linspace(0.0,0.3,4))
test_linspace(Ranges.linspace(0.3,-0.1,5))
test_linspace(Ranges.linspace(0.1,-0.3,5))
test_linspace(Ranges.linspace(0.0,1.0,11))
test_linspace(Ranges.linspace(0.0,1.0,0))
test_linspace(Ranges.linspace(0.0,-1.0,0))
test_linspace(Ranges.linspace(0.0,-1.0,11))
test_linspace(Ranges.linspace(1.0,27.0,1275))
test_linspace(Ranges.linspace(0.0,2.1,4))
test_linspace(Ranges.linspace(0.0,3.3,4))
test_linspace(Ranges.linspace(0.1,3.4,4))
test_linspace(Ranges.linspace(0.0,3.9,4))
test_linspace(Ranges.linspace(0.1,4.0,4))
test_linspace(Ranges.linspace(1.1,3.3,3))
test_linspace(Ranges.linspace(0.3,1.1,9))
test_linspace(Ranges.linspace(0.0,0.0,1))
test_linspace(Ranges.linspace(0.0,0.0,1))

for T = (Float32, Float64,),# BigFloat),
    a = -5:25, s = [-5:-1;1:25;], d = 1:25, n = -1:15
    den   = convert(T,d)
    start = convert(T,a)/den
    stop  = convert(T,(a+(n-1)*s))/den
    test_linspace(Ranges.linspace(start, stop, max(n,0)))
end

# linspace & ranges with very small endpoints
for T = (Float32, Float64)
    z = zero(T)
    u = eps(z)
    @test first(Ranges.linspace(u,u,0)) == u
    @test last(Ranges.linspace(u,u,0)) == u
    @test first(Ranges.linspace(-u,u,0)) == -u
    @test last(Ranges.linspace(-u,u,0)) == u
    @test [Ranges.linspace(-u,u,0);] == []
    @test [Ranges.linspace(-u,-u,1);] == [-u]
    @test [Ranges.linspace(-u,u,2);] == [-u,u]
    @test [Ranges.linspace(-u,u,3);] == [-u,0,u]
    @test first(Ranges.linspace(-u,-u,0)) == -u
    @test last(Ranges.linspace(-u,-u,0)) == -u
    @test first(Ranges.linspace(u,-u,0)) == u
    @test last(Ranges.linspace(u,-u,0)) == -u
    @test [Ranges.linspace(u,-u,0);] == []
    @test [Ranges.linspace(u,u,1);] == [u]
    @test [Ranges.linspace(u,-u,2);] == [u,-u]
    @test [Ranges.linspace(u,-u,3);] == [u,0,-u]
    v = [Ranges.linspace(-u,u,12);]
    @test length(v) == 12
end

# linspace with very large endpoints
for T = (Float32, Float64)
    a = realmax()
    for i = 1:5
        @test [Ranges.linspace(a,a,1);] == [a]
        @test [Ranges.linspace(-a,-a,1);] == [-a]
        b = realmax()
        for j = 1:5
            @test [Ranges.linspace(-a,b,0);] == []
            @test [Ranges.linspace(-a,b,2);] == [-a,b]
            @test [Ranges.linspace(-a,b,3);] == [-a,(b-a)/2,b]
            @test [Ranges.linspace(a,-b,0);] == []
            @test [Ranges.linspace(a,-b,2);] == [a,-b]
            @test [Ranges.linspace(a,-b,3);] == [a,(a-b)/2,-b]
            for c = maxintfloat(T)-3:maxintfloat(T)
                s = Ranges.linspace(-a,b,c)
                @test first(s) == -a
                @test last(s) == b
                c <= typemax(Int) && @test length(s) == c
                @test s.len == c
                s = Ranges.linspace(a,-b,c)
                @test first(s) == a
                @test last(s) == -b
                c <= typemax(Int) && @test length(s) == c
                @test s.len == c
            end
            b = prevfloat(b)
        end
        a = prevfloat(a)
    end
end

# linspace with 1 or 0 elements (whose step length is NaN)
@test issorted(Ranges.linspace(1,1,0))
@test issorted(Ranges.linspace(1,1,1))

# comparing and hashing ranges
let
    Rs = Range[1:2, map(Int32,1:3:17), map(Int64,1:3:17), 1:0, 17:-3:0,
               0.0:0.1:1.0, map(Float32,0.0:0.1:1.0),
               Ranges.linspace(0, 1, 20), map(Float32, Ranges.linspace(0, 1, 20))]
    for r in Rs
        ar = collect(r)
        @test r != ar
        @test !isequal(r,ar)
        for s in Rs
            as = collect(s)
            @test !isequal(r,s) || hash(r)==hash(s)
            @test (r==s) == (ar==as)
        end
    end
end


r = Ranges.linspace(1/3,5/7,6)
@test length(r) == 6
@test r[1] == 1/3
@test abs(r[end] - 5/7) <= eps(5/7)
r = Ranges.linspace(0.25,0.25,1)
@test length(r) == 1
@test_throws ErrorException Ranges.linspace(0.25,0.5,1)

# Issue #11245
let io = IOBuffer()
    show(io, Ranges.linspace(1, 2, 3))
    str = takebuf_string(io)
    @test str == "Ranges.linspace(1.0,2.0,3)"
end

# stringmime/show should display the range or linspace nicely
# to test print_range in range.jl
replstrmime(x) = sprint((io,x) -> show(IOContext(io, limit=true), MIME("text/plain"), x), x)
@test stringmime("text/plain", Ranges.linspace(1,5,7)) == "7-element Ranges.LinSpace{Float64}:\n 1.0,1.66667,2.33333,3.0,3.66667,4.33333,5.0"
@test repr(Ranges.linspace(1,5,7)) == "Ranges.linspace(1.0,5.0,7)"
# next is to test a very large range, which should be fast because print_range
# only examines spacing of the left and right edges of the range, sufficient
# to cover the designated screen size.
@test replstrmime(Ranges.linspace(0,100, 10000)) == "10000-element Ranges.LinSpace{Float64}:\n 0.0,0.010001,0.020002,0.030003,0.040004,…,99.95,99.96,99.97,99.98,99.99,100.0"

@test sprint(io -> show(io,UnitRange(1,2))) == "1:2"
@test sprint(io -> show(io,StepRange(1,2,5))) == "1:2:5"


# Issue 11049 and related
@test promote(Ranges.linspace(0f0, 1f0, 3), Ranges.linspace(0., 5., 2)) ===
    (Ranges.linspace(0., 1., 3), Ranges.linspace(0., 5., 2))
@test convert(Ranges.LinSpace{Float64}, Ranges.linspace(0., 1., 3)) === Ranges.linspace(0., 1., 3)
@test convert(Ranges.LinSpace{Float64}, Ranges.linspace(0f0, 1f0, 3)) === Ranges.linspace(0., 1., 3)

@test promote(Ranges.linspace(0., 1., 3), 0:5) === (Ranges.linspace(0., 1., 3),
                                             Ranges.linspace(0., 5., 6))
@test convert(Ranges.LinSpace{Float64}, 0:5) === Ranges.linspace(0., 5., 6)
@test convert(Ranges.LinSpace{Float64}, 0:1:5) === Ranges.linspace(0., 5., 6)
@test convert(Ranges.LinSpace, 0:5) === Ranges.linspace(0., 5., 6)
@test convert(Ranges.LinSpace, 0:1:5) === Ranges.linspace(0., 5., 6)

function test_range_index(r, s)
    @test typeof(r[s]) == typeof(r)
    @test [r;][s] == [r[s];]
end
test_range_index(Ranges.linspace(0.1, 0.3, 3), 1:2)
test_range_index(Ranges.linspace(0.1, 0.3, 3), 1:0)
test_range_index(Ranges.linspace(1.0, 1.0, 1), 1:1)
test_range_index(Ranges.linspace(1.0, 1.0, 1), 1:0)
test_range_index(Ranges.linspace(1.0, 2.0, 0), 1:0)

function test_linspace_identity{T}(r::Ranges.LinSpace{T}, mr::Ranges.LinSpace{T})
    @test -r == mr
    @test -collect(r) == collect(mr)
    @test isa(-r, Ranges.LinSpace)

    @test 1 + r + (-1) == r
    @test isa(1 + r + (-1), Ranges.LinSpace)
    @test 1 - r - 1 == mr
    @test isa(1 - r - 1, Ranges.LinSpace)

    @test 1 * r * 1 == r
    @test 2 * r * T(0.5) == r
    @test isa(1 * r * 1, Ranges.LinSpace)
    @test r / 1 == r
    @test r / 2 * 2 == r
    @test r / T(0.5) * T(0.5) == r
    @test isa(r / 1, Ranges.LinSpace)

    @test (2 * collect(r) == collect(r * 2) == collect(2 * r) ==
           collect(r * T(2.0)) == collect(T(2.0) * r) ==
           collect(r / T(0.5)) == -collect(mr * T(2.0)))
end

test_linspace_identity(Ranges.linspace(1.0, 27.0, 10), Ranges.linspace(-1.0, -27.0, 10))
test_linspace_identity(Ranges.linspace(1f0, 27f0, 10), Ranges.linspace(-1f0, -27f0, 10))

test_linspace_identity(Ranges.linspace(1.0, 27.0, 0), Ranges.linspace(-1.0, -27.0, 0))
test_linspace_identity(Ranges.linspace(1f0, 27f0, 0), Ranges.linspace(-1f0, -27f0, 0))

test_linspace_identity(Ranges.linspace(1.0, 1.0, 1), Ranges.linspace(-1.0, -1.0, 1))
test_linspace_identity(Ranges.linspace(1f0, 1f0, 1), Ranges.linspace(-1f0, -1f0, 1))

@test reverse(Ranges.linspace(1.0, 27.0, 1275)) == Ranges.linspace(27.0, 1.0, 1275)
@test [reverse(Ranges.linspace(1.0, 27.0, 1275));] ==
    reverse([Ranges.linspace(1.0, 27.0, 1275);])

# PR 12200 and related
for _r in (1:2:100, 1:100, 1f0:2f0:100f0, 1.0:2.0:100.0,
           Ranges.linspace(1, 100, 10), Ranges.linspace(1f0, 100f0, 10))
    float_r = float(_r)
    big_r = big(_r)
    @test typeof(big_r).name === typeof(_r).name
    if eltype(_r) <: AbstractFloat
        @test isa(float_r, typeof(_r))
        @test eltype(big_r) === BigFloat
    else
        @test isa(float_r, Range)
        @test eltype(float_r) <: AbstractFloat
        @test eltype(big_r) === BigInt
    end
end

@test_throws DimensionMismatch Ranges.linspace(1.,5.,5) + Ranges.linspace(1.,5.,6)
@test_throws DimensionMismatch Ranges.linspace(1.,5.,5) - Ranges.linspace(1.,5.,6)
@test_throws DimensionMismatch Ranges.linspace(1.,5.,5) .* Ranges.linspace(1.,5.,6)
@test_throws DimensionMismatch Ranges.linspace(1.,5.,5) ./ Ranges.linspace(1.,5.,6)

# linspace of other types
r = Ranges.linspace(0, 3//10, 4)
@test eltype(r) == Rational{Int}
@test r[2] === 1//10

a, b = 1.0, nextfloat(1.0)
ba, bb = BigFloat(a), BigFloat(b)
r = Ranges.linspace(ba, bb, 3)
@test eltype(r) == BigFloat
@test r[1] == a && r[3] == b
@test r[2] == (ba+bb)/2

a, b = rand(10), rand(10)
ba, bb = big(a), big(b)
r = Ranges.linspace(a, b, 5)
@test r[1] == a && r[5] == b
for i = 2:4
    x = ((5-i)//4)*ba + ((i-1)//4)*bb
    @test r[i] ≈ Float64.(x)
end
