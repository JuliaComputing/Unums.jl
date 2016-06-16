export Bnum, Bbound

immutable Bnum <: Real
    num::BigFloat
    open::Bool
end

immutable Bbound <: Real
    lo::Bnum
    hi::Bnum
end


Bnum(x::Real,open::Bool) = Bnum(BigFloat(x),open)
Bnum(x::Real) = Bnum(x,false)

convert(::Type{Bnum},x::Bnum) = x
convert(::Type{Bnum},x::Real) = Bnum(x)


import Base.MPFR: to_mpfr, ClongMax

-(x::Bnum) = Bnum(-x.num,x.open)


function +(x::Bnum, y::Bnum, r::RoundingMode)
    znum = BigFloat()
    inex = ccall((:mpfr_add,:libmpfr), Int32, (Ptr{BigFloat}, Ptr{BigFloat}, Ptr{BigFloat}, Int32), &znum, &x.num, &y.num, to_mpfr(r)) != 0

    zopen = (x.open | y.open | inex) & (isfinite(x.num) | x.open) & (isfinite(y.num) | y.open)
    Bnum(znum, zopen)
end

function -(x::Bnum, y::Bnum, r::RoundingMode)
    znum = BigFloat()
    inex = ccall((:mpfr_sub,:libmpfr), Int32, (Ptr{BigFloat}, Ptr{BigFloat}, Ptr{BigFloat}, Int32), &znum, &x.num, &y.num, to_mpfr(r)) != 0

    zopen = (x.open | y.open | inex) & (isfinite(x.num) | x.open) & (isfinite(y.num) | y.open)
    Bnum(znum, zopen)
end


function *(x::Bnum, y::Bnum, r::RoundingMode)
    znum = BigFloat()
    inex = ccall((:mpfr_mul,:libmpfr), Int32, (Ptr{BigFloat}, Ptr{BigFloat}, Ptr{BigFloat}, Int32), &znum, &x.num, &y.num, to_mpfr(r)) != 0

    zopen = (x.open | y.open | inex) & (isfinite(x.num) & (x.num!=0) | x.open) & (isfinite(y.num) & (y.num!=0) | y.open)
    Bnum(znum, zopen)
end

function /(x::Bnum, y::Bnum, r::RoundingMode)
    znum = BigFloat()
    inex = ccall((:mpfr_div,:libmpfr), Int32, (Ptr{BigFloat}, Ptr{BigFloat}, Ptr{BigFloat}, Int32), &znum, &x.num, &y.num, to_mpfr(r)) != 0

    zopen = (x.open | y.open | inex) & (isfinite(x.num) & (x.num!=0) | x.open) & (!isfinite(y.num) | y.open)
    Bnum(znum, zopen)
end

function sqrt(x::Bnum, r::RoundingMode)
    znum = BigFloat()
    inex = ccall((:mpfr_sqrt,:libmpfr), Int32, (Ptr{BigFloat}, Ptr{BigFloat}, Int32), &znum, &x.num, to_mpfr(r)) != 0

    zopen = (x.open | inex)
    Bnum(znum, zopen)
end

function ^(x::Bnum, y::ClongMax, r::RoundingMode)
    znum = BigFloat()
    inex = ccall((:mpfr_pow_si, :libmpfr), Int32, (Ptr{BigFloat}, Ptr{BigFloat}, Clong, Int32), &znum, &x.num, y, to_mpfr(r)) != 0

    zopen = y != 0 && (x.open | inex)
    Bnum(znum, zopen)
end



import Base: min, max
# Technically not correct as we need to know the orientation, but assumes
# min is being used for lower bound, and max for upper bound
function min(x::Bnum,y::Bnum)
    if x.num == y.num
        Bnum(x.num,x.open & y.open)
    elseif x.num < y.num
        x
    else
        y
    end
end
function max(x::Bnum,y::Bnum)
    if x.num == y.num
        Bnum(x.num,x.open & y.open)
    elseif x.num > y.num
        x
    else
        y
    end
end

import Base: isnan
isnan(x::Bnum) = isnan(x.num)


function convert(::Type{Bbound},x::Real)
    c = convert(Bnum,x)
    Bbound(c,c)
end

function print(io::IO, b::Bbound)
    nlo = b.lo.num
    flo = Float64(nlo)
    if flo == 0
        # don't print signed zero
        flo = abs(flo)
    end
    
    nhi = b.hi.num
    fhi = Float64(nhi)
    if fhi == 0
        fhi = abs(fhi)
    end
    
    if nlo == nhi
        print(io,nlo==flo?flo:nlo)
    else
        print(io,b.lo.open?'(':'[',nlo==flo?flo:nlo,',',nhi==fhi?fhi:nhi,b.hi.open?')':']')
    end
end

show(io::IO, b::Bbound) = print(io,typeof(b),'\n',b)


-(x::Bbound) = Bbound(-x.hi,-x.lo)

function +(x::Bbound, y::Bbound)
    Bbound(+(x.lo,y.lo,RoundDown), +(x.hi,y.hi,RoundUp))
end
function -(x::Bbound, y::Bbound)
    Bbound(-(x.lo,y.hi,RoundDown), -(x.hi,y.lo,RoundUp))
end


function *(x::Bbound, y::Bbound)
    if isposz(x) # x >= 0
        if isposz(y)
            Bbound(*(x.lo,y.lo,RoundDown),*(x.hi,y.hi,RoundUp))
        elseif isnegz(y)
            Bbound(*(x.hi,y.lo,RoundDown),*(x.lo,y.hi,RoundUp))
        else
            Bbound(*(x.hi,y.lo,RoundDown),*(x.hi,y.hi,RoundUp))
        end
    elseif isnegz(x) # x <= 0
        if isposz(y)
            Bbound(*(x.lo,y.hi,RoundDown),*(x.hi,y.lo,RoundUp))
        elseif isnegz(y)
            Bbound(*(x.hi,y.hi,RoundDown),*(x.lo,y.lo,RoundUp))
        else
            Bbound(*(x.lo,y.hi,RoundDown),*(x.lo,y.lo,RoundUp))
        end
    else
        if isposz(y)
            Bbound(*(x.lo,y.hi,RoundDown),*(x.hi,y.hi,RoundUp))
        elseif isnegz(y)
            Bbound(*(x.hi,y.lo,RoundDown),*(x.lo,y.lo,RoundUp))
        else
            Bbound(min(*(x.lo,y.hi,RoundDown),*(x.hi,y.lo,RoundDown)),
                     max(*(x.lo,y.lo,RoundUp),*(x.hi,y.hi,RoundUp)))
        end
    end
end


function /(x::Bbound, y::Bbound)
    if ispos(y) # b strictly positive
        if isposz(x)
            Bbound(/(x.lo,y.hi,RoundDown), /(x.hi,y.lo,RoundUp))
        elseif isnegz(x)
            Bbound(/(x.lo,y.lo,RoundDown), /(x.hi,y.hi,RoundUp))
        else
            Bbound(/(x.lo,y.lo,RoundDown), /(x.hi,y.lo,RoundUp))
        end
    elseif isneg(y)
        if isposz(x)
            Bbound(/(x.hi,y.hi,RoundDown), /(x.lo,y.lo,RoundUp))
        elseif isnegz(x)
            Bbound(/(x.hi,y.lo,RoundDown), /(x.lo,y.hi,RoundUp))
        else
            Bbound(/(x.hi,y.hi,RoundDown), /(x.lo,y.hi,RoundUp))
        end
    else
        Bbound(NaN)
    end
end


function abs(x::Bbound)
    if isposz(x)
        x
    elseif isnegz(x)
        -x
    else
        Bbound(zero(BigFloat),max(-x.lo,x.hi))
    end
end
        

function sqrt(x::Bbound)
    if isposz(x)
        Bbound(sqrt(x.lo,RoundDown),sqrt(x.hi,RoundUp))
    else
        Bbound(NaN)
    end
end

function ^(x::Bbound,y::Integer)
    if y == 0
        Bbound(1)
    elseif y > 0
        if isodd(y) || isposz(x)
            Bbound(^(x.lo,y,RoundDown),^(x.hi,y,RoundUp))
        elseif isnegz(x)
            Bbound(^(x.hi,y,RoundDown),^(x.lo,y,RoundUp))
        else
            Bbound(zero(BigFloat),^(max(-x.lo,x.hi),y,RoundUp))
        end
    else
        if isodd(y) || isposz(x)
            Bbound(^(x.hi,y,RoundDown),^(x.lo,y,RoundUp))
        elseif isnegz(x)
            Bbound(^(x.lo,y,RoundDown),^(x.hi,y,RoundUp))
        else
            Bbound(NaN)
        end            
    end
end


isposz(x::Bbound) = x.lo.num >= 0
isnegz(x::Bbound) = x.hi.num <= 0
ispos(x::Bbound) = x.lo.num > 0 || x.lo.num == 0 && x.lo.open
isneg(x::Bbound) = x.hi.num < 0 || x.hi.num == 0 && x.hi.open

==(x::Bnum,y::Bnum) = x.num == y.num && x.open == y.open
==(x::Bbound,y::Bbound) = x.lo == y.lo && x.hi == y.hi

<(x::Bbound,y::Bbound) = x.hi.num < y.lo.num || x.hi.num == y.lo.num && (x.hi.open | y.lo.open)
<=(x::Bbound,y::Bbound) = x.hi.num <= y.lo.num

function in(x::Bbound,y::Bbound)
    (y.lo.num < x.lo.num || (y.lo.num == x.lo.num && (!y.lo.open || x.lo.open))) &&
    (x.hi.num < y.hi.num || (x.hi.num == y.hi.num && (!y.hi.open || x.hi.open)))
end
