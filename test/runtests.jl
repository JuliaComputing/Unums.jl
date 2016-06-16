using Unums
using Base.Test

import Unums: Unum00, Ubound00, Ubound22, Ubound34


for U in (Ubound00, Ubound22, Ubound34, Unum00)
    @test U(0) + U(1) == U(1)
    @test U(0) - U(1) == U(-1)
    @test U(0) + U(2) == U(2)
    @test U(0) - U(2) == U(-2)
    @test U(1) + U(1) == U(2)
    @test U(1) - U(1) == U(0)
    @test U(1) + U(2) == U(3)
    @test U(1) - U(2) == U(-1)
    @test U(0) + U(0.5) == U(0.5)
    @test U(0) - U(0.5) == U(-0.5)
    @test U(1) + U(0.5) == U(1.5)
    @test U(1) - U(0.5) == U(0.5)

    @test U(1)*U(1) == U(1)
    @test U(1)*U(2) == U(2)
    @test U(1)/U(2) == U(0.5)
    @test U(1)/U(-2) == U(-0.5)
    @test sqrt(U(1)) == U(1)
    @test sqrt(U(2)) == U(sqrt(2))
    @test U(2)^4 == U(2^4)

    if U <: Unums.Ubound
        @test U(-1,1)^2 == U(0,1)
    end
end


for U in (Unum22, Ubound22)
    u = convert(U, 2.1)
    lo,hi,olo,ohi = limits(u)
    @test lo == 2
    @test hi == 2.125
    @test olo == true
    @test ohi == true
end

@test Unum22(2) in Unum22(2)
@test Unum22(2) in Ubound22(2,3)
@test !(Unum22(2) in Unum22(2.1))

@test isfinite(Unum22(2))
@test isfinite(Unum22(2.1))
@test !isfinite(Unum22(1e300))
