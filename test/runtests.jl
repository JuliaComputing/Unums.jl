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
end
