using HermiteNormalForm
using LinearAlgebra
using Test

@testset "HermiteNormalForm.jl" begin
    for m in 1:20, n in 1:20
        A = rand(-big(3):3, m, n)
        H, U, rk = hnf(A)
        N = HermiteNormalForm.nullspace(A)
        @test size(N, 2) == n - rk
        # FIXME: hnf only works for matrices with all non-zero principles minors
        # by column pivoting for now.
        if rk == min(m, n) && LinearAlgebra.issuccess(lu(Rational.(A)', check=false))
            @test ishnf(H)
            @test abs(det(U.//1)) == 1
            @test iszero(A * N)
            _, _, rk = hnf(N)
            @test rk == min(size(N)...)
        end
    end
end
