using HermiteNormalForm
using LinearAlgebra
using Test

@testset "HermiteNormalForm.jl" begin
    H2 = [
         1   0   0
         6  15   0
        -1   0   0
         2   5   0
         1  -1   1
        -3   0  -2
         3  -5   4
        -4  -9  -1
         2  -3   2
        -5  -9  -2
        -2   3  -2
       ]
    @test !ishnf(H2)
    for m in 1:20, n in 1:20
        A = rand(-big(3):3, m, n)
        H, U, rk = hnf(A)
        N = HermiteNormalForm.nullspace(A)
        @test size(N, 2) == n - rk
        @test ishnf(H)
        @test abs(det(U.//1)) == 1
        @test iszero(A * N)
        _, _, rkn = hnf(N)
        @test size(HermiteNormalForm.nullspace(N), 2) == 0
        @test rkn == min(size(N)...)
    end
end
