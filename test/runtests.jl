using HermiteNormalForm
using LinearAlgebra
using Test

@testset "HermiteNormalForm.jl" begin
    for m = 1:20, n = 1:20
        A = rand(-big(3):3, m, n)
        F = hnf(A)
        if F !== nothing
            H, U = F
            @test ishnf(H)
            @test abs(det(U .// 1)) == 1
            @test LinearAlgebra.issuccess(lu(Rational.(A)', check = false))
        else
            @test !LinearAlgebra.issuccess(lu(Rational.(A)', check = false))
        end
    end

    N = 5
    d = append!(zeros(Int, 11N), -5:5)
    for _ = 1:1000
        A = rand(d, N, N)
        all(iszero, A) && continue
        NS = HermiteNormalForm.nullspace(A)
        @test iszero(size(HermiteNormalForm.nullspace(NS), 2))
        @test all(iszero, A * NS)
    end
end
