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
end
