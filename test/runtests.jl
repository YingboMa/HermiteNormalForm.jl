using HermiteNormalForm
using LinearAlgebra
using Test

@testset "HermiteNormalForm.jl" begin
    for m in 1:20, n in 1:20
        A = rand(-big(3):3, m, n)
        F = hnf(A)
        if F !== nothing
            H, U = F
            @test ishnf(H)
            @test abs(det(U.//1)) == 1
        else
            @test !LinearAlgebra.issuccess(lu(Rational.(A)', check=false))
        end
    end
end
