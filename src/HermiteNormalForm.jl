module HermiteNormalForm

using LinearAlgebra
export hnf!, hnf, ishnf

"""
    ishnf(A)

Check if matrix `A` is in Hermite norm form.
"""
function ishnf(A)
    m, n = size(A)
    for j = 1:n
        for i = 1:min(m, j - 1)
            A[i, j] == 0 || return false
        end
        if m >= j
            A[j, j] > 0 || return false
        end
    end
    for i = 1:min(m, n)
        p = A[i, i]
        for j = 1:min(n, i - 1)
            0 <= A[i, j] < p || return false
        end
    end
    return true
end

"""
    H, U = hnf(A)

`U` is unimodular and `H` is the Hermite normal form of `A`. If `A` is singular,
then it returns `nothing`.
"""
hnf(A) = hnf!(copy(A))

function hnf!(A)
    m, n = size(A)
    T = eltype(A)
    U = Matrix{T}(I, n, n)
    zz = zero(T)
    for k = 1:m
        # Pivot: A[k, k] should not be 0
        if n >= k && A[k, k] == zz
            pivot = k
            for j = k+1:n
                if A[k, j] != zz
                    pivot = j
                end
            end
            pivot == k && return nothing # rank deficient
            for i = 1:m
                A[i, k], A[i, pivot] = A[i, pivot], A[i, k]
            end
            for i = 1:n
                U[i, k], U[i, pivot] = U[i, pivot], U[i, k]
            end
        end
        # [A11   0 ] k-1 rows
        # [A21  A22]
        # Zero out A[k, k+1:n] === A22[1, 2:end] by multiplying
        # [p  -A[k, j]/d]
        # [q   A[k, k]/d]
        for j = k+1:n
            Akk, Akj = A[k, k], A[k, j]
            d, p, q = gcdx(Akk, Akj)
            Akkd, Akjd = div(Akk, d), div(Akj, d)
            for i = 1:m
                Aik, Aij = A[i, k], A[i, j]
                A[i, k] = Aik * p + Aij * q
                A[i, j] = -Aik * Akjd + Aij * Akkd
            end
            for i = 1:n
                Uik, Uij = U[i, k], U[i, j]
                U[i, k] = Uik * p + Uij * q
                U[i, j] = -Uik * Akjd + Uij * Akkd
            end
        end
        n >= k || continue
        # Ensure the positivity of A[k, k]
        if A[k, k] < zz
            @. A[:, k] = -A[:, k]
            @. U[:, k] = -U[:, k]
        end
        # Minimize A[k, 1:k-1] === A21[1, :]
        for j = 1:k-1
            mul = fld(A[k, j], A[k, k])
            @. A[:, j] -= mul * A[:, k]
            @. U[:, j] -= mul * U[:, k]
        end
    end
    A, U
end

macro vp(expr)
    nodes = (Symbol("llvm.loop.vectorize.predicate.enable"), 1)
    if expr.head != :for
        error("Syntax error: loopinfo needs a for loop")
    end
    push!(expr.args[2].args, Expr(:loopinfo, nodes))
    return esc(expr)
end
@inline function pivotCols!(A, K, i, M, piv)
    j = piv
    @inbounds while (A[i, piv] == 0)
        if ((piv += 1) > M)
            return true
        end
    end
    @inbounds if (j != piv)
        @vp for k in axes(A, 1)
            A[k, j], A[k, piv] = A[k, piv], A[k, j]
        end
        if K !== nothing
            @vp for k in axes(K, 1)
                K[k, j], K[k, piv] = K[k, piv], K[k, j]
            end
        end
    end
    return false
end
@inline function zeroSupDiagonal!(A, b, rr, c)
    M, N = size(A)
    @inbounds for j = c+1:N
        Aii = A[rr, c]
        Aij = A[rr, j]
        iszero(Aij) && continue
        if abs(Aii) == 1
            @vp for k = 1:M
                Ack = A[k, c]
                Ajk = A[k, j]
                A[k, c] = Aii * Ack
                A[k, j] = Aii * Ajk - Aij * Ack
            end
            if b !== nothing
                @vp for k in axes(b, 1)
                    bc, bj = b[k, c], b[k, j]
                    b[k, c] = Aii * bc
                    b[k, j] = Aii * bj  - Aij * bc
                end
            end
        else
            r, p, q = gcdx(Aii, Aij)
            Aiir = Base.sdiv_int(Aii, r)
            Aijr = Base.sdiv_int(Aij, r)
            @vp for k = 1:M
                Ack = A[k, c]
                Ajk = A[k, j]
                A[k, c] = p * Ack + q * Ajk
                A[k, j] = Aiir * Ajk - Aijr * Ack
            end
            if b !== nothing
                @vp for k in axes(b, 1)
                    bc, bj = b[k, c], b[k, j]
                    b[k, c] = p * bc + q * bj
                    b[k, j] = bj * Aiir - bc * Aijr
                end
            end
        end
    end
end
@inline function zeroSubDiagonal!(A, b, rr, c)
    N = size(A, 1)
    @inbounds for j = 1:c-1
        Aic = A[rr, c]
        Aij = A[rr, j]
        iszero(Aij) && continue
        g = gcd(Aic, Aij)
        Aicr = Base.sdiv_int(Aic, g)
        Aijr = Base.sdiv_int(Aij, g)
        @vp for k = 1:N
            Ack = A[k, c] * Aijr
            Ajk = A[k, j] * Aicr
            A[k, j] = Ajk - Ack
        end
        if b !== nothing
            @vp for k in axes(b, 1)
                Ack = b[k, c] * Aijr
                Ajk = b[k, j] * Aicr
                b[k, j] = Ajk - Ack
            end
        end
    end
end
@inline function simplify!(E, q = nothing)
    M, N = size(E)
    (N == 0) && return
    dec = 0
    for m = 1:M
        if (m - dec > N)
            break
        end
        if (pivotCols!(E, q, m, N, m - dec))
            dec += 1
            continue
        end
        zeroSupDiagonal!(E, q, m, m - dec)
        zeroSubDiagonal!(E, q, m, m - dec)
    end
end
nullspace(A) = nullspace!(copy(A))
@inline function identity!(C)
    @inbounds @vp for i in eachindex(C)
        C[i] = 0
    end
    @inbounds @vp for i = 1:size(C, 1)
        C[i, i] = 1
    end
end
@inline function _allZero(x)
    @inbounds for i in eachindex(x)
        x[i] == 0 || return false
    end
    return true
end
function nullspace!(B, C = Matrix{Int}(undef, size(B, 2), size(B, 2)))
    M = size(B, 2)
    identity!(C)
    simplify!(B, C)
    Mnew = M
    while (Mnew > 0) && _allZero(@inbounds(view(B, :, Mnew)))
        Mnew -= 1
    end
    @inbounds(view(C, :, Mnew+1:lastindex(C, 1)))
end


function smith!(A)
    H, U = hnf!(A)
    D, L = hnf!(H')
    L', D, U
end
smith(A) = smith!(copy(A))

end
