module HermiteNormalForm

using LinearAlgebra
export hnf!, hnf, ishnf

"""
    ishnf(A)

Check if matrix `A` is in Hermite norm form.
"""
function ishnf(A)
    m, n = size(A)
    for j in 1:n
        for i in 1:min(m, j-1)
            A[i, j] == 0 || return false
        end
        if m >= j
            A[j, j] > 0 || return false
        end
    end
    for i in 1:min(m, n)
        p = A[i, i]
        for j in 1:min(n, i-1)
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
    for k in 1:m
        # Pivot: A[k, k] should not be 0
        if n >= k && A[k, k] == zz
            pivot = k
            for j in k+1:n
                if A[k, j] != zz
                    pivot = j
                end
            end
            pivot == k && return nothing # rank deficient
            for i in 1:m
                A[i, k], A[i, pivot] = A[i, pivot], A[i, k]
            end
            for i in 1:n
                U[i, k], U[i, pivot] = U[i, pivot], U[i, k]
            end
        end
        # [A11   0 ] k-1 rows
        # [A21  A22]
        # Zero out A[k, k+1:n] === A22[1, 2:end] by multiplying
        # [p  -A[k, j]/d]
        # [q   A[k, k]/d]
        for j in k+1:n
            Akk, Akj = A[k, k], A[k, j]
            d, p, q = gcdx(Akk, Akj)
            Akkd, Akjd = div(Akk, d), div(Akj, d)
            for i in 1:m
                Aik, Aij = A[i, k], A[i, j]
                A[i, k] = Aik * p + Aij * q
                A[i, j] = -Aik * Akjd + Aij * Akkd
            end
            for i in 1:n
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
        for j in 1:k-1
            mul = fld(A[k, j], A[k, k])
            @. A[:, j] -= mul * A[:, k]
            @. U[:, j] -= mul * U[:, k]
        end
    end
    A, U
end

end
