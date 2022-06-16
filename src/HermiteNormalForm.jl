module HermiteNormalForm

using LinearAlgebra
export hnf!, hnf, ishnf

"""
    ishnf(A)

Check if matrix `A` is in Hermite norm form.
"""
function ishnf(A)
    LinearAlgebra.istril(A) || return false

    m, n = size(A)
    k = 1
    while k <= min(m, n)
        pivot = find_pivot(A, k)
        if iszero(pivot)
            ki = findfirst(!iszero, @view A[k+1:end, k])
            if ki === nothing
                k = k + 1
                continue
            end
            ki += k
        else
            ki = k
        end

        p = A[ki, k]
        for j in 1:k-1
            0 <= A[ki, j] < p || return false
        end
        k = ki + 1
    end
    return true
end

"""
    H, U, rank = hnf(A)

`U` is unimodular. `H` is the Hermite normal form of `A`. `rank` is the rank of
`A`.

    H, rank = hnf(A, Val(false))
"""
hnf(A, with_transformation = Val(true)) = hnf!(copy(A), with_transformation)

function nullspace(A)
    n = size(A, 2)
    U = Matrix{eltype(A)}(I, n, n)
    H, U, rank = _hnf_like!(copy(A), U, Val(true))
    @view U[:, rank+1:end]
end

Base.@propagate_inbounds function find_pivot(A, k, zz = zero(eltype(A)))
    n = size(A, 2)
    pivot = k
    if n >= k && A[k, k] == zz
        for j in k+1:n
            if A[k, j] != zz
                pivot = j
            end
        end
        pivot == k && (pivot = zero(pivot)) # rank deficient
    end
    return pivot
end

function hnf!(A, ::Val{with_transformation}) where with_transformation
    n = size(A, 2)
    _hnf_like!(A, with_transformation ? Matrix{eltype(A)}(I, n, n) : nothing)
end
function _hnf_like!(A, U = Matrix{T}(I, n, n), ::Val{diagonalize} = Val(false)) where {diagonalize}
    Base.require_one_based_indexing(A)
    m, n = size(A)
    minmn = min(m, n)
    T = eltype(A)
    zz = zero(T)
    rank = 0
    @inbounds for k in 1:m
        pivot = find_pivot(A, k, zz)
        has_subdiag = n >= k
        if iszero(pivot) # rank deficient
            ki = findfirst(!iszero, @view A[k+1:end, k])
            ki === nothing && continue
            ki += k
        else
            ki = k
            if pivot != k
                for i in 1:m
                    A[i, k], A[i, pivot] = A[i, pivot], A[i, k]
                end
                U === nothing || for i in 1:n
                    U[i, k], U[i, pivot] = U[i, pivot], U[i, k]
                end
            end
        end
        rank += k <= minmn

        # [A11   0 ] k-1 rows
        # [A21  A22]
        # Zero out A[k, k+1:n] === A22[1, 2:end] by multiplying
        # [p  -A[k, j]/d]
        # [q   A[k, k]/d]
        for j in k+1:n
            Akk, Akj = A[ki, k], A[ki, j]
            d, p, q = gcdx(Akk, Akj)
            Akkd, Akjd = div(Akk, d), div(Akj, d)
            for i in 1:m
                Aik, Aij = A[i, k], A[i, j]
                A[i, k] = Aik * p + Aij * q
                A[i, j] = -Aik * Akjd + Aij * Akkd
            end
            U === nothing || for i in 1:n
                Uik, Uij = U[i, k], U[i, j]
                U[i, k] = Uik * p + Uij * q
                U[i, j] = -Uik * Akjd + Uij * Akkd
            end
        end

        has_subdiag || continue

        if diagonalize
            # Zero out A[k, 1:k-1] === A21[1, 1:k-1] by doing
            # A[:, j] = A[:, [k j]] * [-A[k, j], A[k, k]]
            #
            # Note that we then have:
            #   A[k, j] = A[k, [k j]] * [-A[k, j], A[k, k]]
            # = A[k, k] * -A[k, j] + A[k, j] * A[k, k] = 0
            for j in 1:k-1
                Akk, Akj = A[ki, k], A[ki, j]
                d = gcd(Akk, Akj)
                Akkd, Akjd = div(Akk, d), div(Akj, d)
                for i in 1:m
                    Aik, Aij = A[i, k], A[i, j]
                    A[i, j] = -Aik * Akjd + Aij * Akkd
                end
                U === nothing || for i in 1:n
                    Uik, Uij = U[i, k], U[i, j]
                    U[i, j] = -Uik * Akjd + Uij * Akkd
                end
            end
        else
            # Ensure the positivity of A[k, k]
            if A[ki, k] < zz
                @. A[:, k] = -A[:, k]
                U === nothing || @. U[:, k] = -U[:, k]
            end
            # Minimize A[k, 1:k-1] === A21[1, :]
            for j in 1:k-1
                mul = fld(A[ki, j], A[ki, k])
                @. A[:, j] -= mul * A[:, k]
                U === nothing || @. U[:, j] -= mul * U[:, k]
            end
        end
    end
    A, U, rank
end

end
