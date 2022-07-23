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
        for j = 1:k-1
            0 <= A[ki, j] < p || return false
        end
        k = ki + 1
    end
    return true
end

struct ElementaryUnimodularReduction{I<:Integer}
    j::I
    k::I
end

# zero out [j, k] using column j with pivot (j, j)
Base.:*(A::AbstractMatrix, e::ElementaryUnimodularReduction) = rmul!(copy(A), e)
function LinearAlgebra.rmul!(A::AbstractMatrix, e::ElementaryUnimodularReduction, B = nothing)
    k = e.j
    j = e.k
    Akj = A[k, j]
    Akk = A[k, k]
    d, p, q = gcdx(Akk, Akj)
    c, s = div(Akj, d), div(Akk, d)
    for i = axes(A, 1)
        Aik, Aij = A[i, k], A[i, j]
        A[i, k] = Aik * p + Aij * q
        A[i, j] = -Aik * c + Aij * s
    end
    B === nothing || for i = axes(B, 1)
        Bik, Bij = B[i, k], B[i, j]
        B[i, k] = Bik * p + Bij * q
        B[i, j] = -Bik * c + Bij * s
    end
    A
end

# zero out [j, k] using row k with pivot (k, k)
Base.:*(e::ElementaryUnimodularReduction, A::AbstractMatrix) = lmul!(e, copy(A))
function LinearAlgebra.lmul!(e::ElementaryUnimodularReduction, A::AbstractMatrix, B = nothing)
    j = e.j
    k = e.k
    Ajk = A[j, k]
    Akk = A[k, k]
    d, p, q = gcdx(Akk, Ajk)
    c, s = div(Ajk, d), div(Akk, d)
    for i = axes(A, 2)
        Aki, Aji = A[k, i], A[j, i]
        A[k, i] = Aki * p + Aji * q
        A[j, i] = -Aki * c + Aji * s
    end
    B === nothing || for i = axes(B, 2)
        Bki, Bji = B[k, i], B[j, i]
        B[k, i] = Bki * p + Bji * q
        B[j, i] = -Bki * c + Bji * s
    end
    A
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
    @view U[:, findall(iszero, eachcol(H))]
end

function simplify_linearsys(A)
    A = copy(A)
    _hnf_like!(A, Matrix{eltype(A)}(I, size(A, 2), size(A, 2)), Val(true))
    A
end

Base.@propagate_inbounds function find_pivot(A, k, zz = zero(eltype(A)))
    n = size(A, 2)
    pivot = k
    if n >= k && A[k, k] == zz
        for j = k+1:n
            if A[k, j] != zz
                pivot = j
            end
        end
        pivot == k && (pivot = zero(pivot)) # rank deficient
    end
    return pivot
end

function hnf!(A, ::Val{with_transformation}) where {with_transformation}
    n = size(A, 2)
    _hnf_like!(A, with_transformation ? Matrix{eltype(A)}(I, n, n) : nothing)
end
function _hnf_like!(
    A,
    U = Matrix{eltype(A)}(I, size(A, 2), size(A, 2)),
    ::Val{diagonalize} = Val(false),
) where {diagonalize}
    Base.require_one_based_indexing(A)
    m, n = size(A)
    minmn = min(m, n)
    T = eltype(A)
    zz = zero(T)
    rank = 0
    @inbounds for k = 1:m
        pivot = find_pivot(A, k, zz)
        has_subdiag = n >= k
        if iszero(pivot) # rank deficient
            ki = findfirst(!iszero, @view A[k+1:end, k])
            ki === nothing && continue
            ki += k
        else
            ki = k
            if pivot != k
                for i = 1:m
                    A[i, k], A[i, pivot] = A[i, pivot], A[i, k]
                end
                U === nothing || for i = 1:n
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
        # E.g.
        # [A[:, k] A[:, j]] [p  -A[k, j]/d] = [p * A[:, k] + q * A[:, j]   -A[k, j]/d * A[:, k] + A[k, k]/d * A[:, j]]
        #                   [q   A[k, k]/d]
        #
        # [p                   q] [A[k, :]  = [p * A[k, :] + q * A[j, :]
        # [-A[j, k]/d  A[k, k]/d]  A[j, :]]  -A[j, k]/d * A[k, :] + A[k, k]/d * A[j, :]]
        # A x = 0
        # U A x = U 0 = 0
        # H x = 0
        for j = k+1:n
            Akk, Akj = A[ki, k], A[ki, j]
            d, p, q = gcdx(Akk, Akj)
            Akkd, Akjd = div(Akk, d), div(Akj, d)
            for i = 1:m
                Aik, Aij = A[i, k], A[i, j]
                A[i, k] = Aik * p + Aij * q
                A[i, j] = -Aik * Akjd + Aij * Akkd
            end
            U === nothing || for i = 1:n
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
            for j = 1:k-1
                Akk, Akj = A[ki, k], A[ki, j]
                d = gcd(Akk, Akj)
                Akkd, Akjd = div(Akk, d), div(Akj, d)
                for i = 1:m
                    Aik, Aij = A[i, k], A[i, j]
                    A[i, j] = -Aik * Akjd + Aij * Akkd
                end
                U === nothing || for i = 1:n
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
            for j = 1:k-1
                mul = fld(A[ki, j], A[ki, k])
                @. A[:, j] -= mul * A[:, k]
                U === nothing || @. U[:, j] -= mul * U[:, k]
            end
        end
    end
    A, U, rank
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
            for k = 1:M
                Ack = A[k, c]
                Ajk = A[k, j]
                A[k, c] = Aii * Ack
                A[k, j] = Aii * Ajk - Aij * Ack
            end
            if b !== nothing
                for k in axes(b, 1)
                    bc, bj = b[k, c], b[k, j]
                    b[k, c] = Aii * bc
                    b[k, j] = Aii * bj - Aij * bc
                end
            end
        else
            r, p, q = gcdx(Aii, Aij)
            #Aiir = Base.sdiv_int(Aii, r)
            #Aijr = Base.sdiv_int(Aij, r)
            Aiir = Base.div(Aii, r)
            Aijr = Base.div(Aij, r)
            for k = 1:M
                Ack = A[k, c]
                Ajk = A[k, j]
                A[k, c] = p * Ack + q * Ajk
                A[k, j] = Aiir * Ajk - Aijr * Ack
            end
            if b !== nothing
                for k in axes(b, 1)
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
        Aicr = Base.div(Aic, g)
        Aijr = Base.div(Aij, g)
        for k = 1:N
            Ack = A[k, c] * Aijr
            Ajk = A[k, j] * Aicr
            A[k, j] = Ajk - Ack
        end
        if b !== nothing
            for k in axes(b, 1)
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
    E
end
nullspace2(A) = nullspace2!(copy(A))
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
function nullspacecore!(B, C)
    M = size(B, 2)
    identity!(C)
    simplify!(B, C)
    Mnew = M
    while (Mnew > 0) && _allZero(@inbounds(view(B, :, Mnew)))
        Mnew -= 1
    end
    return Mnew
end

function nullspace2!(B, C = Matrix{Int}(undef, size(B, 2), size(B, 2)))
    M = nullspacecore!(B, C)
    @inbounds(view(C, :, M+1:lastindex(C, 1)))
end

orthogonalnullbasis(A) = orthogonalnullbasis!(copy(A))
function orthogonalnullbasis!(
    A,
    B = Matrix{eltype(A)}(undef, size(A, 2), size(A, 2)),
    C = similar(B),
)
    M = size(C, 1)
    B[1:size(A, 1), :] .= A
    simplify!(@view(B[1:size(A, 1), :])')
    Mnew = Morig = nullspacecore!(A, C)
    Mnew == M && return @view(B[Morig+1:end, :])
    B[Mnew+1, :] .= @view(C[:, Mnew+1])
    B2 = similar(B)
    i = 0
    while Mnew < M
        B2 .= B
        Mnew = nullspacecore!(@view(B2[1:Mnew+1, :]), C)
        Mnew == M && break
        B[Mnew+1, :] .= @view(C[:, Mnew+1])
        @assert ((i += 1) < 5)
    end
    @view(B[Morig+1:end, :])
end

proj!(out, u, v) = out .-= u'v // (u'u) .* u
orthogonalize(A) = orthogonalize!(copy(A))
function orthogonalize!(O::AbstractMatrix{T}) where {T<:Integer}
    buffer = similar(O, Rational{T}, size(O, 1))
    for i = 1:size(O, 2)
        copyto!(buffer, @view O[:, i])
        for j = 1:i-1
            @views proj!(buffer, O[:, j], O[:, i])
        end
        lm = mapreduce(denominator, lcm, buffer; init = one(T))
        @views @. O[:, i] = numerator(buffer * lm)
    end
    O
end
orthogonalnullbasis2(A) = orthogonalize!(nullspace2(A))



end
