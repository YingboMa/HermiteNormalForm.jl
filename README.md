# HermiteNormalForm

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://YingboMa.github.io/HermiteNormalForm.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://YingboMa.github.io/HermiteNormalForm.jl/dev)
[![Build Status](https://github.com/YingboMa/HermiteNormalForm.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/YingboMa/HermiteNormalForm.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/YingboMa/HermiteNormalForm.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/YingboMa/HermiteNormalForm.jl)

This package implements the Hermite norm form. It can be used to compute the
unimodular extension of an integer matrix. In the following example, we compute
the unimodular matrix `Uinv` which shares rows with the given rectangular matrix
`A`.

Reference: https://math.stackexchange.com/a/4397337/662983

```julia
julia> using HermiteNormalForm, LinearAlgebra

julia> A = [
        -2  -1   0   3   0
        -3   0  -2  -1   1
         3  -1   1  -3  -1
       ];

julia> H, U = hnf(A); @assert all(isone, diag(H))

julia> Uinv = Int.(inv(U.//1))
5×5 Matrix{Int64}:
 -2  -1   0   3   0
 -3   0  -2  -1   1
  3  -1   1  -3  -1
 -7   0  -4   0   0
  3   0   2   1   0

julia> Int(det(Uinv.//1))
1
```
