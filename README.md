# PolySets.jl

[![Build Status](https://github.com/Theozeud/PolySets.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Theozeud/PolySets.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/Theozeud/PolySets.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/Theozeud/PolySet.jl)
![License](https://img.shields.io/badge/license-MIT-blue.svg)


**PolySets.jl** is a Julia package for efficient and vectorized manipulation of univariate polynomial sets.

## Motivation

When working with a large number of univariate polynomials, it is often desirable to avoid iterative, coefficient-by-coefficient manipulations. `PolySets.jl` is designed to provide a fast and structured way to store and operate on collections of polynomials using array-based representations.

The core idea is to store multiple polynomials in a single 2D matrix, where:
- each **row** corresponds to a polynomial,
- each **column** holds the coefficient for a given degree (in increasing order),
- and operations such as evaluation, differentiation, integration, or addition are performed in a fully vectorized way.

This structure is particularly well-suited for numerical applications where many polynomials are manipulated simultaneously.

## Features

- Compact storage of univariate polynomials
- Fast vectorized evaluation using Horner's method
- Vectorized operations (derivation, integration)
- Support for polynomial basis generation (e.g. monomials, Legendre, IntLegendre)
- Interoperable with standard Julia arrays
- Compatible with GPU acceleration (planned)


## How to use `PolySet`

Here's a quick example showing how to create a set of 1000 polynomials of degree 999 and evaluate them on 1000 points:

```julia
using YourPolyPackage  # Replace with the actual package name
using BenchmarkTools
using Plots

# Generate a PolySet of size 1000 × 1000 (1000 polynomials of degree 999)
n_polys = 1000
degree = 999
ps = PolySet(randn(n_polys, degree + 1))  # Random coefficients

# Generate 1000 evaluation points in [-1, 1]
x = range(-1.0, 1.0; length=n_polys)

# Evaluate all polynomials at all points (vectorized broadcasting)
@btime y = ps.(x)  # Returns a (1000, 1000) matrix

# Plot a few polynomials
plot(x, ps[1:5].(x)', legend=false, title="First 5 polynomials")
```

