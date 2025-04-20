# PolySet.jl

[![Build Status](https://github.com/Theozeud/PolySet.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Theozeud/PolySet.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/Theozeud/PolySet.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/Theozeud/PolySet.jl)
![License](https://img.shields.io/badge/license-MIT-blue.svg)


**PolySet.jl** is a Julia package for efficient and vectorized manipulation of univariate polynomial sets.

## ✨ Motivation

When working with a large number of univariate polynomials, it is often desirable to avoid iterative, coefficient-by-coefficient manipulations. `PolySet.jl` is designed to provide a fast and structured way to store and operate on collections of polynomials using array-based representations.

The core idea is to store multiple polynomials in a single 2D matrix, where:
- each **row** corresponds to a polynomial,
- each **column** holds the coefficient for a given degree (in increasing order),
- and operations such as evaluation, differentiation, integration, or addition are performed in a fully vectorized way.

This structure is particularly well-suited for numerical applications, symbolic prototyping, and situations where many polynomials are manipulated simultaneously.

## ⚙️ Features

- Compact storage of univariate polynomials
- Fast vectorized evaluation using Horner's method
- Vectorized differentiation
- Support for polynomial basis generation (e.g. monomials)
- Interoperable with standard Julia arrays
- Compatible with GPU acceleration (planned)
