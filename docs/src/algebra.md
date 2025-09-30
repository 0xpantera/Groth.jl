# Groth Algebra

```@meta
CurrentModule = GrothAlgebra
DocTestSetup = :(using GrothAlgebra)
```

```@contents
Depth = 2
```

## Overview

- Prime-field arithmetic for BN254 and secp256k1.
- Generic polynomial utilities and group abstractions consumed by the curve and proof layers.

## Key Modules

- `src/FiniteFields.jl` – BigInt-backed fields with normalisation, `invmod`, and `powermod` helpers.
- `src/Polynomial.jl` – Degree/leading-coefficient utilities, Horner evaluation, interpolation, derivative, and FFT scaffolding.
- `src/Group.jl` – Generic group element interface with scalar multiplication, w-NAF helpers, and MSM support.

## Implementation Notes

- Added `interpolate_prefix_points` so subset domains recover coefficients before padding for the coset FFT path.
- FFT helpers now gate the dense fallback behind parity assertions.

## Follow-ups

- Evaluate FFT twiddle caching and mixed-radix support once QAP domain alignment lands.

## Finite Fields

```@docs
FiniteFieldElement
BN254Field
BN254ScalarField
Secp256k1Field
```

## Polynomials

```@docs
GrothAlgebra.Polynomial
GrothAlgebra.degree
leading_coefficient
evaluate
roots_of_unity
fft_polynomial_multiply
```

### Example

```@example
using GrothAlgebra

F = BN254ScalarField
p = Polynomial([F(3), F(2), F(1)])  # 1x^2 + 2x + 3
x = F(5)
evaluate(p, x)
```

## Group Interfaces and MSM

```@docs
GrothAlgebra.GroupElem
GrothAlgebra.scalar_mul
GrothAlgebra.scalar_mul_wnaf
GrothAlgebra.multi_scalar_mul
GrothAlgebra.FixedBaseTable
GrothAlgebra.build_fixed_table
GrothAlgebra.mul_fixed
```

### Example

```@example
using GrothAlgebra

struct DummyCurve <: AbstractCurve end

struct DummyPoint <: GroupElem{DummyCurve}
    value::Int
end

Base.zero(::Type{DummyPoint}) = DummyPoint(0)
Base.:+(a::DummyPoint, b::DummyPoint) = DummyPoint(a.value + b.value)
Base.:-(a::DummyPoint) = DummyPoint(-a.value)
DummyPoint(3)

scalar_mul(DummyPoint(3), 5)
```
