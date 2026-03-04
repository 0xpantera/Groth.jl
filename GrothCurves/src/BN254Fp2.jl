# BN254 curve field implementations.
#
# Defines the BN254 base and quadratic extension fields used throughout the
# pairing engine.

# using GrothAlgebra
# using StaticArrays

# BN254Fq is now provided by GrothAlgebra
# Just re-export what we need
const BN254_PRIME = prime(BN254Fq)

"""
    Fp2Element

Element of the quadratic extension field Fp2 = Fp[u]/(u² + 1).
Represented as c0 + c1*u where c0, c1 ∈ Fp.
"""
struct Fp2Element
    coeffs::SVector{2,BN254Fq}

    function Fp2Element(c0::BN254Fq, c1::BN254Fq)
        new(SVector(c0, c1))
    end
end

# Convenient constructors
Fp2Element(c0::Integer, c1::Integer) = Fp2Element(bn254_fq(c0), bn254_fq(c1))
Fp2Element(c0) = Fp2Element(bn254_fq(c0), zero(BN254Fq))

# Access components
Base.getindex(a::Fp2Element, i::Int) = a.coeffs[i]
real(a::Fp2Element) = a.coeffs[1]
imag(a::Fp2Element) = a.coeffs[2]

# Basic operations
Base.zero(::Type{Fp2Element}) = Fp2Element(zero(BN254Fq), zero(BN254Fq))
Base.one(::Type{Fp2Element}) = Fp2Element(one(BN254Fq), zero(BN254Fq))
Base.zero(::Fp2Element) = zero(Fp2Element)
Base.one(::Fp2Element) = one(Fp2Element)

Base.iszero(a::Fp2Element) = iszero(a[1]) && iszero(a[2])
Base.isone(a::Fp2Element) = isone(a[1]) && iszero(a[2])

# Equality
Base.:(==)(a::Fp2Element, b::Fp2Element) = a.coeffs == b.coeffs
Base.isequal(a::Fp2Element, b::Fp2Element) = isequal(a.coeffs, b.coeffs)

# Arithmetic operations
Base.:+(a::Fp2Element, b::Fp2Element) = Fp2Element(a[1] + b[1], a[2] + b[2])
Base.:-(a::Fp2Element, b::Fp2Element) = Fp2Element(a[1] - b[1], a[2] - b[2])
Base.:-(a::Fp2Element) = Fp2Element(-a[1], -a[2])

"""
    *(a::Fp2Element, b::Fp2Element)

Multiply two `Fp2Element`s using the relation `u² = -1`.
The product `(a₀ + a₁u)(b₀ + b₁u)` equals `(a₀b₀ - a₁b₁) + (a₀b₁ + a₁b₀)u`.
"""
function Base.:*(a::Fp2Element, b::Fp2Element)
    # Real part: a0*b0 - a1*b1
    c0 = a[1] * b[1] - a[2] * b[2]
    # Imaginary part: a0*b1 + a1*b0
    c1 = a[1] * b[2] + a[2] * b[1]
    return Fp2Element(c0, c1)
end

# Scalar multiplication
Base.:*(a::Fp2Element, k::Integer) = Fp2Element(a[1] * k, a[2] * k)
Base.:*(k::Integer, a::Fp2Element) = a * k
Base.:*(a::Fp2Element, k::BN254Fq) = Fp2Element(a[1] * k, a[2] * k)
Base.:*(k::BN254Fq, a::Fp2Element) = a * k

"""
    conjugate(a::Fp2Element)

Return the conjugate `a - b*u` of `a + b*u`.
"""
conjugate(a::Fp2Element) = Fp2Element(a[1], -a[2])

"""
    norm(a::Fp2Element)

Compute the norm `a² + b²` of `a + b*u`.
"""
norm(a::Fp2Element) = a[1]^2 + a[2]^2

"""
    inv(a::Fp2Element)

Compute the multiplicative inverse `(a - b*u) / (a² + b²)`.
"""
function Base.inv(a::Fp2Element)
    if iszero(a)
        throw(DivideError())
    end
    n = norm(a)
    n_inv = inv(n)
    conj = conjugate(a)
    return Fp2Element(conj[1] * n_inv, conj[2] * n_inv)
end

Base.:/(a::Fp2Element, b::Fp2Element) = a * inv(b)

"""
    ^(a::Fp2Element, n::Integer)

Raise an `Fp2Element` to the integer power `n` using binary exponentiation.
"""
function Base.:^(a::Fp2Element, n::Integer)
    if n == 0
        return one(Fp2Element)
    elseif n < 0
        return inv(a)^(-n)
    end

    result = one(Fp2Element)
    base = a
    exp = n

    while exp > 0
        if exp & 1 == 1
            result = result * base
        end
        base = base * base
        exp >>= 1
    end

    return result
end

"""
    frobenius(a::Fp2Element)

Apply the Frobenius endomorphism `(a + b*u) ↦ (a - b*u)` (since `u^p = -u`).
"""
frobenius(a::Fp2Element) = conjugate(a)

# Display
function Base.show(io::IO, a::Fp2Element)
    if iszero(a[2])
        print(io, "Fp2(", a[1].value, ")")
    else
        print(io, "Fp2(", a[1].value, " + ", a[2].value, "*u)")
    end
end

# Export types and functions
export BN254Fq, bn254_fq, BN254_PRIME
export Fp2Element, conjugate, norm, frobenius, real, imag
