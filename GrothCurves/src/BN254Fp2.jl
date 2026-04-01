# BN254 curve field implementations.
#
# Defines the BN254 base and quadratic extension fields used throughout the
# pairing engine.

# BN254Fq is now provided by GrothAlgebra
# Just re-export what we need
const BN254_PRIME = prime(BN254Fq)

"""
    Fp2Element

Element of the quadratic extension field `Fp2 = Fp[u] / (u^2 + 1)`.
Represented as `c0 + c1*u` with `c0, c1 in Fp`.
"""
struct Fp2Element
    c0::BN254Fq
    c1::BN254Fq
end

# Convenient constructors
Fp2Element(c0::Integer, c1::Integer) = Fp2Element(bn254_fq(c0), bn254_fq(c1))
Fp2Element(c0::BN254Fq) = Fp2Element(c0, zero(BN254Fq))
Fp2Element(c0::Integer) = Fp2Element(bn254_fq(c0), zero(BN254Fq))

# Access components
function Base.getindex(a::Fp2Element, i::Int)
    if i == 1
        return a.c0
    elseif i == 2
        return a.c1
    end
    throw(BoundsError(a, i))
end
real(a::Fp2Element) = a.c0
imag(a::Fp2Element) = a.c1

# Basic operations
Base.zero(::Type{Fp2Element}) = Fp2Element(zero(BN254Fq), zero(BN254Fq))
Base.one(::Type{Fp2Element}) = Fp2Element(one(BN254Fq), zero(BN254Fq))
Base.zero(::Fp2Element) = zero(Fp2Element)
Base.one(::Fp2Element) = one(Fp2Element)

Base.iszero(a::Fp2Element) = iszero(a.c0) && iszero(a.c1)
Base.isone(a::Fp2Element) = isone(a.c0) && iszero(a.c1)

# Equality
Base.:(==)(a::Fp2Element, b::Fp2Element) = a.c0 == b.c0 && a.c1 == b.c1
Base.isequal(a::Fp2Element, b::Fp2Element) = isequal(a.c0, b.c0) && isequal(a.c1, b.c1)

# Arithmetic operations
Base.:+(a::Fp2Element, b::Fp2Element) = Fp2Element(a.c0 + b.c0, a.c1 + b.c1)
Base.:-(a::Fp2Element, b::Fp2Element) = Fp2Element(a.c0 - b.c0, a.c1 - b.c1)
Base.:-(a::Fp2Element) = Fp2Element(-a.c0, -a.c1)

function square(a::Fp2Element)
    ab = a.c0 * a.c1
    c0 = (a.c0 + a.c1) * (a.c0 - a.c1)
    c1 = ab + ab
    return Fp2Element(c0, c1)
end

"""
    *(a::Fp2Element, b::Fp2Element)

Multiply two `Fp2Element`s using the relation `u² = -1`.
The product `(a₀ + a₁u)(b₀ + b₁u)` equals `(a₀b₀ - a₁b₁) + (a₀b₁ + a₁b₀)u`.
"""
function Base.:*(a::Fp2Element, b::Fp2Element)
    aa = a.c0 * b.c0
    bb = a.c1 * b.c1
    c0 = aa - bb
    c1 = (a.c0 + a.c1) * (b.c0 + b.c1) - aa - bb
    return Fp2Element(c0, c1)
end

# Scalar multiplication
Base.:*(a::Fp2Element, k::Integer) = Fp2Element(a.c0 * k, a.c1 * k)
Base.:*(k::Integer, a::Fp2Element) = a * k
Base.:*(a::Fp2Element, k::BN254Fq) = Fp2Element(a.c0 * k, a.c1 * k)
Base.:*(k::BN254Fq, a::Fp2Element) = a * k

"""
    conjugate(a::Fp2Element)

Return the conjugate `a - b*u` of `a + b*u`.
"""
conjugate(a::Fp2Element) = Fp2Element(a.c0, -a.c1)

"""
    norm(a::Fp2Element)

Compute the norm `a² + b²` of `a + b*u`.
"""
norm(a::Fp2Element) = a.c0^2 + a.c1^2

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
    return Fp2Element(conj.c0 * n_inv, conj.c1 * n_inv)
end

Base.:/(a::Fp2Element, b::Fp2Element) = a * inv(b)

"""
    ^(a::Fp2Element, n::Integer)

Raise an `Fp2Element` to the integer power `n` using binary exponentiation.
"""
function Base.:^(a::Fp2Element, n::Integer)
    if n == 0
        return one(Fp2Element)
    elseif n == 2
        return square(a)
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
        base = square(base)
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
    if iszero(a.c1)
        print(io, "Fp2(", convert(BigInt, a.c0), ")")
    else
        print(io, "Fp2(", convert(BigInt, a.c0), " + ", convert(BigInt, a.c1), "*u)")
    end
end

# Export types and functions
export BN254Fq, bn254_fq, BN254_PRIME
export Fp2Element, conjugate, norm, frobenius, real, imag, square
