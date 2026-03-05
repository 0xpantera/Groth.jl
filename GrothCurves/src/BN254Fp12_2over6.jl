# BN254 Fp12 implementation as quadratic extension of Fp6.
#
# Realises Fp12 = Fp6[w]/(w² - v).

# using StaticArrays

"""
    Fp12Element

Element of dodecic extension
``\\mathbb{F}_{p^{12}} = \\mathbb{F}_{p^6}[w]/(w^2 - v)``.
Represented as ``c_0 + c_1 w`` where ``c_0, c_1 \\in \\mathbb{F}_{p^6}``.
"""
struct Fp12Element
    coeffs::SVector{2,Fp6Element}

    function Fp12Element(c0::Fp6Element, c1::Fp6Element)
        new(SVector(c0, c1))
    end
end

# Convenient constructors
Fp12Element(c0) = Fp12Element(Fp6Element(c0), zero(Fp6Element))

# Access components
Base.getindex(a::Fp12Element, i::Int) = a.coeffs[i]

# Basic operations
Base.zero(::Type{Fp12Element}) = Fp12Element(zero(Fp6Element), zero(Fp6Element))
Base.one(::Type{Fp12Element}) = Fp12Element(one(Fp6Element), zero(Fp6Element))
Base.zero(::Fp12Element) = zero(Fp12Element)
Base.one(::Fp12Element) = one(Fp12Element)

Base.iszero(a::Fp12Element) = iszero(a[1]) && iszero(a[2])
Base.isone(a::Fp12Element) = isone(a[1]) && iszero(a[2])

# Equality
Base.:(==)(a::Fp12Element, b::Fp12Element) = a.coeffs == b.coeffs

# Addition and subtraction
Base.:+(a::Fp12Element, b::Fp12Element) = Fp12Element(a[1] + b[1], a[2] + b[2])
Base.:-(a::Fp12Element, b::Fp12Element) = Fp12Element(a[1] - b[1], a[2] - b[2])
Base.:-(a::Fp12Element) = Fp12Element(-a[1], -a[2])

"""
    *(a::Fp12Element, b::Fp12Element)

Multiply two `Fp12Element`s using the relation `w² = v`.
For `(a₀ + a₁w)(b₀ + b₁w)` the product satisfies
`c₀ = a₀b₀ + a₁b₁v` and `c₁ = a₀b₁ + a₁b₀`.
"""
function Base.:*(a::Fp12Element, b::Fp12Element)
    # v is represented as (0, 1, 0) in Fp6
    v = Fp6Element(zero(Fp2Element), one(Fp2Element), zero(Fp2Element))

    c0 = a[1] * b[1] + a[2] * b[2] * v
    c1 = a[1] * b[2] + a[2] * b[1]

    return Fp12Element(c0, c1)
end

"""
    square(a::Fp12Element)

Square an `Fp12Element` using the optimized quadratic-extension formula.
"""
function square(a::Fp12Element)
    v = Fp6Element(zero(Fp2Element), one(Fp2Element), zero(Fp2Element))

    c0 = square(a[1]) + square(a[2]) * v
    c1 = Fp6Element(2) * a[1] * a[2]

    return Fp12Element(c0, c1)
end

Base.:^(a::Fp12Element, ::Val{2}) = square(a)

"""
    conjugate(a::Fp12Element)

Return the conjugate `(a₀ - a₁w)` of `(a₀ + a₁w)` in Fp12.
"""
conjugate(a::Fp12Element) = Fp12Element(a[1], -a[2])

"""
    inv(a::Fp12Element)

Compute the multiplicative inverse `(a₀ - a₁w) / (a₀² - a₁²v)` in Fp12.
"""
function Base.inv(a::Fp12Element)
    if iszero(a)
        throw(DivideError())
    end

    v = Fp6Element(zero(Fp2Element), one(Fp2Element), zero(Fp2Element))

    # Norm = a0² - a1²*v
    norm = square(a[1]) - square(a[2]) * v
    norm_inv = inv(norm)

    # (a0 - a1*w) / norm
    conj = conjugate(a)
    return Fp12Element(conj[1] * norm_inv, conj[2] * norm_inv)
end

Base.:/(a::Fp12Element, b::Fp12Element) = a * inv(b)

"""
    ^(a::Fp12Element, n::Integer)

Raise an `Fp12Element` to the integer power `n` with binary exponentiation.
"""
function Base.:^(a::Fp12Element, n::Integer)
    if n == 0
        return one(Fp12Element)
    elseif n == 2
        return square(a)
    elseif n < 0
        return inv(a)^(-n)
    end

    result = one(Fp12Element)
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
    frobenius(a::Fp12Element)

Apply the Frobenius endomorphism in Fp12.
For BN254, this map has exploitable structure.
"""
function frobenius(a::Fp12Element)
    # For now, implement basic version
    # Can optimize with precomputed constants later
    return a^BN254_PRIME
end

"""
    frobenius(a::Fp12Element, k::Integer)

Apply the `k`-fold Frobenius endomorphism.
"""
function frobenius(a::Fp12Element, k::Integer)
    result = a
    for _ in 1:k
        result = frobenius(result)
    end
    return result
end

# Scalar multiplication for convenience
Base.:*(k::Integer, a::Fp12Element) = Fp12Element(k * a[1], k * a[2])
Base.:*(a::Fp12Element, k::Integer) = k * a

# Display
function Base.show(io::IO, a::Fp12Element)
    print(io, "Fp12([", a[1], ", ", a[2], "])")
end

# Type alias for GT elements
const GTElement = Fp12Element

# Export types and functions
export Fp12Element, GTElement, square, conjugate, frobenius
