# BN254 Fp12 implementation as quadratic extension of Fp6.
#
# Realises Fp12 = Fp6[w]/(w² - v).

"""
    Fp12Element

Element of dodecic extension `Fp12 = Fp6[w] / (w^2 - v)`.
Represented as `c0 + c1*w` with `c0, c1 in Fp6`.
"""
struct Fp12Element
    c0::Fp6Element
    c1::Fp6Element
end

# Convenient constructors
Fp12Element(c0::Fp6Element) = Fp12Element(c0, zero(Fp6Element))
Fp12Element(c0) = Fp12Element(Fp6Element(c0), zero(Fp6Element))

const FQ12_NONRESIDUE = Fp6Element(zero(Fp2Element), one(Fp2Element), zero(Fp2Element))

"""
    mul_fp6_by_nonresidue(a::Fp6Element)

Multiply `a` by the quadratic-extension nonresidue `v = (0, 1, 0)`.
"""
mul_fp6_by_nonresidue(a::Fp6Element) = Fp6Element(mul_fp2_by_nonresidue(a.c2), a.c0, a.c1)

# Access components
function Base.getindex(a::Fp12Element, i::Int)
    if i == 1
        return a.c0
    elseif i == 2
        return a.c1
    end
    throw(BoundsError(a, i))
end

# Basic operations
Base.zero(::Type{Fp12Element}) = Fp12Element(zero(Fp6Element), zero(Fp6Element))
Base.one(::Type{Fp12Element}) = Fp12Element(one(Fp6Element), zero(Fp6Element))
Base.zero(::Fp12Element) = zero(Fp12Element)
Base.one(::Fp12Element) = one(Fp12Element)

Base.iszero(a::Fp12Element) = iszero(a.c0) && iszero(a.c1)
Base.isone(a::Fp12Element) = isone(a.c0) && iszero(a.c1)

# Equality
Base.:(==)(a::Fp12Element, b::Fp12Element) = a.c0 == b.c0 && a.c1 == b.c1
Base.isequal(a::Fp12Element, b::Fp12Element) = isequal(a.c0, b.c0) && isequal(a.c1, b.c1)

# Addition and subtraction
Base.:+(a::Fp12Element, b::Fp12Element) = Fp12Element(a.c0 + b.c0, a.c1 + b.c1)
Base.:-(a::Fp12Element, b::Fp12Element) = Fp12Element(a.c0 - b.c0, a.c1 - b.c1)
Base.:-(a::Fp12Element) = Fp12Element(-a.c0, -a.c1)

"""
    *(a::Fp12Element, b::Fp12Element)

Multiply two `Fp12Element`s using the relation `w² = v`.
For `(a₀ + a₁w)(b₀ + b₁w)` the product satisfies
`c₀ = a₀b₀ + a₁b₁v` and `c₁ = a₀b₁ + a₁b₀`.
"""
function Base.:*(a::Fp12Element, b::Fp12Element)
    aa = a.c0 * b.c0
    bb = a.c1 * b.c1
    c0 = aa + mul_fp6_by_nonresidue(bb)
    c1 = (a.c0 + a.c1) * (b.c0 + b.c1) - aa - bb
    return Fp12Element(c0, c1)
end

"""
    square(a::Fp12Element)

Square an `Fp12Element` using the optimized quadratic-extension formula.
"""
function square(a::Fp12Element)
    c0_sq = square(a.c0)
    c1_sq = square(a.c1)
    c0 = c0_sq + mul_fp6_by_nonresidue(c1_sq)
    c1 = (a.c0 * a.c1) * 2
    return Fp12Element(c0, c1)
end

Base.:^(a::Fp12Element, ::Val{2}) = square(a)

"""
    conjugate(a::Fp12Element)

Return the conjugate `(a₀ - a₁w)` of `(a₀ + a₁w)` in Fp12.
"""
conjugate(a::Fp12Element) = Fp12Element(a.c0, -a.c1)

"""
    inv(a::Fp12Element)

Compute the multiplicative inverse `(a₀ - a₁w) / (a₀² - a₁²v)` in Fp12.
"""
function Base.inv(a::Fp12Element)
    if iszero(a)
        throw(DivideError())
    end

    norm = square(a.c0) - mul_fp6_by_nonresidue(square(a.c1))
    norm_inv = inv(norm)

    conj = conjugate(a)
    return Fp12Element(conj.c0 * norm_inv, conj.c1 * norm_inv)
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
Base.:*(k::Integer, a::Fp12Element) = Fp12Element(k * a.c0, k * a.c1)
Base.:*(a::Fp12Element, k::Integer) = k * a

# Display
function Base.show(io::IO, a::Fp12Element)
    print(io, "Fp12([", a.c0, ", ", a.c1, "])")
end

# Type alias for GT elements
const GTElement = Fp12Element

# Export types and functions
export Fp12Element, FQ12_NONRESIDUE, GTElement, mul_fp6_by_nonresidue, square, conjugate, frobenius
