# BN254 Fp6 implementation as cubic extension of Fp2.
#
# Realises Fp6 = Fp2[v]/(v³ - ξ) with ξ = 9 + u ∈ Fp2.

# using StaticArrays

"""
    Fp6Element

Element of sextic extension Fp6 = Fp2[v]/(v³ - ξ) where ξ = 9 + u.
Represented as c0 + c1*v + c2*v² where c0, c1, c2 ∈ Fp2.
"""
struct Fp6Element
    coeffs::SVector{3,Fp2Element}

    function Fp6Element(c0::Fp2Element, c1::Fp2Element, c2::Fp2Element)
        new(SVector(c0, c1, c2))
    end
end

# Convenient constructors
Fp6Element(c0, c1, c2) = Fp6Element(Fp2Element(c0), Fp2Element(c1), Fp2Element(c2))
Fp6Element(c0) = Fp6Element(Fp2Element(c0), zero(Fp2Element), zero(Fp2Element))

# Access components
Base.getindex(a::Fp6Element, i::Int) = a.coeffs[i]

# Basic operations
Base.zero(::Type{Fp6Element}) = Fp6Element(zero(Fp2Element), zero(Fp2Element), zero(Fp2Element))
Base.one(::Type{Fp6Element}) = Fp6Element(one(Fp2Element), zero(Fp2Element), zero(Fp2Element))
Base.zero(::Fp6Element) = zero(Fp6Element)
Base.one(::Fp6Element) = one(Fp6Element)

Base.iszero(a::Fp6Element) = all(iszero, a.coeffs)
Base.isone(a::Fp6Element) = isone(a[1]) && iszero(a[2]) && iszero(a[3])

# Equality
Base.:(==)(a::Fp6Element, b::Fp6Element) = a.coeffs == b.coeffs

# Addition and subtraction
Base.:+(a::Fp6Element, b::Fp6Element) = Fp6Element(a[1] + b[1], a[2] + b[2], a[3] + b[3])
Base.:-(a::Fp6Element, b::Fp6Element) = Fp6Element(a[1] - b[1], a[2] - b[2], a[3] - b[3])
Base.:-(a::Fp6Element) = Fp6Element(-a[1], -a[2], -a[3])

"""
    *(a::Fp6Element, b::Fp6Element)

Multiply two `Fp6Element`s using the relation `v³ = ξ = 9 + u`.
The product expands as `(a₀ + a₁v + a₂v²)(b₀ + b₁v + b₂v²)` with Karatsuba
style recombinations.
"""
function Base.:*(a::Fp6Element, b::Fp6Element)
    # ξ = 9 + u in Fp2 (non-residue for cubic extension)
    ξ = Fp2Element(9, 1)

    # Use Karatsuba-like formula to reduce multiplications
    v0 = a[1] * b[1]
    v1 = a[2] * b[2]
    v2 = a[3] * b[3]

    # Compute cross terms
    c0 = v0 + ((a[2] + a[3]) * (b[2] + b[3]) - v1 - v2) * ξ
    c1 = (a[1] + a[2]) * (b[1] + b[2]) - v0 - v1 + v2 * ξ
    c2 = (a[1] + a[3]) * (b[1] + b[3]) - v0 - v2 + v1

    return Fp6Element(c0, c1, c2)
end

"""
    square(a::Fp6Element)

Square an `Fp6Element` using the optimized sextic formula.
"""
function square(a::Fp6Element)
    ξ = Fp2Element(9, 1)

    c0 = a[1]^2 + Fp2Element(2) * a[2] * a[3] * ξ
    c1 = Fp2Element(2) * a[1] * a[2] + a[3]^2 * ξ
    c2 = Fp2Element(2) * a[1] * a[3] + a[2]^2

    return Fp6Element(c0, c1, c2)
end

Base.:^(a::Fp6Element, ::Val{2}) = square(a)

"""
    inv(a::Fp6Element)

Compute the multiplicative inverse in Fp6.
"""
function Base.inv(a::Fp6Element)
    if iszero(a)
        throw(DivideError())
    end

    ξ = Fp2Element(9, 1)

    # For a = a0 + a1*v + a2*v² in Fp6
    # We compute the inverse using the norm map to Fp2

    # Compute intermediate values
    a0_sq = a[1]^2
    a1_sq = a[2]^2
    a2_sq = a[3]^2

    # v0 = a0² - a1*a2*ξ
    v0 = a0_sq - a[2] * a[3] * ξ

    # v1 = a2²*ξ - a0*a1
    v1 = a2_sq * ξ - a[1] * a[2]

    # v2 = a1² - a0*a2
    v2 = a1_sq - a[1] * a[3]

    # Norm = a0*v0 + a1*v2*ξ + a2*v1*ξ
    norm = a[1] * v0 + a[2] * v2 * ξ + a[3] * v1 * ξ
    norm_inv = inv(norm)

    # The inverse is (v0, v1, v2) / norm
    return Fp6Element(v0 * norm_inv, v1 * norm_inv, v2 * norm_inv)
end

Base.:/(a::Fp6Element, b::Fp6Element) = a * inv(b)

"""
    ^(a::Fp6Element, n::Integer)

Raise an `Fp6Element` to `n` with binary exponentiation.
"""
function Base.:^(a::Fp6Element, n::Integer)
    if n == 0
        return one(Fp6Element)
    elseif n == 2
        return square(a)
    elseif n < 0
        return inv(a)^(-n)
    end

    result = one(Fp6Element)
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

# Scalar multiplication for convenience
Base.:*(k::Integer, a::Fp6Element) = Fp6Element(k * a[1], k * a[2], k * a[3])
Base.:*(a::Fp6Element, k::Integer) = k * a

# Display
function Base.show(io::IO, a::Fp6Element)
    print(io, "Fp6([", a[1], ", ", a[2], ", ", a[3], "])")
end

# Export types and functions
export Fp6Element, square
