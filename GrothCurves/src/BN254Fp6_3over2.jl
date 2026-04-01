# BN254 Fp6 implementation as cubic extension of Fp2.
#
# Realises Fp6 = Fp2[v]/(v³ - ξ) with ξ = 9 + u ∈ Fp2.

"""
    Fp6Element

Element of sextic extension `Fp6 = Fp2[v] / (v^3 - xi)` with `xi = 9 + u`.
Represented as `c0 + c1*v + c2*v^2` with `c0, c1, c2 in Fp2`.
"""
struct Fp6Element
    c0::Fp2Element
    c1::Fp2Element
    c2::Fp2Element

    function Fp6Element(c0::Fp2Element, c1::Fp2Element, c2::Fp2Element)
        new(c0, c1, c2)
    end
end

# Convenient constructors
Fp6Element(c0, c1, c2) = Fp6Element(Fp2Element(c0), Fp2Element(c1), Fp2Element(c2))
Fp6Element(c0::Fp2Element) = Fp6Element(c0, zero(Fp2Element), zero(Fp2Element))
Fp6Element(c0) = Fp6Element(Fp2Element(c0), zero(Fp2Element), zero(Fp2Element))

const FQ6_NONRESIDUE = Fp2Element(9, 1)

"""
    mul_fp2_by_nonresidue(a::Fp2Element)

Multiply `a` by the sextic-twist nonresidue `ξ = 9 + u`.
"""
function mul_fp2_by_nonresidue(a::Fp2Element)
    c0_x2 = a.c0 + a.c0
    c0_x4 = c0_x2 + c0_x2
    c0_x8 = c0_x4 + c0_x4
    c1_x2 = a.c1 + a.c1
    c1_x4 = c1_x2 + c1_x2
    c1_x8 = c1_x4 + c1_x4
    return Fp2Element(c0_x8 + a.c0 - a.c1, c1_x8 + a.c1 + a.c0)
end

# Access components
function Base.getindex(a::Fp6Element, i::Int)
    if i == 1
        return a.c0
    elseif i == 2
        return a.c1
    elseif i == 3
        return a.c2
    end
    throw(BoundsError(a, i))
end

# Basic operations
Base.zero(::Type{Fp6Element}) = Fp6Element(zero(Fp2Element), zero(Fp2Element), zero(Fp2Element))
Base.one(::Type{Fp6Element}) = Fp6Element(one(Fp2Element), zero(Fp2Element), zero(Fp2Element))
Base.zero(::Fp6Element) = zero(Fp6Element)
Base.one(::Fp6Element) = one(Fp6Element)

Base.iszero(a::Fp6Element) = iszero(a.c0) && iszero(a.c1) && iszero(a.c2)
Base.isone(a::Fp6Element) = isone(a.c0) && iszero(a.c1) && iszero(a.c2)

# Equality
Base.:(==)(a::Fp6Element, b::Fp6Element) = a.c0 == b.c0 && a.c1 == b.c1 && a.c2 == b.c2
Base.isequal(a::Fp6Element, b::Fp6Element) = isequal(a.c0, b.c0) && isequal(a.c1, b.c1) && isequal(a.c2, b.c2)

# Addition and subtraction
Base.:+(a::Fp6Element, b::Fp6Element) = Fp6Element(a.c0 + b.c0, a.c1 + b.c1, a.c2 + b.c2)
Base.:-(a::Fp6Element, b::Fp6Element) = Fp6Element(a.c0 - b.c0, a.c1 - b.c1, a.c2 - b.c2)
Base.:-(a::Fp6Element) = Fp6Element(-a.c0, -a.c1, -a.c2)

"""
    *(a::Fp6Element, b::Fp6Element)

Multiply two `Fp6Element`s using the relation `v³ = ξ = 9 + u`.
The product expands as `(a₀ + a₁v + a₂v²)(b₀ + b₁v + b₂v²)` with Karatsuba
style recombinations.
"""
function Base.:*(a::Fp6Element, b::Fp6Element)
    v0 = a.c0 * b.c0
    v1 = a.c1 * b.c1
    v2 = a.c2 * b.c2

    c0 = v0 + mul_fp2_by_nonresidue((a.c1 + a.c2) * (b.c1 + b.c2) - v1 - v2)
    c1 = (a.c0 + a.c1) * (b.c0 + b.c1) - v0 - v1 + mul_fp2_by_nonresidue(v2)
    c2 = (a.c0 + a.c2) * (b.c0 + b.c2) - v0 - v2 + v1

    return Fp6Element(c0, c1, c2)
end

"""
    square(a::Fp6Element)

Square an `Fp6Element` using the optimized sextic formula.
"""
function square(a::Fp6Element)
    c1c2 = a.c1 * a.c2
    c0c1 = a.c0 * a.c1
    c0c2 = a.c0 * a.c2
    c0 = square(a.c0) + mul_fp2_by_nonresidue(c1c2 + c1c2)
    c1 = (c0c1 + c0c1) + mul_fp2_by_nonresidue(square(a.c2))
    c2 = (c0c2 + c0c2) + square(a.c1)
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

    a0_sq = square(a.c0)
    a1_sq = square(a.c1)
    a2_sq = square(a.c2)

    v0 = a0_sq - mul_fp2_by_nonresidue(a.c1 * a.c2)
    v1 = mul_fp2_by_nonresidue(a2_sq) - a.c0 * a.c1
    v2 = a1_sq - a.c0 * a.c2

    norm = a.c0 * v0 + mul_fp2_by_nonresidue(a.c1 * v2 + a.c2 * v1)
    norm_inv = inv(norm)

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
Base.:*(k::Integer, a::Fp6Element) = Fp6Element(k * a.c0, k * a.c1, k * a.c2)
Base.:*(a::Fp6Element, k::Integer) = k * a

# Display
function Base.show(io::IO, a::Fp6Element)
    print(io, "Fp6([", a.c0, ", ", a.c1, ", ", a.c2, "])")
end

# Export types and functions
export Fp6Element, FQ6_NONRESIDUE, mul_fp2_by_nonresidue, square
