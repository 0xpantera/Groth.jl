"""
BN254 curve field implementations.

The BN254 curve (also known as alt-bn128) is defined over:
- Base field Fp with p = 21888242871839275222246405745257275088696311157297823662689037894645226208583
- Extension field Fp2 = Fp[u]/(u² + 1)
"""

using GrothAlgebra
using StaticArrays

# BN254Field is now provided by GrothAlgebra
# Just re-export what we need
const BN254_PRIME = prime(BN254Field)

# Note: SimpleBN254Field.jl is kept for reference but not used
# It contains an alternative BigInt-based implementation created
# during the migration from UInt256 to BigInt

"""
    Fp2Element

Element of the quadratic extension field Fp2 = Fp[u]/(u² + 1).
Represented as c0 + c1*u where c0, c1 ∈ Fp.
"""
struct Fp2Element
    coeffs::SVector{2, BN254Field}
    
    function Fp2Element(c0::BN254Field, c1::BN254Field)
        new(SVector(c0, c1))
    end
end

# Convenient constructors
Fp2Element(c0::Integer, c1::Integer) = Fp2Element(bn254_field(c0), bn254_field(c1))
Fp2Element(c0) = Fp2Element(bn254_field(c0), zero(BN254Field))

# Access components
Base.getindex(a::Fp2Element, i::Int) = a.coeffs[i]
real(a::Fp2Element) = a.coeffs[1]
imag(a::Fp2Element) = a.coeffs[2]

# Basic operations
Base.zero(::Type{Fp2Element}) = Fp2Element(zero(BN254Field), zero(BN254Field))
Base.one(::Type{Fp2Element}) = Fp2Element(one(BN254Field), zero(BN254Field))
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

Multiplication in Fp2 using the relation u² = -1.
(a0 + a1*u) * (b0 + b1*u) = (a0*b0 - a1*b1) + (a0*b1 + a1*b0)*u
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
Base.:*(a::Fp2Element, k::BN254Field) = Fp2Element(a[1] * k, a[2] * k)
Base.:*(k::BN254Field, a::Fp2Element) = a * k

"""
    conjugate(a::Fp2Element)

Conjugate of a + b*u is a - b*u.
"""
conjugate(a::Fp2Element) = Fp2Element(a[1], -a[2])

"""
    norm(a::Fp2Element)

Norm of a + b*u is a² + b².
"""
norm(a::Fp2Element) = a[1]^2 + a[2]^2

"""
    inv(a::Fp2Element)

Multiplicative inverse using the formula:
(a + b*u)⁻¹ = (a - b*u) / (a² + b²)
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

Exponentiation using binary method.
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

Frobenius endomorphism: (a + b*u) -> (a - b*u) since u^p = -u in Fp2.
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
export BN254Field, bn254_field, BN254_PRIME
export Fp2Element, conjugate, norm, frobenius, real, imag