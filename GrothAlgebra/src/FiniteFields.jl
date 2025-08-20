"""
Finite field implementations using BigInt for arbitrary precision.

This module provides finite field arithmetic for various prime fields,
using BigInt internally to avoid LLVM issues with large integer types.
"""

# Abstract type for all finite fields
abstract type FiniteFieldElement end

# Helper functions that must be implemented by concrete types
prime(::Type{<:FiniteFieldElement}) = error("prime() not implemented")
field_name(::Type{<:FiniteFieldElement}) = "FiniteField"

# ===========================================
# BN254 Field Implementation
# ===========================================

struct BN254Field <: FiniteFieldElement
    value::BigInt
    
    # Inner constructor ensures normalization
    function BN254Field(value::BigInt, normalized::Bool=false)
        if normalized
            new(value)
        else
            p = prime(BN254Field)
            new(mod(value, p))
        end
    end
end

# BN254 prime
prime(::Type{BN254Field}) = parse(BigInt, "21888242871839275222246405745257275088696311157297823662689037894645226208583")
field_name(::Type{BN254Field}) = "BN254"

# Constructors
BN254Field(x::Integer) = BN254Field(BigInt(x))
bn254_field(x) = BN254Field(x)

# ===========================================
# Secp256k1 Field Implementation  
# ===========================================

struct Secp256k1Field <: FiniteFieldElement
    value::BigInt
    
    function Secp256k1Field(value::BigInt, normalized::Bool=false)
        if normalized
            new(value)
        else
            p = prime(Secp256k1Field)
            new(mod(value, p))
        end
    end
end

# Secp256k1 prime: 2^256 - 2^32 - 977
prime(::Type{Secp256k1Field}) = parse(BigInt, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F", base=16)
field_name(::Type{Secp256k1Field}) = "Secp256k1"

# Constructors
Secp256k1Field(x::Integer) = Secp256k1Field(BigInt(x))
secp256k1_field(x) = Secp256k1Field(x)

# ===========================================
# Generic Field Operations
# ===========================================

# Zero and one
Base.zero(::Type{T}) where {T<:FiniteFieldElement} = T(BigInt(0), true)
Base.one(::Type{T}) where {T<:FiniteFieldElement} = T(BigInt(1), true)
Base.zero(x::T) where {T<:FiniteFieldElement} = zero(T)
Base.one(x::T) where {T<:FiniteFieldElement} = one(T)

# Predicates
Base.iszero(x::FiniteFieldElement) = x.value == 0
Base.isone(x::FiniteFieldElement) = x.value == 1
is_zero(x::FiniteFieldElement) = iszero(x)
is_one(x::FiniteFieldElement) = isone(x)
is_unity(x::FiniteFieldElement) = isone(x)

# Comparison
Base.:(==)(x::T, y::T) where {T<:FiniteFieldElement} = x.value == y.value
Base.isequal(x::T, y::T) where {T<:FiniteFieldElement} = x.value == y.value

# Arithmetic operations
function Base.:+(x::T, y::T) where {T<:FiniteFieldElement}
    p = prime(T)
    result = mod(x.value + y.value, p)
    return T(result, true)
end

function Base.:-(x::T, y::T) where {T<:FiniteFieldElement}
    p = prime(T)
    result = mod(x.value - y.value, p)
    return T(result, true)
end

function Base.:-(x::T) where {T<:FiniteFieldElement}
    if iszero(x)
        return x
    end
    p = prime(T)
    return T(p - x.value, true)
end

function Base.:*(x::T, y::T) where {T<:FiniteFieldElement}
    p = prime(T)
    result = mod(x.value * y.value, p)
    return T(result, true)
end

function Base.:*(x::T, n::Integer) where {T<:FiniteFieldElement}
    p = prime(T)
    result = mod(x.value * n, p)
    return T(result, true)
end

Base.:*(n::Integer, x::T) where {T<:FiniteFieldElement} = x * n

function Base.inv(x::T) where {T<:FiniteFieldElement}
    if iszero(x)
        throw(DivideError())
    end
    p = prime(T)
    # Use built-in modular inverse
    return T(invmod(x.value, p), true)
end

Base.:/(x::T, y::T) where {T<:FiniteFieldElement} = x * inv(y)

function Base.:^(x::T, n::Integer) where {T<:FiniteFieldElement}
    if n == 0
        return one(T)
    elseif n < 0
        return inv(x)^(-n)
    end
    p = prime(T)
    # Use built-in modular exponentiation
    return T(powermod(x.value, n, p), true)
end

# Display
function Base.show(io::IO, x::T) where {T<:FiniteFieldElement}
    print(io, field_name(T), "(", x.value, ")")
end

# Conversion
Base.convert(::Type{BigInt}, x::FiniteFieldElement) = x.value
Base.convert(::Type{Integer}, x::FiniteFieldElement) = x.value

# Make types compatible with FieldElem interface
const FieldElem = FiniteFieldElement

# Export everything
export FiniteFieldElement, FieldElem
export BN254Field, bn254_field  
export Secp256k1Field, secp256k1_field
export prime, is_zero, is_one, is_unity