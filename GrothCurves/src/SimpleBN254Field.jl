"""
Simplified BN254 field implementation using BigInt to avoid LLVM issues.
"""

# BN254 prime
const BN254_PRIME_VALUE = BigInt("21888242871839275222246405745257275088696311157297823662689037894645226208583")

"""
    BN254FieldElement

A field element in the BN254 prime field.
"""
struct BN254FieldElement
    value::BigInt
    
    function BN254FieldElement(value::BigInt, normalized::Bool=false)
        if normalized
            new(value)
        else
            new(mod(value, BN254_PRIME_VALUE))
        end
    end
end

# Constructors
BN254FieldElement(value::Integer) = BN254FieldElement(BigInt(value))

# Basic operations
Base.zero(::Type{BN254FieldElement}) = BN254FieldElement(BigInt(0), true)
Base.one(::Type{BN254FieldElement}) = BN254FieldElement(BigInt(1), true)
Base.iszero(a::BN254FieldElement) = a.value == 0
Base.isone(a::BN254FieldElement) = a.value == 1

# Arithmetic
Base.:+(a::BN254FieldElement, b::BN254FieldElement) = BN254FieldElement(mod(a.value + b.value, BN254_PRIME_VALUE), true)
Base.:-(a::BN254FieldElement, b::BN254FieldElement) = BN254FieldElement(mod(a.value - b.value, BN254_PRIME_VALUE), true)
Base.:-(a::BN254FieldElement) = BN254FieldElement(mod(-a.value, BN254_PRIME_VALUE), true)
Base.:*(a::BN254FieldElement, b::BN254FieldElement) = BN254FieldElement(mod(a.value * b.value, BN254_PRIME_VALUE), true)

function Base.inv(a::BN254FieldElement)
    if iszero(a)
        throw(DivideError())
    end
    # Extended Euclidean algorithm
    return BN254FieldElement(invmod(a.value, BN254_PRIME_VALUE), true)
end

Base.:/(a::BN254FieldElement, b::BN254FieldElement) = a * inv(b)
Base.:^(a::BN254FieldElement, n::Integer) = BN254FieldElement(powermod(a.value, n, BN254_PRIME_VALUE), true)

# Comparison
Base.:(==)(a::BN254FieldElement, b::BN254FieldElement) = a.value == b.value

# Display
Base.show(io::IO, a::BN254FieldElement) = print(io, "BN254(", a.value, ")")

# Make it compatible with FieldElem interface
const BN254Field = BN254FieldElement
bn254_field(x) = BN254FieldElement(x)

# Export
export BN254FieldElement, BN254Field, bn254_field, BN254_PRIME_VALUE