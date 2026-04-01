"""
    FiniteFieldElement

Abstract supertype for prime-field elements.

Concrete field types may choose their own internal representation, but they must
support the internal backend hooks defined in this file.
"""
abstract type FiniteFieldElement end

using Primes: isprime

# Helper functions that must be implemented by concrete types
prime(::Type{<:FiniteFieldElement}) = error("prime() not implemented")
field_name(::Type{<:FiniteFieldElement}) = "FiniteField"

# Backend hooks

"""
    canonical_bigint(x::FiniteFieldElement)

Return the canonical non-negative integer representative of a field element.
The current default backend stores this directly as `BigInt`, but future
backends may override this hook.
"""
canonical_bigint(x::FiniteFieldElement) = getfield(x, :value)

"""
    construct_normalized(::Type{T}, value::BigInt) where T<:FiniteFieldElement

Construct a field element from an already-normalized canonical representative.
"""
construct_normalized(::Type{T}, value::BigInt) where {T<:FiniteFieldElement} = T(value, true)

"""
    field_zero_raw(::Type{T}) where T<:FiniteFieldElement

Canonical raw zero for the field backend.
"""
field_zero_raw(::Type{T}) where {T<:FiniteFieldElement} = BigInt(0)

"""
    field_one_raw(::Type{T}) where T<:FiniteFieldElement

Canonical raw one for the field backend.
"""
field_one_raw(::Type{T}) where {T<:FiniteFieldElement} = BigInt(1)

"""
    field_equal(x::T, y::T) where T<:FiniteFieldElement

Backend-aware equality for field elements.
"""
field_equal(x::T, y::T) where {T<:FiniteFieldElement} = canonical_bigint(x) == canonical_bigint(y)

"""
    field_add(x::T, y::T) where T<:FiniteFieldElement

Backend-aware field addition.
"""
function field_add(x::T, y::T) where {T<:FiniteFieldElement}
    result = mod(canonical_bigint(x) + canonical_bigint(y), prime(T))
    return construct_normalized(T, result)
end

"""
    field_sub(x::T, y::T) where T<:FiniteFieldElement

Backend-aware field subtraction.
"""
function field_sub(x::T, y::T) where {T<:FiniteFieldElement}
    result = mod(canonical_bigint(x) - canonical_bigint(y), prime(T))
    return construct_normalized(T, result)
end

"""
    field_neg(x::T) where T<:FiniteFieldElement

Backend-aware field negation.
"""
function field_neg(x::T) where {T<:FiniteFieldElement}
    iszero(x) && return x
    return construct_normalized(T, prime(T) - canonical_bigint(x))
end

"""
    field_mul(x::T, y::T) where T<:FiniteFieldElement

Backend-aware field multiplication.
"""
function field_mul(x::T, y::T) where {T<:FiniteFieldElement}
    result = mod(canonical_bigint(x) * canonical_bigint(y), prime(T))
    return construct_normalized(T, result)
end

"""
    field_mul_int(x::T, n::Integer) where T<:FiniteFieldElement

Backend-aware scalar multiplication by an integer.
"""
function field_mul_int(x::T, n::Integer) where {T<:FiniteFieldElement}
    result = mod(canonical_bigint(x) * n, prime(T))
    return construct_normalized(T, result)
end

"""
    field_inv(x::T) where T<:FiniteFieldElement

Backend-aware multiplicative inverse.
"""
function field_inv(x::T) where {T<:FiniteFieldElement}
    iszero(x) && throw(DivideError())
    return construct_normalized(T, invmod(canonical_bigint(x), prime(T)))
end

"""
    field_pow(x::T, n::Integer) where T<:FiniteFieldElement

Backend-aware exponentiation.
"""
function field_pow(x::T, n::Integer) where {T<:FiniteFieldElement}
    if n == 0
        return one(T)
    elseif n < 0
        return inv(x)^(-n)
    end
    return construct_normalized(T, powermod(canonical_bigint(x), n, prime(T)))
end

# ===========================================
# Generic Galois Field (Prime Only)
# ===========================================

const _galois_field_prime_cache = Dict{Int,Bool}()

_checked_prime_modulus(p) = throw(ArgumentError("GaloisField modulus must be an integer > 1"))

function _checked_prime_modulus(p::Integer)
    p_int = try
        Int(p)
    catch
        throw(ArgumentError("GaloisField modulus must fit in Int"))
    end
    if p_int <= 1 || p_int != p
        throw(ArgumentError("GaloisField modulus must be an integer > 1"))
    end
    is_prime = get!(_galois_field_prime_cache, p_int) do
        isprime(p_int)
    end
    is_prime || throw(ArgumentError("GaloisField modulus must be prime"))
    return p_int
end

"""
    GaloisField{P}

Prime-field element in GF(P) using the current default `BigInt` backend.
Currently supports only prime fields (P must be prime).
"""
struct GaloisField{P} <: FiniteFieldElement
    value::BigInt

    function GaloisField{P}(value::BigInt, normalized::Bool=false) where {P}
        _checked_prime_modulus(P)
        if normalized
            new(value)
        else
            p = BigInt(P)
            new(mod(value, p))
        end
    end
end

prime(::Type{GaloisField{P}}) where {P} = BigInt(P)
field_name(::Type{GaloisField{P}}) where {P} = "GF($(P))"

GaloisField{P}(x::Integer) where {P} = GaloisField{P}(BigInt(x))

"""
    galois_field(p::Integer)

Return the prime field type `GaloisField{p}` after validating `p` is prime.
"""
function galois_field(p::Integer)
    p_int = _checked_prime_modulus(p)
    return GaloisField{p_int}
end

# ===========================================
# BN254 Field Implementation
# ===========================================

"""
    BN254Fq

Element of the BN254 base field using the current default `BigInt` backend.
"""
struct BN254Fq <: FiniteFieldElement
    value::BigInt

    # Inner constructor ensures normalization
    function BN254Fq(value::BigInt, normalized::Bool=false)
        if normalized
            new(value)
        else
            p = prime(BN254Fq)
            new(mod(value, p))
        end
    end
end

# BN254 prime
prime(::Type{BN254Fq}) = parse(BigInt, "21888242871839275222246405745257275088696311157297823662689037894645226208583")
field_name(::Type{BN254Fq}) = "BN254"

# Constructors
BN254Fq(x::Integer) = BN254Fq(BigInt(x))
bn254_fq(x) = BN254Fq(x)

# ===========================================
# Secp256k1 Field Implementation  
# ===========================================

"""
    Secp256k1Field

Element of the secp256k1 base field using the current default `BigInt`
backend.
"""
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
Base.zero(::Type{T}) where {T<:FiniteFieldElement} = construct_normalized(T, field_zero_raw(T))
Base.one(::Type{T}) where {T<:FiniteFieldElement} = construct_normalized(T, field_one_raw(T))
Base.zero(x::T) where {T<:FiniteFieldElement} = zero(T)
Base.one(x::T) where {T<:FiniteFieldElement} = one(T)

# Predicates
Base.iszero(x::FiniteFieldElement) = canonical_bigint(x) == 0
Base.isone(x::FiniteFieldElement) = canonical_bigint(x) == 1
is_zero(x::FiniteFieldElement) = iszero(x)
is_one(x::FiniteFieldElement) = isone(x)
is_unity(x::FiniteFieldElement) = isone(x)

# Comparison
Base.:(==)(x::T, y::T) where {T<:FiniteFieldElement} = field_equal(x, y)
Base.isequal(x::T, y::T) where {T<:FiniteFieldElement} = field_equal(x, y)

# Arithmetic operations
Base.:+(x::T, y::T) where {T<:FiniteFieldElement} = field_add(x, y)
Base.:-(x::T, y::T) where {T<:FiniteFieldElement} = field_sub(x, y)
Base.:-(x::T) where {T<:FiniteFieldElement} = field_neg(x)
Base.:*(x::T, y::T) where {T<:FiniteFieldElement} = field_mul(x, y)
Base.:*(x::T, n::Integer) where {T<:FiniteFieldElement} = field_mul_int(x, n)

Base.:*(n::Integer, x::T) where {T<:FiniteFieldElement} = x * n

Base.inv(x::T) where {T<:FiniteFieldElement} = field_inv(x)

Base.:/(x::T, y::T) where {T<:FiniteFieldElement} = x * inv(y)

Base.:^(x::T, n::Integer) where {T<:FiniteFieldElement} = field_pow(x, n)

# Display
function Base.show(io::IO, x::T) where {T<:FiniteFieldElement}
    print(io, field_name(T), "(", canonical_bigint(x), ")")
end

# Conversion
Base.convert(::Type{BigInt}, x::FiniteFieldElement) = canonical_bigint(x)

# Make types compatible with FieldElem interface
const FieldElem = FiniteFieldElement

# Export everything
export FiniteFieldElement, FieldElem
export GaloisField, galois_field
export BN254Fq, bn254_fq
export Secp256k1Field, secp256k1_field
export prime, is_zero, is_one, is_unity

# ===========================================
# BN254 Scalar Field (Fr) Implementation
# ===========================================

"""
    BN254Fr

Element of the BN254 scalar field (Fr) used for curve scalars.
"""
struct BN254Fr <: FiniteFieldElement
    value::BigInt

    function BN254Fr(value::BigInt, normalized::Bool=false)
        if normalized
            new(value)
        else
            p = prime(BN254Fr)
            new(mod(value, p))
        end
    end
end

prime(::Type{BN254Fr}) = parse(BigInt, "21888242871839275222246405745257275088548364400416034343698204186575808495617")
field_name(::Type{BN254Fr}) = "BN254Fr"

BN254Fr(x::Integer) = BN254Fr(BigInt(x))
bn254_fr(x) = BN254Fr(x)

export BN254Fr, bn254_fr
