"""
    FiniteFieldElement

Abstract supertype for prime-field elements.

Concrete field types may choose their own internal representation, but they must
support the internal backend hooks defined in this file.
"""
abstract type FiniteFieldElement end

using Primes: isprime
using StaticArrays: MVector

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
# Montgomery prime-field helpers
# ===========================================

abstract type MontgomeryFiniteFieldElement <: FiniteFieldElement end

const MONT_WORD_BITS = 64
const MONT_LIMB_COUNT = 4
const MONT_WORD_MASK = (BigInt(1) << MONT_WORD_BITS) - 1
const MONT_WORD_MASK_U128 = (UInt128(1) << MONT_WORD_BITS) - UInt128(1)
const MONT_RADIX = BigInt(1) << (MONT_WORD_BITS * MONT_LIMB_COUNT)
const CANONICAL_ONE_LIMBS = (UInt64(1), UInt64(0), UInt64(0), UInt64(0))

montgomery_modulus_limbs(::Type{<:MontgomeryFiniteFieldElement}) = error("montgomery_modulus_limbs() not implemented")
montgomery_r2_limbs(::Type{<:MontgomeryFiniteFieldElement}) = error("montgomery_r2_limbs() not implemented")
montgomery_one_limbs(::Type{<:MontgomeryFiniteFieldElement}) = error("montgomery_one_limbs() not implemented")
montgomery_n0inv(::Type{<:MontgomeryFiniteFieldElement}) = error("montgomery_n0inv() not implemented")

@inline zero_limbs() = (UInt64(0), UInt64(0), UInt64(0), UInt64(0))
@inline low_u64(x::UInt128) = UInt64(x & MONT_WORD_MASK_U128)

function bigint_to_limbs(x::BigInt)::NTuple{4,UInt64}
    y = x
    return ntuple(MONT_LIMB_COUNT) do _
        limb = UInt64(y & MONT_WORD_MASK)
        y >>= MONT_WORD_BITS
        limb
    end
end

function limbs_to_bigint(limbs::NTuple{4,UInt64})::BigInt
    acc = BigInt(0)
    for i in MONT_LIMB_COUNT:-1:1
        acc <<= MONT_WORD_BITS
        acc += limbs[i]
    end
    return acc
end

function montgomery_constants(modulus::BigInt)
    modulus_limbs = bigint_to_limbs(modulus)
    n0inv = UInt64(mod(-invmod(BigInt(modulus_limbs[1]), BigInt(1) << MONT_WORD_BITS), BigInt(1) << MONT_WORD_BITS))
    one_limbs = bigint_to_limbs(mod(MONT_RADIX, modulus))
    r2_limbs = bigint_to_limbs(mod(MONT_RADIX * MONT_RADIX, modulus))
    return modulus_limbs, one_limbs, r2_limbs, n0inv
end

@inline function limbs_gte(lhs::NTuple{4,UInt64}, rhs::NTuple{4,UInt64})::Bool
    @inbounds for i in MONT_LIMB_COUNT:-1:1
        lhs_i = lhs[i]
        rhs_i = rhs[i]
        lhs_i == rhs_i && continue
        return lhs_i > rhs_i
    end
    return true
end

@inline function subborrow_u64(a::UInt64, b::UInt64, borrow::UInt64)
    rhs = UInt128(b) + UInt128(borrow)
    if UInt128(a) >= rhs
        return UInt64(UInt128(a) - rhs), UInt64(0)
    end
    return UInt64((UInt128(1) << MONT_WORD_BITS) + UInt128(a) - rhs), UInt64(1)
end

@inline function add_limbs(lhs::NTuple{4,UInt64}, rhs::NTuple{4,UInt64})::Tuple{NTuple{4,UInt64},UInt64}
    out = MVector{4,UInt64}(ntuple(_ -> UInt64(0), 4))
    carry = UInt64(0)
    @inbounds for i in 1:MONT_LIMB_COUNT
        wide = UInt128(lhs[i]) + UInt128(rhs[i]) + UInt128(carry)
        out[i] = low_u64(wide)
        carry = UInt64(wide >> MONT_WORD_BITS)
    end
    return (out[1], out[2], out[3], out[4]), carry
end

@inline function sub_limbs(lhs::NTuple{4,UInt64}, rhs::NTuple{4,UInt64})::Tuple{NTuple{4,UInt64},UInt64}
    out = MVector{4,UInt64}(ntuple(_ -> UInt64(0), 4))
    borrow = UInt64(0)
    @inbounds for i in 1:MONT_LIMB_COUNT
        out[i], borrow = subborrow_u64(lhs[i], rhs[i], borrow)
    end
    return (out[1], out[2], out[3], out[4]), borrow
end

@inline function add_mod(lhs::NTuple{4,UInt64}, rhs::NTuple{4,UInt64}, modulus::NTuple{4,UInt64})::NTuple{4,UInt64}
    sum_limbs, _ = add_limbs(lhs, rhs)
    if limbs_gte(sum_limbs, modulus)
        return first(sub_limbs(sum_limbs, modulus))
    end
    return sum_limbs
end

@inline function sub_mod(lhs::NTuple{4,UInt64}, rhs::NTuple{4,UInt64}, modulus::NTuple{4,UInt64})::NTuple{4,UInt64}
    diff_limbs, borrow = sub_limbs(lhs, rhs)
    if borrow != 0
        return first(add_limbs(diff_limbs, modulus))
    end
    return diff_limbs
end

@inline function mul_4x4(lhs::NTuple{4,UInt64}, rhs::NTuple{4,UInt64})
    t = MVector{10,UInt64}(ntuple(_ -> UInt64(0), 10))
    @inbounds for i in 1:MONT_LIMB_COUNT
        carry = UInt128(0)
        lhs_i = lhs[i]
        for j in 1:MONT_LIMB_COUNT
            idx = i + j - 1
            prod = UInt128(lhs_i) * UInt128(rhs[j])
            sum = UInt128(t[idx]) + UInt128(low_u64(prod)) + carry
            t[idx] = low_u64(sum)
            carry = UInt128(prod >> MONT_WORD_BITS) + (sum >> MONT_WORD_BITS)
        end
        idx = i + MONT_LIMB_COUNT
        sum = UInt128(t[idx]) + carry
        t[idx] = low_u64(sum)
        carry = sum >> MONT_WORD_BITS
        while carry != 0
            idx += 1
            sum = UInt128(t[idx]) + carry
            t[idx] = low_u64(sum)
            carry = sum >> MONT_WORD_BITS
        end
    end
    return t
end

@inline function montgomery_reduce(::Type{T}, t::MVector{10,UInt64}) where {T<:MontgomeryFiniteFieldElement}
    modulus = montgomery_modulus_limbs(T)
    n0inv = montgomery_n0inv(T)

    @inbounds for i in 1:MONT_LIMB_COUNT
        m = t[i] * n0inv
        carry = UInt128(0)
        for j in 1:MONT_LIMB_COUNT
            idx = i + j - 1
            prod = UInt128(m) * UInt128(modulus[j])
            sum = UInt128(t[idx]) + UInt128(low_u64(prod)) + carry
            t[idx] = low_u64(sum)
            carry = UInt128(prod >> MONT_WORD_BITS) + (sum >> MONT_WORD_BITS)
        end
        idx = i + MONT_LIMB_COUNT
        sum = UInt128(t[idx]) + carry
        t[idx] = low_u64(sum)
        carry = sum >> MONT_WORD_BITS
        while carry != 0
            idx += 1
            sum = UInt128(t[idx]) + carry
            t[idx] = low_u64(sum)
            carry = sum >> MONT_WORD_BITS
        end
    end

    reduced = (t[5], t[6], t[7], t[8])
    if limbs_gte(reduced, modulus)
        return first(sub_limbs(reduced, modulus))
    end
    return reduced
end

@inline function montgomery_mul_limbs(::Type{T}, lhs::NTuple{4,UInt64}, rhs::NTuple{4,UInt64}) where {T<:MontgomeryFiniteFieldElement}
    return montgomery_reduce(T, mul_4x4(lhs, rhs))
end

@inline montgomery_square_limbs(::Type{T}, lhs::NTuple{4,UInt64}) where {T<:MontgomeryFiniteFieldElement} =
    montgomery_mul_limbs(T, lhs, lhs)

@inline function montgomery_encode_limbs(::Type{T}, canonical::NTuple{4,UInt64}) where {T<:MontgomeryFiniteFieldElement}
    return montgomery_mul_limbs(T, canonical, montgomery_r2_limbs(T))
end

@inline function montgomery_decode_limbs(::Type{T}, montgomery::NTuple{4,UInt64}) where {T<:MontgomeryFiniteFieldElement}
    return montgomery_mul_limbs(T, montgomery, CANONICAL_ONE_LIMBS)
end

@inline montgomery_storage(x::MontgomeryFiniteFieldElement) = getfield(x, :limbs)

function normalize_canonical(::Type{T}, value::BigInt, normalized::Bool)::BigInt where {T<:MontgomeryFiniteFieldElement}
    if normalized && value >= 0 && value < prime(T)
        return value
    end
    return mod(value, prime(T))
end

construct_montgomery(::Type{T}, limbs::NTuple{4,UInt64}) where {T<:MontgomeryFiniteFieldElement} = T(limbs)

function construct_normalized(::Type{T}, value::BigInt) where {T<:MontgomeryFiniteFieldElement}
    canonical = normalize_canonical(T, value, true)
    return construct_montgomery(T, montgomery_encode_limbs(T, bigint_to_limbs(canonical)))
end

canonical_bigint(x::T) where {T<:MontgomeryFiniteFieldElement} =
    limbs_to_bigint(montgomery_decode_limbs(T, montgomery_storage(x)))

@inline canonical_limbs(x::T) where {T<:MontgomeryFiniteFieldElement} =
    montgomery_decode_limbs(T, montgomery_storage(x))

Base.iszero(x::T) where {T<:MontgomeryFiniteFieldElement} = montgomery_storage(x) == zero_limbs()
Base.isone(x::T) where {T<:MontgomeryFiniteFieldElement} = montgomery_storage(x) == montgomery_one_limbs(T)

field_equal(x::T, y::T) where {T<:MontgomeryFiniteFieldElement} = montgomery_storage(x) == montgomery_storage(y)

function field_add(x::T, y::T) where {T<:MontgomeryFiniteFieldElement}
    return construct_montgomery(T, add_mod(montgomery_storage(x), montgomery_storage(y), montgomery_modulus_limbs(T)))
end

function field_sub(x::T, y::T) where {T<:MontgomeryFiniteFieldElement}
    return construct_montgomery(T, sub_mod(montgomery_storage(x), montgomery_storage(y), montgomery_modulus_limbs(T)))
end

function field_neg(x::T) where {T<:MontgomeryFiniteFieldElement}
    iszero(x) && return x
    return construct_montgomery(T, first(sub_limbs(montgomery_modulus_limbs(T), montgomery_storage(x))))
end

function field_mul(x::T, y::T) where {T<:MontgomeryFiniteFieldElement}
    return construct_montgomery(T, montgomery_mul_limbs(T, montgomery_storage(x), montgomery_storage(y)))
end

field_mul_int(x::T, n::Integer) where {T<:MontgomeryFiniteFieldElement} = field_mul(x, T(n))

function field_pow_nonnegative(x::T, n::Integer) where {T<:MontgomeryFiniteFieldElement}
    result = one(T)
    base = x
    exponent = n
    while exponent > 0
        if isodd(exponent)
            result = result * base
        end
        exponent >>= 1
        exponent == 0 && break
        base = base * base
    end
    return result
end

function field_inv(x::T) where {T<:MontgomeryFiniteFieldElement}
    iszero(x) && throw(DivideError())
    return construct_normalized(T, invmod(canonical_bigint(x), prime(T)))
end

function field_pow(x::T, n::Integer) where {T<:MontgomeryFiniteFieldElement}
    if n == 0
        return one(T)
    elseif n < 0
        return field_pow_nonnegative(inv(x), -n)
    end
    return field_pow_nonnegative(x, n)
end

function Base.getproperty(x::T, name::Symbol) where {T<:MontgomeryFiniteFieldElement}
    if name === :value
        return canonical_bigint(x)
    end
    return getfield(x, name)
end

function Base.propertynames(::T, private::Bool=false) where {T<:MontgomeryFiniteFieldElement}
    private && return (:limbs, :value)
    return (:value,)
end

# ===========================================
# BN254 Field Implementation
# ===========================================

const BN254_FQ_PRIME = parse(BigInt, "21888242871839275222246405745257275088696311157297823662689037894645226208583")
const BN254_FQ_MODULUS_LIMBS, BN254_FQ_ONE_LIMBS, BN254_FQ_R2_LIMBS, BN254_FQ_N0INV =
    montgomery_constants(BN254_FQ_PRIME)

"""
    BN254Fq

Element of the BN254 base field backed by a 4-limb Montgomery representation.
"""
struct BN254Fq <: MontgomeryFiniteFieldElement
    limbs::NTuple{4,UInt64}
end

# BN254 prime
prime(::Type{BN254Fq}) = BN254_FQ_PRIME
field_name(::Type{BN254Fq}) = "BN254"
montgomery_modulus_limbs(::Type{BN254Fq}) = BN254_FQ_MODULUS_LIMBS
montgomery_r2_limbs(::Type{BN254Fq}) = BN254_FQ_R2_LIMBS
montgomery_one_limbs(::Type{BN254Fq}) = BN254_FQ_ONE_LIMBS
montgomery_n0inv(::Type{BN254Fq}) = BN254_FQ_N0INV

# Constructors
BN254Fq(value::BigInt, normalized::Bool=false) = construct_normalized(BN254Fq, normalize_canonical(BN254Fq, value, normalized))
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

const BN254_FR_PRIME = parse(BigInt, "21888242871839275222246405745257275088548364400416034343698204186575808495617")
const BN254_FR_MODULUS_LIMBS, BN254_FR_ONE_LIMBS, BN254_FR_R2_LIMBS, BN254_FR_N0INV =
    montgomery_constants(BN254_FR_PRIME)

"""
    BN254Fr

Element of the BN254 scalar field (Fr) used for curve scalars.
"""
struct BN254Fr <: MontgomeryFiniteFieldElement
    limbs::NTuple{4,UInt64}
end

prime(::Type{BN254Fr}) = BN254_FR_PRIME
field_name(::Type{BN254Fr}) = "BN254Fr"
montgomery_modulus_limbs(::Type{BN254Fr}) = BN254_FR_MODULUS_LIMBS
montgomery_r2_limbs(::Type{BN254Fr}) = BN254_FR_R2_LIMBS
montgomery_one_limbs(::Type{BN254Fr}) = BN254_FR_ONE_LIMBS
montgomery_n0inv(::Type{BN254Fr}) = BN254_FR_N0INV

BN254Fr(value::BigInt, normalized::Bool=false) = construct_normalized(BN254Fr, normalize_canonical(BN254Fr, value, normalized))
BN254Fr(x::Integer) = BN254Fr(BigInt(x))
bn254_fr(x) = BN254Fr(x)

export BN254Fr, bn254_fr
