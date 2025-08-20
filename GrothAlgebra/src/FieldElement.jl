"""
Field element arithmetic for prime fields.

This module provides type-stable field element operations using compile-time
prime specification via Val{P} parameters.
"""

using BitIntegers

# Error types
struct FieldError <: Exception
    msg::String
end

struct DifferentFieldsError <: Exception
    msg::String
end

# secp256k1 prime: p = 2^256 - 2^32 - 2^9 - 2^8 - 2^7 - 2^6 - 2^4 - 1
const SECP256K1_PRIME = BigInt("0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F")

"""
    FieldElem{P}

A prime field element where P is a Val{prime} compile-time parameter.
The underlying value is stored as a UInt256 and is always normalized to [0, P-1].
"""
struct FieldElem{P}
    value::UInt256

    # Inner constructor that assumes value is already normalized
    function FieldElem{P}(value::UInt256, ::Val{:normalized}) where P
        prime = _get_prime(Val{P}())
        if value >= prime
            throw(FieldError("Value $value must be less than prime $prime"))
        end
        new{P}(value)
    end
end

# Helper function to extract prime from Val{P}
@inline _get_prime(::Val{P}) where P = P

# Type alias for secp256k1 field elements
const Secp256k1Field = FieldElem{SECP256K1_PRIME}

"""
    FieldElem{P}(value::Integer) where P

Create a field element from an integer, normalizing it to [0, P-1].
"""
function FieldElem{P}(value::Integer) where P
    prime = _get_prime(Val{P}())
    normalized = mod(value, prime)
    FieldElem{P}(UInt256(normalized))
end

"""
    FieldElem{P}(value::UInt256) where P

Create a field element from a UInt256, normalizing it to [0, P-1].
"""
function FieldElem{P}(value::UInt256) where P
    prime = _get_prime(Val{P}())
    normalized = value % prime
    # Use the inner constructor with normalized flag
    FieldElem{P}(normalized, Val(:normalized))
end

# Convenience constructors
"""
    secp256k1_field(value::Integer)

Create a secp256k1 field element from an integer.
"""
secp256k1_field(value::Integer) = Secp256k1Field(value)

# Basic properties
"""
    characteristic(::Type{FieldElem{P}}) where P

Return the characteristic (prime) of the field.
"""
@inline characteristic(::Type{FieldElem{P}}) where P = _get_prime(Val{P}())

"""
    characteristic(x::FieldElem{P}) where P

Return the characteristic (prime) of the field element.
"""
@inline characteristic(x::FieldElem{P}) where P = characteristic(typeof(x))

"""
    prime_field(::Type{FieldElem{P}}) where P

Return the prime of the field (same as characteristic).
"""
@inline prime_field(::Type{FieldElem{P}}) where P = characteristic(FieldElem{P})

# Comparison operations
"""
    Base.:(==)(x::FieldElem{P}, y::FieldElem{P}) where P

Test equality of two field elements in the same field.
"""
Base.:(==)(x::FieldElem{P}, y::FieldElem{P}) where P = x.value == y.value

"""
    Base.:(==)(x::FieldElem{P}, y::FieldElem{Q}) where {P, Q}

Field elements from different fields are never equal.
"""
Base.:(==)(x::FieldElem{P}, y::FieldElem{Q}) where {P, Q} = false

"""
    Base.isequal(x::FieldElem{P}, y::FieldElem{P}) where P

Test equality of two field elements in the same field.
"""
Base.isequal(x::FieldElem{P}, y::FieldElem{P}) where P = x.value == y.value

# Zero and one elements
"""
    Base.zero(::Type{FieldElem{P}}) where P

Return the zero element of the field.
"""
@inline Base.zero(::Type{FieldElem{P}}) where P = FieldElem{P}(UInt256(0), Val(:normalized))

"""
    Base.zero(x::FieldElem{P}) where P

Return the zero element of the same field as x.
"""
@inline Base.zero(x::FieldElem{P}) where P = zero(typeof(x))

"""
    Base.one(::Type{FieldElem{P}}) where P

Return the one element of the field.
"""
@inline Base.one(::Type{FieldElem{P}}) where P = FieldElem{P}(UInt256(1), Val(:normalized))

"""
    Base.one(x::FieldElem{P}) where P

Return the one element of the same field as x.
"""
@inline Base.one(x::FieldElem{P}) where P = one(typeof(x))

# Predicates
"""
    Base.iszero(x::FieldElem{P}) where P

Test if the field element is zero.
"""
@inline Base.iszero(x::FieldElem{P}) where P = x.value == UInt256(0)

"""
    Base.isone(x::FieldElem{P}) where P

Test if the field element is one.
"""
@inline Base.isone(x::FieldElem{P}) where P = x.value == UInt256(1)

"""
    is_zero(x::FieldElem{P}) where P

Test if the field element is zero.
"""
@inline is_zero(x::FieldElem{P}) where P = iszero(x)

"""
    is_one(x::FieldElem{P}) where P

Test if the field element is one.
"""
@inline is_one(x::FieldElem{P}) where P = isone(x)

"""
    is_unity(x::FieldElem{P}) where P

Test if the field element is the multiplicative identity (one).
"""
@inline is_unity(x::FieldElem{P}) where P = isone(x)

# Arithmetic operations
"""
    Base.:+(x::FieldElem{P}, y::FieldElem{P}) where P

Add two field elements in the same field.
"""
function Base.:+(x::FieldElem{P}, y::FieldElem{P}) where P
    prime = _get_prime(Val{P}())
    result = x.value + y.value
    if result >= prime
        result -= prime
    end
    FieldElem{P}(result, Val(:normalized))
end

"""
    Base.:-(x::FieldElem{P}, y::FieldElem{P}) where P

Subtract two field elements in the same field.
"""
function Base.:-(x::FieldElem{P}, y::FieldElem{P}) where P
    prime = _get_prime(Val{P}())
    if x.value >= y.value
        result = x.value - y.value
    else
        result = prime - (y.value - x.value)
    end
    FieldElem{P}(result, Val(:normalized))
end

"""
    Base.:-(x::FieldElem{P}) where P

Negate a field element.
"""
function Base.:-(x::FieldElem{P}) where P
    if iszero(x)
        return x
    end
    prime = _get_prime(Val{P}())
    FieldElem{P}(prime - x.value, Val(:normalized))
end

"""
    Base.:*(x::FieldElem{P}, y::FieldElem{P}) where P

Multiply two field elements in the same field.
"""
function Base.:*(x::FieldElem{P}, y::FieldElem{P}) where P
    prime = _get_prime(Val{P}())
    # Use wider arithmetic to avoid overflow
    result = (UInt512(x.value) * UInt512(y.value)) % UInt512(prime)
    FieldElem{P}(UInt256(result), Val(:normalized))
end

"""
    Base.:*(x::FieldElem{P}, n::Integer) where P

Multiply a field element by an integer (scalar multiplication).
"""
function Base.:*(x::FieldElem{P}, n::Integer) where P
    prime = _get_prime(Val{P}())
    n_normalized = mod(n, prime)
    result = (UInt512(x.value) * UInt512(n_normalized)) % UInt512(prime)
    FieldElem{P}(UInt256(result), Val(:normalized))
end

"""
    Base.:*(n::Integer, x::FieldElem{P}) where P

Multiply an integer by a field element (scalar multiplication).
"""
Base.:*(n::Integer, x::FieldElem{P}) where P = x * n

"""
    Base.inv(x::FieldElem{P}) where P

Compute the multiplicative inverse of a field element using Fermat's little theorem.
For prime p, a^(-1) ≡ a^(p-2) (mod p).
This is constant-time for fixed prime P.
"""
function Base.inv(x::FieldElem{P}) where P
    if iszero(x)
        throw(DivideError())
    end
    prime = _get_prime(Val{P}())
    # Use Fermat's little theorem: a^(-1) ≡ a^(p-2) (mod p)
    exponent = prime - UInt256(2)
    result = _pow_mod(x.value, exponent, prime)
    FieldElem{P}(result, Val(:normalized))
end

"""
    Base.:/(x::FieldElem{P}, y::FieldElem{P}) where P

Divide two field elements in the same field.
"""
Base.:/(x::FieldElem{P}, y::FieldElem{P}) where P = x * inv(y)

"""
    Base.:^(x::FieldElem{P}, n::Integer) where P

Raise a field element to an integer power.
"""
function Base.:^(x::FieldElem{P}, n::Integer) where P
    if n == 0
        return one(x)
    elseif n < 0
        return inv(x)^(-n)
    end

    prime = _get_prime(Val{P}())
    # Reduce exponent modulo p-1 (Fermat's little theorem)
    n_reduced = mod(n, prime - UInt256(1))
    result = _pow_mod(x.value, UInt256(n_reduced), prime)
    FieldElem{P}(result, Val(:normalized))
end

"""
    is_odd(x::FieldElem{P}) where P

Test if the field element has an odd underlying value.
"""
@inline is_odd(x::FieldElem{P}) where P = isodd(x.value)

# Helper function for modular exponentiation using iterative square-and-multiply
function _pow_mod(base::UInt256, exp::UInt256, mod::UInt256)
    if exp == UInt256(0)
        return UInt256(1)
    end

    if mod == UInt256(1)
        return UInt256(0)
    end

    result = UInt256(1)
    base = base % mod
    exp_copy = exp

    while exp_copy > UInt256(0)
        if exp_copy & UInt256(1) == UInt256(1)
            result = (UInt512(result) * UInt512(base)) % UInt512(mod) |> UInt256
        end
        exp_copy >>= 1
        if exp_copy > UInt256(0)  # Avoid unnecessary squaring on last iteration
            base = (UInt512(base) * UInt512(base)) % UInt512(mod) |> UInt256
        end
    end

    return result
end

# Display
function Base.show(io::IO, x::FieldElem{P}) where P
    prime = _get_prime(Val{P}())
    if prime == SECP256K1_PRIME
        print(io, "0x$(string(x.value, base=16, pad=2)) mod secp256k1")
    else
        # For other primes, show a shorter format
        prime_str = string(prime)
        if length(prime_str) > 20
            prime_short = "$(prime_str[1:8])...$(prime_str[end-7:end])"
        else
            prime_short = prime_str
        end
        print(io, "0x$(string(x.value, base=16, pad=2)) mod $(prime_short)")
    end
end

# Conversion to/from integers
"""
    Base.convert(::Type{Integer}, x::FieldElem{P}) where P

Convert a field element to its underlying integer value.
"""
Base.convert(::Type{Integer}, x::FieldElem{P}) where P = x.value

"""
    Base.convert(::Type{UInt256}, x::FieldElem{P}) where P

Convert a field element to its underlying UInt256 value.
"""
Base.convert(::Type{UInt256}, x::FieldElem{P}) where P = x.value

# Promotion rules for mixed arithmetic with integers
"""
    Base.promote_rule(::Type{FieldElem{P}}, ::Type{<:Integer}) where P

Promote integer to field element for mixed arithmetic operations.
"""
Base.promote_rule(::Type{FieldElem{P}}, ::Type{<:Integer}) where P = FieldElem{P}

"""
    Base.convert(::Type{FieldElem{P}}, n::Integer) where P

Convert an integer to a field element.
"""
Base.convert(::Type{FieldElem{P}}, n::Integer) where P = FieldElem{P}(n)

# Broadcasting support
"""
    Base.Broadcast.broadcastable(x::FieldElem{P}) where P

Make field elements broadcastable for vectorized operations.
"""
Base.Broadcast.broadcastable(x::FieldElem{P}) where P = Ref(x)

# Promoted arithmetic operations for mixed FieldElem + Integer arithmetic
"""
    Base.:+(x::FieldElem{P}, n::Integer) where P

Add a field element and an integer using promotion.
"""
function Base.:+(x::FieldElem{P}, n::Integer) where P
    y = convert(FieldElem{P}, n)
    return x + y
end

"""
    Base.:+(n::Integer, x::FieldElem{P}) where P

Add an integer and a field element using promotion.
"""
Base.:+(n::Integer, x::FieldElem{P}) where P = x + n

"""
    Base.:-(x::FieldElem{P}, n::Integer) where P

Subtract an integer from a field element using promotion.
"""
function Base.:-(x::FieldElem{P}, n::Integer) where P
    y = convert(FieldElem{P}, n)
    return x - y
end

"""
    Base.:-(n::Integer, x::FieldElem{P}) where P

Subtract a field element from an integer using promotion.
"""
function Base.:-(n::Integer, x::FieldElem{P}) where P
    y = convert(FieldElem{P}, n)
    return y - x
end

"""
    Base.:/(x::FieldElem{P}, n::Integer) where P

Divide a field element by an integer using promotion.
"""
function Base.:/(x::FieldElem{P}, n::Integer) where P
    y = convert(FieldElem{P}, n)
    return x / y
end
