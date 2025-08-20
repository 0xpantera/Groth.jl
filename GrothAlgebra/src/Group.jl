"""
Abstract group interface and generic algorithms.

This module provides curve-agnostic group operations and generic algorithms
that work with any concrete group element implementation.
"""

"""
    AbstractCurve

Abstract type for curve parameter structs. Each concrete curve should
define a struct that inherits from this type to tag group elements.
"""
abstract type AbstractCurve end

"""
    GroupElem{C<:AbstractCurve}

Abstract type for group elements over curve C. Concrete implementations
should inherit from this type and implement the required operations.

Required operations for concrete implementations:
- `+(::GroupElem{C}, ::GroupElem{C})` - Group addition
- `-(::GroupElem{C})` - Group negation
- `zero(::Type{GroupElem{C}})` - Group identity element
- `iszero(::GroupElem{C})` - Test for identity element
"""
abstract type GroupElem{C<:AbstractCurve} end

# Generic group operations

"""
    Base.zero(x::GroupElem{C}) where C

Return the identity element of the same group as x.
"""
@inline Base.zero(x::GroupElem{C}) where C = zero(typeof(x))

"""
    Base.iszero(x::GroupElem{C}) where C

Test if the group element is the identity element.
Concrete implementations MUST override this method to avoid infinite recursion.
"""
Base.iszero(x::GroupElem{C}) where C = throw(MethodError(iszero, (x,), "Concrete GroupElem types must implement iszero to avoid infinite recursion with =="))

"""
    Base.:-(x::GroupElem{C}, y::GroupElem{C}) where C

Subtract two group elements: x - y = x + (-y).
"""
Base.:-(x::GroupElem{C}, y::GroupElem{C}) where C = x + (-y)

"""
    group_identity(x::GroupElem{C}) where C

Return the identity element of the same group as x.
"""
@inline group_identity(x::GroupElem{C}) where C = zero(x)



"""
    Base.inv(x::GroupElem{C}) where C

Return the inverse of a group element (same as negation for additive groups).
"""
@inline Base.inv(x::GroupElem{C}) where C = -x

# Generic scalar multiplication using double-and-add

"""
    scalar_mul(P::GroupElem{C}, k::Integer) where C

Generic scalar multiplication using binary double-and-add algorithm.
Computes k * P where k is an integer and P is a group element.

This is a generic implementation that works for any GroupElem{C} that
implements addition and doubling. Concrete implementations may override
this for better performance.
"""
function scalar_mul(P::GroupElem{C}, k::Integer) where C
    if k == 0
        return zero(P)
    elseif k == 1
        return P
    elseif k < 0
        return scalar_mul(-P, -k)
    end

    # Double-and-add algorithm
    Q = zero(P)
    bits = reverse(digits(k, base=2))

    for bit in bits
        Q = Q + Q  # Double (dispatches to curve-specific addition)
        if bit == 1
            Q = Q + P  # Add (dispatches to curve-specific addition)
        end
    end

    return Q
end

"""
    scalar_mul(P::GroupElem{C}, k::UInt256) where C

Optimized scalar multiplication for UInt256 scalars.
"""
function scalar_mul(P::GroupElem{C}, k::UInt256) where C
    if k == UInt256(0)
        return zero(P)
    elseif k == UInt256(1)
        return P
    end

    # Double-and-add with bit scanning
    Q = zero(P)
    temp_P = P
    temp_k = k

    while temp_k > UInt256(0)
        if temp_k & UInt256(1) == UInt256(1)
            Q = Q + temp_P
        end
        temp_P = temp_P + temp_P  # Double
        temp_k >>= 1
    end

    return Q
end

"""
    Base.:*(k::Integer, P::GroupElem{C}) where C

Scalar multiplication: k * P.
"""
Base.:*(k::Integer, P::GroupElem{C}) where C = scalar_mul(P, k)

"""
    Base.:*(P::GroupElem{C}, k::Integer) where C

Scalar multiplication: P * k.
"""
Base.:*(P::GroupElem{C}, k::Integer) where C = scalar_mul(P, k)

# Generic multi-scalar multiplication (Straus algorithm)

"""
    multi_scalar_mul(points::Vector{GroupElem{C}}, scalars::Vector{<:Integer}) where C

Generic multi-scalar multiplication using Straus algorithm.
Computes ∑ᵢ scalars[i] * points[i] efficiently.

This is more efficient than computing each scalar multiplication separately
when there are multiple point-scalar pairs.
"""
function multi_scalar_mul(points::Vector{<:GroupElem{C}}, scalars::Vector{<:Integer}) where C
    if length(points) != length(scalars)
        throw(ArgumentError("Points and scalars must have the same length"))
    end

    if isempty(points)
        # Need a way to get the zero element - use the first point if available
        throw(ArgumentError("Cannot compute multi-scalar multiplication of empty vectors"))
    end

    if length(points) == 1
        return scalar_mul(points[1], scalars[1])
    end

    # Straus algorithm (simultaneous multiple point scalar multiplication)
    max_bits = maximum(s -> s == 0 ? 0 : ceil(Int, log2(abs(s) + 1)), scalars)

    if max_bits == 0
        return zero(points[1])
    end

    result = zero(points[1])

    # Process bits from most significant to least significant
    for bit_pos in (max_bits-1):-1:0
        result = result + result  # Double

        # Add points whose scalar has bit set at this position
        for i in 1:length(points)
            if scalars[i] != 0 && (abs(scalars[i]) >> bit_pos) & 1 == 1
                if scalars[i] > 0
                    result = result + points[i]
                else
                    result = result + (-points[i])
                end
            end
        end
    end

    return result
end

# w-NAF (windowed Non-Adjacent Form) utilities

"""
    wnaf_encode(k::Integer, w::Int=4)

Encode an integer in w-NAF (windowed Non-Adjacent Form).
Returns a vector of signed integers where each element is either 0 or an odd
integer in the range [-(2^(w-1)-1), 2^(w-1)-1].

w-NAF reduces the number of additions in scalar multiplication by allowing
precomputed odd multiples.
"""
function wnaf_encode(k::Integer, w::Int=4)
    if w < 2
        throw(ArgumentError("Window size must be at least 2"))
    end

    if k == 0
        return [0]
    end

    naf = Int[]
    temp_k = abs(k)

    while temp_k > 0
        if temp_k & 1 == 1  # k is odd
            # Extract the bottom w bits
            mask = (1 << w) - 1
            digit = temp_k & mask

            # If digit is too large, make it negative
            if digit >= (1 << (w-1))
                digit -= (1 << w)
            end

            push!(naf, digit)
            temp_k -= digit
        else
            push!(naf, 0)
        end

        temp_k >>= 1
    end

    return k >= 0 ? naf : -naf
end

"""
    scalar_mul_wnaf(P::GroupElem{C}, k::Integer, w::Int=4) where C

Scalar multiplication using w-NAF encoding.
This can be more efficient than binary double-and-add for large scalars,
especially when precomputed odd multiples are available.
"""
function scalar_mul_wnaf(P::GroupElem{C}, k::Integer, w::Int=4) where C
    if k == 0
        return zero(P)
    elseif k == 1
        return P
    elseif k == -1
        return -P
    end

    # Precompute odd multiples: P, 3P, 5P, ..., (2^(w-1)-1)P
    max_odd = (1 << (w - 1)) - 1
    precomputed = Vector{typeof(P)}(undef, max_odd)
    precomputed[1] = P  # 1P

    if max_odd > 1
        P2 = P + P  # 2P
        for i in 2:max_odd
            precomputed[i] = precomputed[i-1] + P2  # (2i-1)P
        end
    end

    # Encode k in w-NAF
    naf = wnaf_encode(k, w)

    # Compute scalar multiplication
    result = zero(P)

    for i in length(naf):-1:1
        result = result + result  # Double

        digit = naf[i]
        if digit > 0
            result = result + precomputed[div(digit + 1, 2)]
        elseif digit < 0
            result = result + (-precomputed[div(-digit + 1, 2)])
        end
    end

    return result
end

# Helper functions

"""
    is_on_curve(P::GroupElem{C}) where C

Test if a group element satisfies the curve equation.
This is a generic fallback that should be overridden by concrete implementations.
"""
function is_on_curve(P::GroupElem{C}) where C
    # Generic fallback - concrete implementations should override
    return true
end

"""
    double(P::GroupElem{C}) where C

Double a group element: 2P = P + P.
This is a convenience function that uses the addition operation.
"""
@inline double(P::GroupElem{C}) where C = P + P

"""
    triple(P::GroupElem{C}) where C

Triple a group element: 3P = P + P + P.
"""
@inline triple(P::GroupElem{C}) where C = P + P + P

# Utility functions for working with group elements

"""
    Base.isfinite(P::GroupElem{C}) where C

Test if a group element is finite (not the point at infinity).
For additive groups, this is equivalent to !iszero(P).
"""
@inline Base.isfinite(P::GroupElem{C}) where C = !iszero(P)

"""
    order(P::GroupElem{C}) where C

Compute the order of a group element (the smallest positive integer k such that kP = 0).
This is a generic implementation that may be slow for large orders.
"""
function order(P::GroupElem{C}) where C
    if iszero(P)
        return 1
    end

    Q = P
    k = 1

    while !iszero(Q)
        Q = Q + P
        k += 1

        # Safety check to avoid infinite loops
        if k > 1000000
            throw(OverflowError("Order computation exceeded maximum iterations"))
        end
    end

    return k
end
