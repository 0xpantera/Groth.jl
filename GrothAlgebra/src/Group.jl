# Abstract group interface and generic algorithms.
#
# Provides curve-agnostic group operations and generic algorithms that work
# with any concrete group element implementation.

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

Compute `k * P` with the binary double-and-add algorithm.

This generic fallback works for any `GroupElem{C}` that implements addition
and doubling; concrete curves can override it for performance.
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

Compute `k * P` for `UInt256` scalars using a bit-scanning loop.
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

Multiply `P` by the integer scalar `k`.
"""
Base.:*(k::Integer, P::GroupElem{C}) where C = scalar_mul(P, k)

"""
    Base.:*(P::GroupElem{C}, k::Integer) where C

Multiply `P` by the integer scalar `k`.
"""
Base.:*(P::GroupElem{C}, k::Integer) where C = scalar_mul(P, k)

# Generic multi-scalar multiplication (Straus algorithm)

"""
    multi_scalar_mul(points::Vector{GroupElem{C}}, scalars::Vector{<:Integer}) where C

Compute `∑ᵢ scalars[i] * points[i]` using the Straus multi-scalar algorithm.

This shared doublings path is typically faster than independent scalar
multiplications when handling several point–scalar pairs.
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

Encode `k` in windowed non-adjacent form (w-NAF).

The returned digits are zeros or odd integers in
`[-(2^(w-1)-1), 2^(w-1)-1]`, enabling efficient use of precomputed odd
multiples.
"""
function wnaf_encode(k::Integer, w::Int=4)
    if w < 2
        throw(ArgumentError("Window size must be at least 2"))
    end

    if k == 0
        return [0]
    end

    naf = Int[]
    temp_k = BigInt(abs(k))

    while temp_k > 0
        if isodd(temp_k)  # k is odd
            # Extract the bottom w bits
            mask = BigInt((1 << w) - 1)
            digit_big = temp_k & mask
            digit = Int(digit_big)  # safe: digit < 2^w

            # If digit is too large, make it negative
            if digit >= (1 << (w-1))
                digit -= (1 << w)
            end

            push!(naf, digit)
            temp_k -= BigInt(digit)
        else
            push!(naf, 0)
        end

        temp_k >>= 1
    end

    return k >= 0 ? naf : -naf
end

"""
    scalar_mul_wnaf(P::GroupElem{C}, k::Integer, w::Int=4) where C

Multiply `P` by `k` using w-NAF digits with window `w`.

This path reuses precomputed odd multiples and often beats binary
double-and-add for large scalars.
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

# Fixed-base precomputation (w-NAF table)

"""
    FixedBaseTable{G}

Store precomputed odd multiples for a fixed base point, enabling fast w-NAF
scalar multiplications against the same base.
"""
struct FixedBaseTable{G<:GroupElem}
    window::Int
    precomp::Vector{G}  # [1P, 3P, 5P, ..., (2^(w-1)-1)P]
end

"""
    build_fixed_table(base::G; window::Int=5) where {G<:GroupElem}

Build a fixed-base table of odd multiples for `base` using the chosen window.
"""
function build_fixed_table(base::G; window::Int=5) where {G<:GroupElem}
    if window < 2
        throw(ArgumentError("Window size must be at least 2"))
    end
    max_odd = (1 << (window - 1)) - 1
    precomp = Vector{G}(undef, max_odd)
    precomp[1] = base
    if max_odd > 1
        twoP = base + base
        for i in 2:max_odd
            precomp[i] = precomp[i-1] + twoP
        end
    end
    return FixedBaseTable{G}(window, precomp)
end

"""
    mul_fixed(table::FixedBaseTable{G}, k::Integer) where {G<:GroupElem}

Multiply the table’s base by `k` using its precomputed odd multiples.
"""
function mul_fixed(table::FixedBaseTable{G}, k::Integer) where {G<:GroupElem}
    # Handle trivial cases
    if k == 0
        return zero(table.precomp[1])
    elseif k == 1
        return table.precomp[1]
    elseif k == -1
        return -table.precomp[1]
    end

    naf = wnaf_encode(k, table.window)
    acc = zero(table.precomp[1])
    for i in length(naf):-1:1
        acc = acc + acc
        d = naf[i]
        if d > 0
            acc = acc + table.precomp[div(d + 1, 2)]
        elseif d < 0
            acc = acc + (-table.precomp[div(-d + 1, 2)])
        end
    end
    return acc
end

"""
    batch_mul(table::FixedBaseTable{G}, scalars::Vector{<:Integer}) where {G<:GroupElem}

Compute `kᵢ · base` for each scalar in `scalars` via a shared fixed-base table.
"""
function batch_mul(table::FixedBaseTable{G}, scalars::Vector{<:Integer}) where {G<:GroupElem}
    return [mul_fixed(table, ki) for ki in scalars]
end

# Helper functions

"""
    is_on_curve(P::GroupElem{C}) where C

Return whether `P` satisfies the curve equation.
Concrete curve types should override this fallback.
"""
function is_on_curve(P::GroupElem{C}) where C
    # Generic fallback - concrete implementations should override
    return true
end

"""
    double(P::GroupElem{C}) where C

Return `2P` by adding `P` to itself.
"""
@inline double(P::GroupElem{C}) where C = P + P

"""
    triple(P::GroupElem{C}) where C

Return `3P` via repeated addition.
"""
@inline triple(P::GroupElem{C}) where C = P + P + P

# Utility functions for working with group elements

"""
    Base.isfinite(P::GroupElem{C}) where C

Return `true` if `P` is finite (not the point at infinity).
"""
@inline Base.isfinite(P::GroupElem{C}) where C = !iszero(P)

"""
    order(P::GroupElem{C}) where C

Return the (potentially slow to compute) order of `P`, the least positive `k`
such that `kP = 0`.
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
