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

_scalar_mul_window(::Type{<:GroupElem}, ::Int) = 0

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

    w = _scalar_mul_window(typeof(P), _bit_length(k))
    if w > 0
        return scalar_mul_wnaf(P, k, w)
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

# Generic multi-scalar multiplication

const MSM_PIPPENGER_THRESHOLD = 8

"""
    multi_scalar_mul(points::Vector{GroupElem{C}}, scalars::Vector{<:Integer}) where C

Compute `∑ᵢ scalars[i] * points[i]` using a variable-base MSM backend.

Small inputs use a simple Straus-style bit scan, while larger inputs use a
Pippenger-style bucket method.
"""
function multi_scalar_mul(points::Vector{<:GroupElem{C}}, scalars::Vector{S}) where {C,S<:Integer}
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

    max_bits = maximum(_bit_length, scalars)

    if max_bits == 0
        return zero(points[1])
    end

    if length(points) < MSM_PIPPENGER_THRESHOLD
        return _multi_scalar_mul_straus(points, scalars, max_bits)
    end

    return _multi_scalar_mul_pippenger(points, scalars, max_bits)
end

function _multi_scalar_mul_straus(points::Vector{<:GroupElem{C}}, scalars::Vector{S}, max_bits::Int) where {C,S<:Integer}
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

function _multi_scalar_mul_pippenger(points::Vector{G}, scalars::Vector{S}, max_bits::Int) where {C,G<:GroupElem{C},S<:Integer}
    return _multi_scalar_mul_pippenger(points, scalars, max_bits, _pippenger_window(G, length(points)))
end

function _multi_scalar_mul_pippenger(points::Vector{G}, scalars::Vector{S}, max_bits::Int, window::Int) where {C,G<:GroupElem{C},S<:Integer}
    normalized_points = Vector{G}(undef, length(points))
    normalized_scalars = Vector{S}(undef, length(scalars))
    nonzero_count = 0

    @inbounds for i in eachindex(points, scalars)
        scalar = scalars[i]
        iszero(scalar) && continue
        nonzero_count += 1
        if scalar < 0
            normalized_points[nonzero_count] = -points[i]
            normalized_scalars[nonzero_count] = abs(scalar)
        else
            normalized_points[nonzero_count] = points[i]
            normalized_scalars[nonzero_count] = scalar
        end
    end

    if nonzero_count == 0
        return zero(points[1])
    elseif nonzero_count == 1
        return scalar_mul(normalized_points[1], normalized_scalars[1])
    elseif nonzero_count < MSM_PIPPENGER_THRESHOLD
        resize!(normalized_points, nonzero_count)
        resize!(normalized_scalars, nonzero_count)
        return _multi_scalar_mul_straus(normalized_points, normalized_scalars, max_bits)
    end

    resize!(normalized_points, nonzero_count)
    resize!(normalized_scalars, nonzero_count)

    bucket_count = (1 << window) - 1
    num_windows = cld(max_bits, window)
    zero_point = zero(points[1])
    buckets = fill(zero_point, bucket_count)
    mask = (one(S) << window) - one(S)
    result = zero_point

    for window_index in (num_windows-1):-1:0
        for _ in 1:window
            result = result + result
        end

        fill!(buckets, zero_point)
        shift = window_index * window

        @inbounds for i in eachindex(normalized_points, normalized_scalars)
            digit = _window_digit(normalized_scalars[i], shift, mask)
            if digit != 0
                buckets[digit] = buckets[digit] + normalized_points[i]
            end
        end

        # Reconstruct ∑ bucket[i] * i from the descending running sums.
        running_sum = zero_point
        @inbounds for bucket_index in bucket_count:-1:1
            running_sum = running_sum + buckets[bucket_index]
            result = result + running_sum
        end
    end

    return result
end

"""
    multi_scalar_mul(points::Vector{G}, scalars::Vector{BN254Fr}) where {C,G<:GroupElem{C}}

Compute `∑ᵢ scalars[i] * points[i]` for BN254 scalar-field elements without
round-tripping the scalars through `BigInt`.

The dispatcher decodes each scalar to canonical limbs once, then reuses the
same Straus/Pippenger split as the generic integer MSM path.
"""
function multi_scalar_mul(points::Vector{G}, scalars::Vector{BN254Fr}) where {C,G<:GroupElem{C}}
    if length(points) != length(scalars)
        throw(ArgumentError("Points and scalars must have the same length"))
    end

    isempty(points) && throw(ArgumentError("Cannot compute multi-scalar multiplication of empty vectors"))

    if length(points) == 1
        return scalar_mul(points[1], scalars[1])
    end

    scalar_limbs = Vector{NTuple{4,UInt64}}(undef, length(scalars))
    max_bits = 0
    @inbounds for i in eachindex(scalars)
        limbs = canonical_limbs(scalars[i])
        scalar_limbs[i] = limbs
        max_bits = max(max_bits, _limbs_bit_length(limbs))
    end

    if max_bits == 0
        return zero(points[1])
    end

    if length(points) < MSM_PIPPENGER_THRESHOLD
        return _multi_scalar_mul_straus(points, scalar_limbs, max_bits)
    end

    return _multi_scalar_mul_pippenger(points, scalar_limbs, max_bits, _pippenger_window(G, length(points)))
end

function _multi_scalar_mul_straus(points::Vector{<:GroupElem{C}}, scalars::Vector{NTuple{4,UInt64}}, max_bits::Int) where C
    if max_bits == 0
        return zero(points[1])
    end

    result = zero(points[1])

    for bit_pos in (max_bits - 1):-1:0
        result = result + result

        @inbounds for i in eachindex(points, scalars)
            if _limbs_testbit(scalars[i], bit_pos)
                result = result + points[i]
            end
        end
    end

    return result
end

function _multi_scalar_mul_pippenger(points::Vector{G}, scalars::Vector{NTuple{4,UInt64}}, max_bits::Int, window::Int) where {C,G<:GroupElem{C}}
    normalized_points = Vector{G}(undef, length(points))
    normalized_scalars = Vector{NTuple{4,UInt64}}(undef, length(scalars))
    nonzero_count = 0

    @inbounds for i in eachindex(points, scalars)
        limbs = scalars[i]
        _limbs_iszero(limbs) && continue
        nonzero_count += 1
        normalized_points[nonzero_count] = points[i]
        normalized_scalars[nonzero_count] = limbs
    end

    if nonzero_count == 0
        return zero(points[1])
    elseif nonzero_count == 1
        return scalar_mul(normalized_points[1], _bn254fr_from_canonical_limbs(normalized_scalars[1]))
    elseif nonzero_count < MSM_PIPPENGER_THRESHOLD
        resize!(normalized_points, nonzero_count)
        resize!(normalized_scalars, nonzero_count)
        return _multi_scalar_mul_straus(normalized_points, normalized_scalars, max_bits)
    end

    resize!(normalized_points, nonzero_count)
    resize!(normalized_scalars, nonzero_count)

    bucket_count = (1 << window) - 1
    num_windows = cld(max_bits, window)
    zero_point = zero(points[1])
    buckets = fill(zero_point, bucket_count)
    result = zero_point

    for window_index in (num_windows - 1):-1:0
        for _ in 1:window
            result = result + result
        end

        fill!(buckets, zero_point)
        shift = window_index * window

        @inbounds for i in eachindex(normalized_points, normalized_scalars)
            digit = _limbs_window_digit(normalized_scalars[i], shift, window)
            if digit != 0
                buckets[digit] = buckets[digit] + normalized_points[i]
            end
        end

        running_sum = zero_point
        @inbounds for bucket_index in bucket_count:-1:1
            running_sum = running_sum + buckets[bucket_index]
            result = result + running_sum
        end
    end

    return result
end

@inline _pippenger_window(::Type{<:GroupElem}, size::Int) = _default_pippenger_window(size)

@inline function _default_pippenger_window(size::Int)
    if size < 32
        return 3
    end

    return ((ndigits(size, base=2) - 1) * 69) ÷ 100 + 2
end

@inline _window_digit(s::T, shift::Int, mask::T) where {T<:Integer} = Int((s >> shift) & mask)

@inline function _bit_length(s::Integer)
    s == 0 && return 0
    return ndigits(abs(s), base=2)
end

@inline _limbs_iszero(limbs::NTuple{4,UInt64}) = limbs == zero_limbs()
@inline _limbs_isone(limbs::NTuple{4,UInt64}) = limbs == (UInt64(1), UInt64(0), UInt64(0), UInt64(0))

@inline function _limbs_bit_length(limbs::NTuple{4,UInt64})
    @inbounds for i in MONT_LIMB_COUNT:-1:1
        limb = limbs[i]
        limb == 0 && continue
        return (i - 1) * MONT_WORD_BITS + (MONT_WORD_BITS - leading_zeros(limb))
    end
    return 0
end

@inline _bit_length(s::BN254Fr) = _limbs_bit_length(canonical_limbs(s))

@inline function _limbs_testbit(limbs::NTuple{4,UInt64}, bit::Int)
    bit < 0 && return false
    limb_index = (bit >>> 6) + 1
    limb_index > MONT_LIMB_COUNT && return false
    offset = bit & 63
    return ((limbs[limb_index] >> offset) & UInt64(1)) == UInt64(1)
end

@inline function _limbs_window_digit(limbs::NTuple{4,UInt64}, shift::Int, width::Int)
    digit = 0
    @inbounds for bit_offset in 0:(width - 1)
        bit = shift + bit_offset
        limb_index = (bit >>> 6) + 1
        limb_index > MONT_LIMB_COUNT && break
        offset = bit & 63
        digit |= Int((limbs[limb_index] >> offset) & UInt64(1)) << bit_offset
    end
    return digit
end

@inline function _limbs_sub_small(limbs::NTuple{4,UInt64}, small::UInt64)
    out = MVector{4,UInt64}(undef)
    @inbounds for i in 1:4
        out[i] = limbs[i]
    end
    borrow = small
    @inbounds for i in 1:MONT_LIMB_COUNT
        borrow == 0 && break
        prev = out[i]
        out[i] = prev - borrow
        borrow = prev < borrow ? UInt64(1) : UInt64(0)
    end
    borrow == 0 || throw(ArgumentError("underflow while subtracting small digit from scalar limbs"))
    return (out[1], out[2], out[3], out[4])
end

@inline function _limbs_add_small(limbs::NTuple{4,UInt64}, small::UInt64)
    out = MVector{4,UInt64}(undef)
    @inbounds for i in 1:4
        out[i] = limbs[i]
    end
    carry = small
    @inbounds for i in 1:MONT_LIMB_COUNT
        carry == 0 && break
        prev = out[i]
        out[i] = prev + carry
        carry = out[i] < prev ? UInt64(1) : UInt64(0)
    end
    carry == 0 || throw(ArgumentError("overflow while adding small digit to scalar limbs"))
    return (out[1], out[2], out[3], out[4])
end

@inline function _limbs_shift_right_one(limbs::NTuple{4,UInt64})
    out = MVector{4,UInt64}(undef)
    carry = UInt64(0)
    @inbounds for i in MONT_LIMB_COUNT:-1:1
        limb = limbs[i]
        out[i] = (limb >> 1) | (carry << 63)
        carry = limb & UInt64(1)
    end
    return (out[1], out[2], out[3], out[4])
end

@inline _bn254fr_from_canonical_limbs(limbs::NTuple{4,UInt64}) =
    construct_montgomery(BN254Fr, montgomery_encode_limbs(BN254Fr, limbs))

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
    wnaf_encode(k::BN254Fr, w::Int=4)

Encode a BN254 scalar-field element in windowed non-adjacent form (w-NAF)
without converting through `BigInt`.
"""
function wnaf_encode(k::BN254Fr, w::Int=4)
    if w < 2
        throw(ArgumentError("Window size must be at least 2"))
    end

    limbs = canonical_limbs(k)
    _limbs_iszero(limbs) && return [0]

    naf = Int[]
    width_mask = UInt64((1 << w) - 1)

    while !_limbs_iszero(limbs)
        if isodd(limbs[1])
            digit = Int(limbs[1] & width_mask)
            if digit >= (1 << (w - 1))
                digit -= (1 << w)
            end
            push!(naf, digit)
            if digit >= 0
                limbs = _limbs_sub_small(limbs, UInt64(digit))
            else
                limbs = _limbs_add_small(limbs, UInt64(-digit))
            end
        else
            push!(naf, 0)
        end

        limbs = _limbs_shift_right_one(limbs)
    end

    return naf
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

"""
    scalar_mul_wnaf(P::GroupElem{C}, k::BN254Fr, w::Int=4) where C

Multiply `P` by a BN254 scalar-field element using the limb-native w-NAF path.

This mirrors the integer API while keeping the scalar in canonical limb form
throughout recoding and evaluation.
"""
function scalar_mul_wnaf(P::GroupElem{C}, k::BN254Fr, w::Int=4) where C
    if iszero(k)
        return zero(P)
    elseif isone(k)
        return P
    end

    max_odd = (1 << (w - 1)) - 1
    precomputed = Vector{typeof(P)}(undef, max_odd)
    precomputed[1] = P

    if max_odd > 1
        P2 = P + P
        for i in 2:max_odd
            precomputed[i] = precomputed[i - 1] + P2
        end
    end

    naf = wnaf_encode(k, w)
    result = zero(P)

    for i in length(naf):-1:1
        result = result + result

        digit = naf[i]
        if digit > 0
            result = result + precomputed[div(digit + 1, 2)]
        elseif digit < 0
            result = result + (-precomputed[div(-digit + 1, 2)])
        end
    end

    return result
end

"""
    scalar_mul(P::GroupElem{C}, k::BN254Fr) where C

Multiply `P` by a BN254 scalar-field element using the tuned BN254 scalar path.

The dispatcher prefers w-NAF when a curve-specific window is available and
otherwise falls back to a limb-native double-and-add scan.
"""
function scalar_mul(P::GroupElem{C}, k::BN254Fr) where C
    if iszero(k)
        return zero(P)
    elseif isone(k)
        return P
    end

    bits = _bit_length(k)
    w = _scalar_mul_window(typeof(P), bits)
    if w > 0
        return scalar_mul_wnaf(P, k, w)
    end

    limbs = canonical_limbs(k)
    result = zero(P)
    for bit_pos in (bits - 1):-1:0
        result = result + result
        if _limbs_testbit(limbs, bit_pos)
            result = result + P
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

function mul_fixed(table::FixedBaseTable{G}, k::BN254Fr) where {G<:GroupElem}
    if iszero(k)
        return zero(table.precomp[1])
    elseif isone(k)
        return table.precomp[1]
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

function batch_mul(table::FixedBaseTable{G}, scalars::Vector{BN254Fr}) where {G<:GroupElem}
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
