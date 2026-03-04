# Polynomial arithmetic over prime fields.
#
# Provides polynomial evaluation, interpolation, and FFT-based operations for
# field elements.

"""
    Polynomial{F<:FieldElem}

A polynomial with coefficients in a field F.
Coefficients are stored in ascending order of powers: a₀ + a₁x + a₂x² + ...
"""
struct Polynomial{F<:FieldElem}
    coeffs::Vector{F}

    # Inner constructor that normalizes the polynomial by removing leading zeros
    function Polynomial{F}(coeffs::Vector{F}) where F
        # Remove leading zeros
        while length(coeffs) > 1 && iszero(coeffs[end])
            pop!(coeffs)
        end

        # Handle the zero polynomial case
        if isempty(coeffs)
            coeffs = [zero(F)]
        end

        new{F}(coeffs)
    end
end

"""
    Polynomial(coeffs::Vector{F}) where F<:FieldElem

Create a polynomial from a vector of field element coefficients.
"""
Polynomial(coeffs::Vector{F}) where F<:FieldElem = Polynomial{F}(coeffs)

"""
    Polynomial(coeffs::Vector{<:Integer}, ::Type{F}) where F<:FieldElem

Create a polynomial from integer coefficients, converting them to field elements.
"""
function Polynomial(coeffs::Vector{<:Integer}, ::Type{F}) where F<:FieldElem
    field_coeffs = [F(c) for c in coeffs]
    Polynomial{F}(field_coeffs)
end

# Basic properties

"""
    degree(p::Polynomial{F}) where F

Return the degree of the polynomial.
"""
function degree(p::Polynomial{F}) where F
    if length(p.coeffs) == 1 && iszero(p.coeffs[1])
        return -1  # Degree of zero polynomial is -1 by convention
    end
    return length(p.coeffs) - 1
end

"""
    leading_coefficient(p::Polynomial{F}) where F

Return the leading coefficient of the polynomial.
"""
function leading_coefficient(p::Polynomial{F}) where F
    return p.coeffs[end]
end

"""
    is_zero(p::Polynomial{F}) where F

Test if the polynomial is the zero polynomial.
"""
function is_zero(p::Polynomial{F}) where F
    return length(p.coeffs) == 1 && iszero(p.coeffs[1])
end

"""
    is_constant(p::Polynomial{F}) where F

Test if the polynomial is constant (degree 0).
"""
function is_constant(p::Polynomial{F}) where F
    return degree(p) <= 0
end

"""
    is_monic(p::Polynomial{F}) where F

Test if the polynomial is monic (leading coefficient is 1).
"""
function is_monic(p::Polynomial{F}) where F
    return isone(leading_coefficient(p))
end

# Zero and constant polynomials

"""
    Base.zero(::Type{Polynomial{F}}) where F

Return the zero polynomial.
"""
Base.zero(::Type{Polynomial{F}}) where F = Polynomial{F}([zero(F)])

"""
    Base.zero(p::Polynomial{F}) where F

Return the zero polynomial of the same type as p.
"""
Base.zero(p::Polynomial{F}) where F = zero(typeof(p))

"""
    Base.one(::Type{Polynomial{F}}) where F

Return the constant polynomial 1.
"""
Base.one(::Type{Polynomial{F}}) where F = Polynomial{F}([one(F)])

"""
    Base.one(p::Polynomial{F}) where F

Return the constant polynomial 1 of the same type as p.
"""
Base.one(p::Polynomial{F}) where F = one(typeof(p))

"""
    constant_polynomial(c::F) where F<:FieldElem

Create a constant polynomial with value c.
"""
constant_polynomial(c::F) where F<:FieldElem = Polynomial{F}([c])

"""
    monomial(n::Integer, ::Type{F}) where F<:FieldElem

Create the monomial x^n.
"""
function monomial(n::Integer, ::Type{F}) where F<:FieldElem
    if n < 0
        throw(ArgumentError("Monomial degree must be non-negative"))
    end
    coeffs = [zero(F) for _ in 1:(n + 1)]
    coeffs[end] = one(F)
    return Polynomial{F}(coeffs)
end

# Comparison

"""
    Base.:(==)(p::Polynomial{F}, q::Polynomial{F}) where F

Test equality of two polynomials.
"""
function Base.:(==)(p::Polynomial{F}, q::Polynomial{F}) where F
    if length(p.coeffs) != length(q.coeffs)
        return false
    end
    return p.coeffs == q.coeffs
end

"""
    Base.isequal(p::Polynomial{F}, q::Polynomial{F}) where F

Test equality of two polynomials.
"""
Base.isequal(p::Polynomial{F}, q::Polynomial{F}) where F = p == q

# Arithmetic operations

"""
    Base.:+(p::Polynomial{F}, q::Polynomial{F}) where F

Add two polynomials.
"""
function Base.:+(p::Polynomial{F}, q::Polynomial{F}) where F
    max_len = max(length(p.coeffs), length(q.coeffs))
    result = [zero(F) for _ in 1:max_len]

    for i in 1:length(p.coeffs)
        result[i] += p.coeffs[i]
    end

    for i in 1:length(q.coeffs)
        result[i] += q.coeffs[i]
    end

    return Polynomial{F}(result)
end

"""
    Base.:-(p::Polynomial{F}, q::Polynomial{F}) where F

Subtract two polynomials.
"""
function Base.:-(p::Polynomial{F}, q::Polynomial{F}) where F
    max_len = max(length(p.coeffs), length(q.coeffs))
    result = [zero(F) for _ in 1:max_len]

    for i in 1:length(p.coeffs)
        result[i] += p.coeffs[i]
    end

    for i in 1:length(q.coeffs)
        result[i] -= q.coeffs[i]
    end

    return Polynomial{F}(result)
end

"""
    Base.:-(p::Polynomial{F}) where F

Negate a polynomial.
"""
function Base.:-(p::Polynomial{F}) where F
    return Polynomial{F}([-c for c in p.coeffs])
end

"""
    Base.:*(p::Polynomial{F}, q::Polynomial{F}) where F

Multiply two polynomials.
"""
function Base.:*(p::Polynomial{F}, q::Polynomial{F}) where F
    if is_zero(p) || is_zero(q)
        return zero(Polynomial{F})
    end

    result_len = length(p.coeffs) + length(q.coeffs) - 1
    result = [zero(F) for _ in 1:result_len]

    for i in 1:length(p.coeffs)
        for j in 1:length(q.coeffs)
            result[i + j - 1] += p.coeffs[i] * q.coeffs[j]
        end
    end

    return Polynomial{F}(result)
end

"""
    Base.:*(p::Polynomial{F}, c::F) where F

Multiply a polynomial by a scalar.
"""
function Base.:*(p::Polynomial{F}, c::F) where F
    if iszero(c)
        return zero(Polynomial{F})
    end
    return Polynomial{F}([coeff * c for coeff in p.coeffs])
end

"""
    Base.:*(c::F, p::Polynomial{F}) where F

Multiply a scalar by a polynomial.
"""
Base.:*(c::F, p::Polynomial{F}) where F = p * c

"""
    Base.:^(p::Polynomial{F}, n::Integer) where F

Raise a polynomial to an integer power.
"""
function Base.:^(p::Polynomial{F}, n::Integer) where F
    if n < 0
        throw(ArgumentError("Polynomial exponent must be non-negative"))
    elseif n == 0
        return one(Polynomial{F})
    elseif n == 1
        return p
    end

    # Binary exponentiation
    result = one(Polynomial{F})
    base = p
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

# Evaluation

"""
    evaluate(p::Polynomial{F}, x::F) where F

Evaluate the polynomial at a given point using Horner's method.
"""
function evaluate(p::Polynomial{F}, x::F) where F
    if is_zero(p)
        return zero(F)
    end

    # Horner's method: p(x) = a₀ + x(a₁ + x(a₂ + ... + x(aₙ)))
    result = p.coeffs[end]
    for i in (length(p.coeffs)-1):-1:1
        result = result * x + p.coeffs[i]
    end

    return result
end

"""
    (p::Polynomial{F})(x::F) where F

Evaluate the polynomial at a given point (callable syntax).
"""
(p::Polynomial{F})(x::F) where F = evaluate(p, x)

"""
    evaluate(p::Polynomial{F}, points::Vector{F}) where F

Evaluate the polynomial at multiple points.
"""
function evaluate(p::Polynomial{F}, points::Vector{F}) where F
    return [evaluate(p, x) for x in points]
end

# Interpolation

"""
    interpolate(points::Vector{F}, values::Vector{F}) where F<:FieldElem

Construct a polynomial that passes through the given points using Lagrange interpolation.
Returns the unique polynomial of degree at most n-1 that passes through the n points.
"""
function interpolate(points::Vector{F}, values::Vector{F}) where F<:FieldElem
    if length(points) != length(values)
        throw(ArgumentError("Points and values must have the same length"))
    end

    n = length(points)
    if n == 0
        throw(ArgumentError("Need at least one point for interpolation"))
    end

    if n == 1
        return constant_polynomial(values[1])
    end

    # Lagrange interpolation
    result = zero(Polynomial{F})

    for i in 1:n
        # Construct the i-th Lagrange basis polynomial
        basis = one(Polynomial{F})
        for j in 1:n
            if i != j
                # basis *= (x - points[j]) / (points[i] - points[j])
                numerator = Polynomial{F}([-points[j], one(F)])  # x - points[j]
                denominator = points[i] - points[j]
                basis = basis * numerator * inv(denominator)
            end
        end

        result = result + values[i] * basis
    end

    return result
end

# Derivative

"""
    derivative(p::Polynomial{F}) where F

Compute the formal derivative of the polynomial.
"""
function derivative(p::Polynomial{F}) where F
    if degree(p) <= 0
        return zero(Polynomial{F})
    end

    n = length(p.coeffs)
    result_coeffs = [zero(F) for _ in 1:(n - 1)]

    for i in 2:n
        result_coeffs[i - 1] = F(i - 1) * p.coeffs[i]
    end

    return Polynomial{F}(result_coeffs)
end

# Display

function Base.show(io::IO, p::Polynomial{F}) where F
    if is_zero(p)
        print(io, "Polynomial{$(F)}(0)")
        return
    end

    # Show type and degree information
    deg = degree(p)
    field_name = string(F)
    if occursin("Secp256k1Field", field_name)
        field_name = "Secp256k1Field"
    end

    print(io, "Polynomial{$(field_name)}(degree = $(deg))")

    # For small degree polynomials, also show the expression
    if deg <= 3 && length(p.coeffs) <= 4
        print(io, ": ")
        first_term = true

        for (i, coeff) in enumerate(p.coeffs)
            if iszero(coeff)
                continue
            end

            power = i - 1

            if !first_term
                print(io, " + ")
            end
            first_term = false

            if power == 0
                print(io, "c₀")
            elseif power == 1
                if isone(coeff)
                    print(io, "x")
                else
                    print(io, "c₁⋅x")
                end
            else
                if isone(coeff)
                    print(io, "x^$power")
                else
                    print(io, "c₊⋅x^$power")
                end
            end
        end
    end
end

# -----------------------------------------------------------------------------
# FFT helpers and roots of unity
# -----------------------------------------------------------------------------

const BN254_TWO_ADICITY = 28
const BN254_PRIMITIVE_ROOT = BN254Fr(parse(BigInt, "19103219067921713944291392827692070036145651957329286315305642004821462161904"), true)

max_power_of_two(::Type{F}) where F = throw(ArgumentError("No two-adicity metadata for field $(F)"))
primitive_root_constant(::Type{F}) where F = throw(ArgumentError("No primitive root configured for field $(F)"))

max_power_of_two(::Type{BN254Fr}) = BN254_TWO_ADICITY
primitive_root_constant(::Type{BN254Fr}) = BN254_PRIMITIVE_ROOT

"""
    primitive_root_of_unity(F::Type{<:FieldElem}, log_n::Int)

Return a primitive 2^log_n-th root of unity for field type `F`.
"""
function primitive_root_of_unity(F::Type{<:FieldElem}, log_n::Int)
    log_n >= 0 || throw(ArgumentError("log_n must be non-negative"))
    max_pow = max_power_of_two(F)
    log_n <= max_pow || throw(ArgumentError("Requested 2^$log_n-th root exceeds two-adicity $max_pow for field $(F)"))
    exponent = 1 << (max_pow - log_n)
    return primitive_root_constant(F)^exponent
end

"""
    roots_of_unity(n::Integer, ::Type{F}) where F<:FieldElem

Return the list `[1, ω, ω^2, …]` of n-th roots of unity (power-of-two only).
"""
function roots_of_unity(n::Integer, ::Type{F}) where F<:FieldElem
    n > 0 || throw(ArgumentError("n must be positive"))
    ispow2(n) || throw(ArgumentError("Only power-of-two sizes are supported"))
    log_n = trailing_zeros(n)
    ω = primitive_root_of_unity(F, log_n)
    roots = Vector{F}(undef, n)
    roots[1] = one(F)
    for i in 2:n
        roots[i] = roots[i-1] * ω
    end
    return roots
end

struct EvaluationDomain{F<:FieldElem}
    size::Int
    log_size::Int
    generator::F
    generator_inv::F
    size_inv::F
    offset::F
    offset_inv::F
    offset_pow_size::F
end

@inline function _ensure_offset(::Type{F}, offset) where F<:FieldElem
    offset === nothing && return one(F)
    offset isa F && !iszero(offset) && return offset
    candidate = F(offset)
    iszero(candidate) && throw(ArgumentError("Coset offset must be non-zero"))
    return candidate
end

function EvaluationDomain(::Type{F}, size::Int; offset=nothing) where F<:FieldElem
    size > 0 || throw(ArgumentError("Domain size must be positive"))
    ispow2(size) || throw(ArgumentError("Domain size must be a power of two"))
    log_size = trailing_zeros(size)
    ω = primitive_root_of_unity(F, log_size)
    off = _ensure_offset(F, offset)
    EvaluationDomain{F}(size, log_size, ω, inv(ω), inv(F(size)), off, inv(off), off^size)
end

function get_coset(domain::EvaluationDomain{F}, offset) where F<:FieldElem
    off = _ensure_offset(F, offset)
    new_off = domain.offset * off
    EvaluationDomain{F}(
        domain.size,
        domain.log_size,
        domain.generator,
        domain.generator_inv,
        domain.size_inv,
        new_off,
        inv(new_off),
        new_off^domain.size,
    )
end

coset_offset(domain::EvaluationDomain) = domain.offset
coset_offset_inv(domain::EvaluationDomain) = domain.offset_inv
coset_offset_pow_size(domain::EvaluationDomain) = domain.offset_pow_size

function bitreverse!(values::Vector)
    n = length(values)
    n > 1 || return values
    log_n = trailing_zeros(n)
    width = UInt(sizeof(Int) * 8)
    shift = width - UInt(log_n)
    for i in 0:(n - 1)
        j = Int(Base.bitreverse(UInt(i)) >>> shift)
        if i < j
            values[i + 1], values[j + 1] = values[j + 1], values[i + 1]
        end
    end
    return values
end

function scale_powers!(values::Vector{F}, base::F) where F<:FieldElem
    isone(base) && return values
    pow = one(F)
    for i in eachindex(values)
        values[i] *= pow
        pow *= base
    end
    return values
end

function ntt!(values::Vector{F}, domain::EvaluationDomain{F}; inverse::Bool=false) where F<:FieldElem
    length(values) == domain.size || throw(ArgumentError("Vector length mismatch for NTT domain"))
    bitreverse!(values)
    log_n = domain.log_size
    root = inverse ? domain.generator_inv : domain.generator
    for s in 1:log_n
        m = 1 << s
        half = m >>> 1
        w_m = root^(1 << (log_n - s))
        k = 1
        while k <= domain.size
            w = one(F)
            for j in 0:(half - 1)
                u = values[k + j]
                t = w * values[k + j + half]
                values[k + j] = u + t
                values[k + j + half] = u - t
                w *= w_m
            end
            k += m
        end
    end
    if inverse
        inv_size = domain.size_inv
        for i in 1:domain.size
            values[i] *= inv_size
        end
    end
    return values
end

function fft(coeffs::Vector{F}, domain::EvaluationDomain{F}) where F<:FieldElem
    length(coeffs) <= domain.size || throw(ArgumentError(
        "fft: coefficient vector length $(length(coeffs)) exceeds domain size $(domain.size); " *
        "choose a larger EvaluationDomain or reduce polynomial degree",
    ))
    padded = Vector{F}(undef, domain.size)
    fill!(padded, zero(F))
    @inbounds for i in 1:length(coeffs)
        padded[i] = coeffs[i]
    end
    scale_powers!(padded, domain.offset)
    ntt!(padded, domain)
    return padded
end

function ifft(evals::Vector{F}, domain::EvaluationDomain{F}) where F<:FieldElem
    values = copy(evals)
    ntt!(values, domain; inverse=true)
    scale_powers!(values, domain.offset_inv)
    return values
end

function interpolate_fft(domain::EvaluationDomain{F}, values::Vector{F}) where F<:FieldElem
    length(values) <= domain.size || throw(ArgumentError("Too many evaluation points for domain"))
    buffer = Vector{F}(undef, domain.size)
    fill!(buffer, zero(F))
    for i in 1:length(values)
        buffer[i] = values[i]
    end
    coeffs = ifft(buffer, domain)
    return Polynomial{F}(coeffs)
end

"""
    fft_polynomial_multiply(p::Polynomial{F}, q::Polynomial{F}) where F

Multiply two polynomials using FFT when supported for the coefficient field.
Falls back to naive multiplication if FFT metadata is unavailable.
"""
function fft_polynomial_multiply(p::Polynomial{F}, q::Polynomial{F}) where F
    try
        total_len = length(p.coeffs) + length(q.coeffs) - 1
        total_len <= 0 && return zero(p)
        size = 1 << max(1, ceil(Int, log2(total_len)))
        size >= total_len || (size <<= 1)
        domain = EvaluationDomain(F, size)
        a = fft(p.coeffs, domain)
        b = fft(q.coeffs, domain)
        for i in 1:size
            a[i] *= b[i]
        end
        coeffs = ifft(a, domain)
        result_coeffs = coeffs[1:total_len]
        return Polynomial{F}(result_coeffs)
    catch err
        if err isa ArgumentError
            return p * q
        else
            rethrow(err)
        end
    end
end
