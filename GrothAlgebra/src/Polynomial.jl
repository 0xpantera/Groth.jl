"""
Polynomial arithmetic over prime fields.

This module provides polynomial operations for field elements, including
evaluation, interpolation, and FFT-based operations.
"""

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

# FFT helpers (basic implementation)

"""
    fft_polynomial_multiply(p::Polynomial{F}, q::Polynomial{F}) where F

Multiply two polynomials using FFT (placeholder for future implementation).
Currently falls back to naive multiplication.
"""
function fft_polynomial_multiply(p::Polynomial{F}, q::Polynomial{F}) where F
    # TODO: Implement FFT-based multiplication
    # For now, fall back to naive multiplication
    return p * q
end

"""
    roots_of_unity(n::Integer, ::Type{F}) where F<:FieldElem

Find the n-th roots of unity in field F (placeholder for future implementation).
"""
function roots_of_unity(n::Integer, ::Type{F}) where F<:FieldElem
    # TODO: Implement root finding for FFT
    throw(ArgumentError("FFT roots of unity not yet implemented"))
end
