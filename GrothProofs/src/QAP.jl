"""
QAP (Quadratic Arithmetic Program) conversion from R1CS.

Converts an R1CS into polynomial form using Lagrange interpolation.
"""

using GrothAlgebra
using GrothCurves

# R1CS is included by the parent module

"""
    QAP{F}

Quadratic Arithmetic Program representation.
"""
struct QAP{F}
    # Number of variables
    num_vars::Int
    # Number of constraints
    num_constraints::Int
    # Number of public inputs
    num_public::Int
    # Evaluation domain (power-of-two roots of unity)
    domain::EvaluationDomain{F}
    # Active constraint points (prefix subset of the domain)
    points::Vector{F}
    # Coset domain used for FFT-based evaluation
    coset_domain::EvaluationDomain{F}
    # Reciprocal evaluations of the vanishing polynomial on the coset domain
    vanishing_coset_inv::Vector{F}
    # Target polynomial ``t(x) = x^N - 1`` for the full evaluation domain.
    t::Polynomial{F}
    # Polynomials for each variable
    u::Vector{Polynomial{F}}  # Left polynomials (from L matrix)
    v::Vector{Polynomial{F}}  # Right polynomials (from R matrix)
    w::Vector{Polynomial{F}}  # Output polynomials (from O matrix)
end

"""
    get_roots_of_unity(n::Int, ::Type{F}) where F

Construct a power-of-two evaluation domain and return the first `n` domain
points.
"""
function get_roots_of_unity(n::Int, ::Type{F}) where F
    n > 0 || throw(ArgumentError("Number of roots must be positive"))
    size = 1
    while size < n
        size <<= 1
    end
    domain = EvaluationDomain(F, size)
    points = Vector{F}(undef, n)
    points[1] = one(F)
    ω = domain.generator
    for i in 2:n
        points[i] = points[i-1] * ω
    end
    return points, domain
end

default_coset_offset(::Type{BN254Fr}) = bn254_fr(5)
default_coset_offset(::Type{F}) where F = throw(ArgumentError("No coset offset configured for field"))

"""
    r1cs_to_qap(r1cs::R1CS{F}) where F

Convert an R1CS to QAP form using Lagrange interpolation.
"""
function r1cs_to_qap(r1cs::R1CS{F}) where F
    n = r1cs.num_constraints
    m = r1cs.num_vars
    domain_slots = n + r1cs.num_public

    # Match arkworks' domain shape: active constraints followed by public-input
    # selector slots, rounded up to a power of two.
    domain_points, domain = get_roots_of_unity(domain_slots, F)
    points = domain_points[1:n]
    coset_domain = get_coset(domain, default_coset_offset(F))

    # Full-domain vanishing polynomial t(x) = x^N - 1.
    t_coeffs = fill(zero(F), domain.size + 1)
    t_coeffs[1] = -one(F)
    t_coeffs[end] = one(F)
    t = Polynomial{F}(t_coeffs)

    # Initialize polynomial arrays
    u = Vector{Polynomial{F}}(undef, m)
    v = Vector{Polynomial{F}}(undef, m)
    w = Vector{Polynomial{F}}(undef, m)

    # For each variable, interpolate its polynomial from the constraint values
    for j in 1:m
        u_samples = fill(zero(F), domain.size)
        v_samples = fill(zero(F), domain.size)
        w_samples = fill(zero(F), domain.size)

        @inbounds for i in 1:n
            u_samples[i] = r1cs.L[i, j]
            v_samples[i] = r1cs.R[i, j]
            w_samples[i] = r1cs.O[i, j]
        end
        if j <= r1cs.num_public
            u_samples[n + j] = one(F)
        end

        u[j] = interpolate_fft(domain, u_samples)
        v[j] = interpolate_fft(domain, v_samples)
        w[j] = interpolate_fft(domain, w_samples)
    end

    coset_vanishing = coset_offset_pow_size(coset_domain) - one(F)
    iszero(coset_vanishing) && error("Coset vanishing polynomial evaluates to zero; choose a different offset")
    vanishing_coset_inv = fill(inv(coset_vanishing), coset_domain.size)

    return QAP{F}(m, n, r1cs.num_public, domain, points, coset_domain, vanishing_coset_inv, t, u, v, w)
end

"""
    evaluate_qap(qap::QAP{F}, witness::Witness{F}, x::F) where F

Evaluate the QAP polynomials at point x with the given witness.
Returns ``(u(x), v(x), w(x))`` where each term is a witness-weighted linear
combination of column polynomials.
"""
function evaluate_qap(qap::QAP{F}, witness::Witness{F}, x::F) where F
    w_vals = witness.values

    # Compute linear combinations
    u_eval = zero(F)
    v_eval = zero(F)
    w_eval = zero(F)

    for i in 1:qap.num_vars
        u_eval += w_vals[i] * evaluate(qap.u[i], x)
        v_eval += w_vals[i] * evaluate(qap.v[i], x)
        w_eval += w_vals[i] * evaluate(qap.w[i], x)
    end

    return (u_eval, v_eval, w_eval)
end

function combined_qap_polynomials(qap::QAP{F}, witness::Witness{F}) where F
    w_vals = witness.values

    u_poly = zero(Polynomial{F})
    v_poly = zero(Polynomial{F})
    w_poly = zero(Polynomial{F})

    for i in 1:qap.num_vars
        u_poly = u_poly + w_vals[i] * qap.u[i]
        v_poly = v_poly + w_vals[i] * qap.v[i]
        w_poly = w_poly + w_vals[i] * qap.w[i]
    end

    return u_poly, v_poly, w_poly
end

function compute_h_polynomial_dense(qap::QAP{F}, u_poly::Polynomial{F}, v_poly::Polynomial{F}, w_poly::Polynomial{F}) where F
    p_poly = fft_polynomial_multiply(u_poly, v_poly) - w_poly
    return polynomial_division(p_poly, qap.t)
end

function compute_h_polynomial_dense(qap::QAP{F}, witness::Witness{F}) where F
    u_poly, v_poly, w_poly = combined_qap_polynomials(qap, witness)
    return compute_h_polynomial_dense(qap, u_poly, v_poly, w_poly)
end

function compute_h_polynomial_coset(qap::QAP{F}, u_poly::Polynomial{F}, v_poly::Polynomial{F}, w_poly::Polynomial{F}) where F
    coset = qap.coset_domain
    u_eval = fft(u_poly.coeffs, coset)
    v_eval = fft(v_poly.coeffs, coset)
    w_eval = fft(w_poly.coeffs, coset)
    h_eval = similar(u_eval)
    vanishing_inv = qap.vanishing_coset_inv
    for i in eachindex(u_eval)
        numerator = u_eval[i] * v_eval[i] - w_eval[i]
        h_eval[i] = numerator * vanishing_inv[i]
    end
    coeffs = GrothAlgebra.ifft!(h_eval, coset)
    return Polynomial{F}(coeffs)
end

function compute_h_polynomial_coset(qap::QAP{F}, witness::Witness{F}) where F
    u_poly, v_poly, w_poly = combined_qap_polynomials(qap, witness)
    return compute_h_polynomial_coset(qap, u_poly, v_poly, w_poly)
end

function compute_h_polynomial_checked(qap::QAP{F}, witness::Witness{F}) where F
    u_poly, v_poly, w_poly = combined_qap_polynomials(qap, witness)
    dense_result = compute_h_polynomial_dense(qap, u_poly, v_poly, w_poly)
    coset_result = compute_h_polynomial_coset(qap, u_poly, v_poly, w_poly)

    coset_result == dense_result || error("Coset Groth16 path diverged from dense fallback")

    return coset_result
end

"""
    compute_h_polynomial(qap::QAP{F}, witness::Witness{F}; use_coset::Bool=true) where F

Compute the quotient polynomial ``h(x)``:

```math
h(x) = \\frac{u(x)v(x) - w(x)}{t(x)}
```

By default this computes both the dense quotient and the coset quotient and
asserts that they agree. Pass `use_coset=false` for the dense quotient only.
The prover's hot path calls the internal coset-only helper directly.
"""
function compute_h_polynomial(qap::QAP{F}, witness::Witness{F}; use_coset::Bool=true) where F
    return use_coset ? compute_h_polynomial_checked(qap, witness) : compute_h_polynomial_dense(qap, witness)
end

"""
    polynomial_division(dividend::Polynomial{F}, divisor::Polynomial{F}) where F

Divide polynomial dividend by divisor, returning the quotient.
Assumes the division is exact (no remainder).
"""
function polynomial_division(dividend::Polynomial{F}, divisor::Polynomial{F}) where F
    if is_zero(divisor)
        error("Division by zero polynomial")
    end

    # Working copy of dividend
    remainder = copy(dividend.coeffs)
    quotient = F[]

    divisor_degree = degree(divisor)
    divisor_lead = leading_coefficient(divisor)

    while length(remainder) >= length(divisor.coeffs) && !is_zero(Polynomial{F}(remainder))
        # Get the leading coefficient of current remainder
        coeff = remainder[end] / divisor_lead
        push!(quotient, coeff)

        # Subtract divisor * coeff from remainder
        for i in 0:divisor_degree
            remainder[end-i] -= coeff * divisor.coeffs[end-i]
        end

        # Remove the leading term
        pop!(remainder)
    end

    # Reverse quotient coefficients (we built it backwards)
    reverse!(quotient)

    # Check that remainder is zero (exact division)
    if !all(iszero, remainder)
        error("Polynomial division has non-zero remainder")
    end

    return Polynomial{F}(quotient)
end

# Export types and functions
export QAP, r1cs_to_qap, evaluate_qap, compute_h_polynomial
export get_roots_of_unity, polynomial_division
