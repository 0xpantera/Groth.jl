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
    # Evaluation domain (roots of unity)
    domain::Vector{F}
    # Target polynomial t(x) = ∏(x - ωⁱ)
    t::Polynomial{F}
    # Polynomials for each variable
    u::Vector{Polynomial{F}}  # Left polynomials (from L matrix)
    v::Vector{Polynomial{F}}  # Right polynomials (from R matrix)
    w::Vector{Polynomial{F}}  # Output polynomials (from O matrix)
end

"""
    get_roots_of_unity(n::Int, ::Type{F}) where F

Get n roots of unity in field F.
For simplicity, we use {1, 2, 3, ..., n} as evaluation points.
"""
function get_roots_of_unity(n::Int, ::Type{F}) where F
    return [F(i) for i in 1:n]
end

"""
    r1cs_to_qap(r1cs::R1CS{F}) where F

Convert an R1CS to QAP form using Lagrange interpolation.
"""
function r1cs_to_qap(r1cs::R1CS{F}) where F
    n = r1cs.num_constraints
    m = r1cs.num_vars
    
    # Get evaluation domain
    domain = get_roots_of_unity(n, F)
    
    # Compute target polynomial t(x) = ∏(x - ωⁱ)
    t = Polynomial{F}([one(F)])  # Start with polynomial 1
    for ω in domain
        # Multiply by (x - ω)
        factor = Polynomial{F}([-ω, one(F)])
        t = t * factor
    end
    
    # Initialize polynomial arrays
    u = Vector{Polynomial{F}}(undef, m)
    v = Vector{Polynomial{F}}(undef, m)
    w = Vector{Polynomial{F}}(undef, m)
    
    # For each variable, interpolate its polynomial from the constraint values
    for j in 1:m
        # Extract column j from each matrix
        u_points = [r1cs.L[i, j] for i in 1:n]
        v_points = [r1cs.R[i, j] for i in 1:n]
        w_points = [r1cs.O[i, j] for i in 1:n]
        
        # Interpolate polynomials
        u[j] = interpolate(domain, u_points)
        v[j] = interpolate(domain, v_points)
        w[j] = interpolate(domain, w_points)
    end
    
    return QAP{F}(m, n, r1cs.num_public, domain, t, u, v, w)
end

"""
    evaluate_qap(qap::QAP{F}, witness::Witness{F}, x::F) where F

Evaluate the QAP polynomials at point x with the given witness.
Returns (u(x), v(x), w(x)) where each is the linear combination of polynomials.
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

"""
    compute_h_polynomial(qap::QAP{F}, witness::Witness{F}) where F

Compute the quotient polynomial h(x) = (u(x) * v(x) - w(x)) / t(x).
"""
function compute_h_polynomial(qap::QAP{F}, witness::Witness{F}) where F
    w_vals = witness.values
    
    # Compute linear combinations of polynomials
    u_poly = zero(Polynomial{F})
    v_poly = zero(Polynomial{F})
    w_poly = zero(Polynomial{F})
    
    for i in 1:qap.num_vars
        u_poly = u_poly + w_vals[i] * qap.u[i]
        v_poly = v_poly + w_vals[i] * qap.v[i]
        w_poly = w_poly + w_vals[i] * qap.w[i]
    end
    
    # Compute p(x) = u(x) * v(x) - w(x)
    p_poly = u_poly * v_poly - w_poly
    
    # Divide by target polynomial
    h_poly = polynomial_division(p_poly, qap.t)
    
    return h_poly
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
            remainder[end - i] -= coeff * divisor.coeffs[end - i]
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