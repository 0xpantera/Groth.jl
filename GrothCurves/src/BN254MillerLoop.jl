"""
BN254 Miller loop implementation for the optimal ate pairing.

This module implements the Miller loop algorithm which is the core
computation for the BN254 pairing. The ate pairing uses a shorter
loop than the Tate pairing for better performance.

References:
- Aranha et al., "The Realm of the Pairings", 2013
  https://eprint.iacr.org/2013/722.pdf
  Sections 4.3-4.4 for projective doubling/addition and M-twist line placement
  
- LambdaClass, "How we implemented the BN254 Ate pairing in lambdaworks", 2024
  https://blog.lambdaclass.com/how-we-implemented-the-bn254-ate-pairing-in-lambdaworks/
  Shows the exact BN254 coordinate mapping and line evaluation formulas
  
- Arkworks implementation: ark-bn254
  https://github.com/arkworks-rs/algebra
  Reference for sparse 0-1-4 multiplication (mul_by_014)
"""

using GrothAlgebra
using StaticArrays

# Import what we need
import GrothCurves: conjugate, BN254_PRIME

# BN254 parameter u = 4965661367192848881
# For optimal ate pairing with Frobenius corrections, we iterate over u only
# The factor of 6 and addition of 2 are handled by the two Frobenius correction lines
const BN254_U = BigInt(4965661367192848881)
const BN254_U_IS_NEGATIVE = false  # u > 0 for BN254
const ATE_LOOP_COUNT = BN254_U  # Just u, not 6u+2 (corrections handle the rest)

# G2 curve parameter: b' = 3/(u+9) for the twist curve y² = x³ + b'
const G2_B_TWIST = Fp2Element(3, 0) / Fp2Element(9, 1)  # 3/(u+9)

# ============================================================================
# Frobenius Constants for G2 on the M-twist
# ============================================================================

# ξ = 9 + u ∈ Fp2 (non-residue for the sextic twist)
# Using the same constant that Fp6 uses internally
const XI_FP2 = Fp2Element(bn254_field(9), bn254_field(1))

# Precompute gamma constants for G2 Frobenius on the M-twist
# γ1[j] = ξ^((p-1)j/6), γ2[j] = γ1[j]*conj(γ1[j]), γ3[j] = γ1[j]*γ2[j]
function compute_g2_frobenius_constants()
    pow = div(BN254_PRIME - 1, 6)
    
    γ1 = Vector{Fp2Element}(undef, 5)
    for j in 1:5
        γ1[j] = XI_FP2^(pow * BigInt(j))
    end
    
    γ2 = Vector{Fp2Element}(undef, 5)
    for j in 1:5
        γ2[j] = γ1[j] * conjugate(γ1[j])
    end
    
    γ3 = Vector{Fp2Element}(undef, 5)
    for j in 1:5
        γ3[j] = γ1[j] * γ2[j]
    end
    
    return γ1, γ2, γ3
end

# Compute once at module load
const (G2_GAMMA1, G2_GAMMA2, G2_GAMMA3) = compute_g2_frobenius_constants()

"""
    LineCoeffs

Raw coefficients for a line function in the Miller loop, before evaluation at P.
For M-twist, the line has form: a*x + b*y + c where x,y are indeterminates.
"""
struct LineCoeffs
    a::Fp2Element  # x coefficient
    b::Fp2Element  # y coefficient  
    c::Fp2Element  # constant term
end

"""
    doubling_step(T::G2Point) -> (G2Point, LineCoeffs)

Compute the doubling step in the Miller loop.
Returns 2T and the raw line coefficients for the tangent line at T.

Uses M-twist formulas for BN254 with Jacobian coordinates.
"""
function doubling_step(T::G2Point)
    # Handle point at infinity
    if iszero(T)
        return (T, LineCoeffs(zero(Fp2Element), zero(Fp2Element), one(Fp2Element)))
    end
    
    # Get coordinates of T in Jacobian form
    X, Y, Z = x_coord(T), y_coord(T), z_coord(T)
    
    # Precomputations
    X2 = X^2
    X3 = X2 * X
    Y2 = Y^2
    Y4 = Y2^2
    Z2 = Z^2
    Z3 = Z2 * Z
    Z6 = Z2^3  # More efficient than Z3 * Z3
    
    # Step 1: Compute line coefficients for Jacobian coordinates BEFORE mutating T
    # For point T = (X:Y:Z) where x = X/Z², y = Y/Z³
    # The correct line coefficients (up to scaling) are:
    a = Fp2Element(3) * X2 * Z2                          # x coefficient  
    b = -Fp2Element(2) * Y * Z3                          # y coefficient
    c = -X3 + Fp2Element(2) * G2_B_TWIST * Z6            # constant term
    
    # Step 2: NOW compute 2T
    S = Fp2Element(4) * X * Y2
    M = Fp2Element(3) * X2  
    T_val = M^2 - Fp2Element(2) * S
    
    # New coordinates for 2T
    X_2T = T_val
    Y_2T = M * (S - T_val) - Fp2Element(8) * Y4
    Z_2T = Fp2Element(2) * Y * Z
    
    result_point = G2Point(X_2T, Y_2T, Z_2T)
    coeffs = LineCoeffs(a, b, c)
    
    return (result_point, coeffs)
end

"""
    addition_step(T::G2Point, Q::G2Point) -> (G2Point, LineCoeffs)

Compute the addition step in the Miller loop.
Returns T+Q and the raw line coefficients for the line through T and Q.

This implements the first equation from page 14 of "The Realm of the Pairings".
Uses mixed addition with Q in affine coordinates.
"""
function addition_step(T::G2Point, Q::G2Point)
    # Handle special cases
    if iszero(T)
        return (Q, LineCoeffs(zero(Fp2Element), zero(Fp2Element), one(Fp2Element)))
    end
    if iszero(Q)
        return (T, LineCoeffs(zero(Fp2Element), zero(Fp2Element), one(Fp2Element)))
    end
    
    # If T == Q, use doubling instead
    if T == Q
        return doubling_step(T)
    end
    
    # Get coordinates
    X, Y, Z = x_coord(T), y_coord(T), z_coord(T)
    xQ, yQ = to_affine(Q)  # Q in affine (common for pairing)
    
    # Mixed addition formulas (T projective, Q affine)
    Z2 = Z^2
    Z3 = Z2 * Z
    
    # Bring Q to same projective denominator as T
    U2 = xQ * Z2
    S2 = yQ * Z3
    
    # Differences
    H = U2 - X  # x-coordinate difference
    R = S2 - Y  # y-coordinate difference
    
    # Check if points are equal (shouldn't happen if we checked T == Q above)
    if iszero(H) && iszero(R)
        return doubling_step(T)
    end
    
    # Step 1: Compute line coefficients BEFORE mutating T
    # For Jacobian coordinates where T = (X:Y:Z) with x = X/Z², y = Y/Z³
    # and Q = (xQ, yQ) in affine
    # The line through T and Q has coefficients (with proper Z scaling):
    a = R                            # coefficient for x
    b = -Z * H                       # coefficient for y (with Z scaling)
    c = Z * (H * yQ) - R * xQ        # constant term (with Z scaling)
    
    # Step 2: NOW compute T + Q (mutate the point)
    H2 = H^2
    H3 = H2 * H
    
    X3 = R^2 - H3 - Fp2Element(2) * X * H2
    Y3 = R * (X * H2 - X3) - Y * H3
    Z3 = Z * H
    
    result_point = G2Point(X3, Y3, Z3)
    coeffs = LineCoeffs(a, b, c)
    
    return (result_point, coeffs)
end

"""
    evaluate_line(coeffs::LineCoeffs, P::G1Point) -> Fp12Element

Evaluate the line function at point P and embed into Fp12.
Uses sparse 0-1-4 embedding matching arkworks' mul_by_014.

The line a*x + b*y + c is evaluated at P = (x_P, y_P) to give
a*x_P + b*y_P + c, which is then embedded into Fp12.
"""
function evaluate_line(coeffs::LineCoeffs, P::G1Point)
    # Get affine coordinates of P
    P_x, P_y = if iszero(P)
        (zero(BN254Field), zero(BN254Field))
    else
        to_affine(P)
    end
    
    # Convert to Fp2 for computation
    xP = Fp2Element(P_x, zero(BN254Field))
    yP = Fp2Element(P_y, zero(BN254Field))
    
    # M-twist embedding for ψ: x→xv, y→yw
    # Line ℓ(ψ(P)) = (a*xP*v + c) + (b*yP)*w
    # This gives us Fp12((a*xP*v + c), b*yP)
    # In Fp6 coordinates:
    #   c0 = (c, a*xP, 0) since v = (0,1,0)
    #   c1 = (b*yP, 0, 0)
    
    # Build sparse Fp12 element with correct M-twist mapping
    c0_fp6 = Fp6Element(coeffs.c, coeffs.a * xP, zero(Fp2Element))
    c1_fp6 = Fp6Element(coeffs.b * yP, zero(Fp2Element), zero(Fp2Element))
    
    return Fp12Element(c0_fp6, c1_fp6)
end

"""
    frobenius_g2(Q::G2Point, power::Int) -> G2Point

Apply Frobenius endomorphism to a G2 point.
For G2 over Fp2, this conjugates coordinates and applies twist constants.
"""
function frobenius_g2(Q::G2Point, power::Int)
    if iszero(Q)
        return Q
    end
    
    # Get affine coordinates
    x, y = to_affine(Q)
    
    if power == 1
        # π: (x,y) → (x̄·γ₁,₂, ȳ·γ₁,₃)
        # Apply Frobenius (conjugation in Fp2) and multiply by gamma constants
        x_final = conjugate(x) * G2_GAMMA1[2]
        y_final = conjugate(y) * G2_GAMMA1[3]
        
        return G2Point(x_final, y_final)
        
    elseif power == 2
        # π²: (x,y) → (x·γ₂,₂, y·γ₂,₃)
        # No conjugation needed for even powers
        x_final = x * G2_GAMMA2[2]
        y_final = y * G2_GAMMA2[3]
        
        return G2Point(x_final, y_final)
        
    elseif power == 3
        # π³: (x,y) → (x̄·γ₃,₂, ȳ·γ₃,₃)
        # Conjugation needed for odd powers
        x_final = conjugate(x) * G2_GAMMA3[2]
        y_final = conjugate(y) * G2_GAMMA3[3]
        
        return G2Point(x_final, y_final)
        
    else
        # For other powers, compose
        result = Q
        for _ in 1:power
            result = frobenius_g2(result, 1)
        end
        return result
    end
end

"""
    miller_loop(P::G1Point, Q::G2Point) -> Fp12Element

Compute the Miller loop for the optimal ate pairing.
This is the main computational component of the pairing.

The loop iterates over u (not 6u+2), with two Frobenius correction lines
at the end to achieve the full optimal ate pairing.
"""
function miller_loop(P::G1Point, Q::G2Point)
    # Handle special cases
    if iszero(P) || iszero(Q)
        return one(Fp12Element)
    end
    
    # Initialize
    T = Q
    f = one(Fp12Element)
    
    # Get binary representation of loop count (just u for optimal ate with corrections)
    # We process from MSB to LSB (excluding the most significant 1)
    loop_bits = digits(Bool, ATE_LOOP_COUNT, base=2)
    
    # Start from the second-most significant bit
    for i in (length(loop_bits)-1):-1:1
        # Double step (always performed)
        T_new, line_coeffs = doubling_step(T)
        T = T_new
        
        # Square f and multiply by line function evaluated at P
        f = f^2 * evaluate_line(line_coeffs, P)
        
        # Addition step (if bit is 1)
        if loop_bits[i]
            T_new, line_coeffs = addition_step(T, Q)
            T = T_new
            f = f * evaluate_line(line_coeffs, P)
        end
    end
    
    # For BN254 ate pairing, we need two Frobenius-based correction steps
    # These are essential for bilinearity to work correctly
    
    # Compute Frobenius images of Q
    # π(Q) applies Frobenius to the x,y coordinates
    Q_pi = frobenius_g2(Q, 1)   # π(Q)
    Q_pi2 = frobenius_g2(Q, 2)  # π²(Q)
    
    # Correction step 1: Line through T and π(Q)
    T1, line1 = addition_step(T, Q_pi)
    f = f * evaluate_line(line1, P)
    
    # Correction step 2: Line through T1 and -π²(Q)
    neg_Q_pi2 = -Q_pi2  # Use proper group negation
    T2, line2 = addition_step(T1, neg_Q_pi2)
    f = f * evaluate_line(line2, P)
    
    return f
end

# Export functions
export LineCoeffs, doubling_step, addition_step, evaluate_line, miller_loop
export frobenius_g2