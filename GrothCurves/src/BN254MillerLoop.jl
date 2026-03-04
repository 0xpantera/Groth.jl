# BN254 Miller loop implementation for the optimal ate pairing.
#
# Implements the Miller loop core of the BN254 pairing. The ate loop is shorter
# than the Tate loop, improving performance.
#
# References:
# - Aranha et al., "The Realm of the Pairings", 2013.
# - LambdaClass, "How we implemented the BN254 Ate pairing in lambdaworks",
#   2024.
# - Arkworks `ark-bn254` for sparse multiplication helpers.

# using GrothAlgebra
# using StaticArrays

# Import what we need
import GrothCurves: conjugate, BN254_PRIME

# BN254 parameter u = 4965661367192848881
# For optimal ate pairing with Frobenius corrections, we iterate over u only
# The factor of 6 and addition of 2 are handled by the two Frobenius correction lines
const BN254_U = BigInt(4965661367192848881)
const BN254_U_IS_NEGATIVE = false  # u > 0 for BN254
# NAF representation of u for the ate loop (from arkworks)
# This is the signed binary representation where digits are in {0, 1, -1}
const ATE_LOOP_COUNT_NAF = Int8[
    0, 0, 0, 1, 0, 1, 0, -1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1, 0, -1, 0, 0, 0, 1, 0, -1, 0, 0, 0,
    0, -1, 0, 0, 1, 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0, 0, -1, 0, 1, 0, -1, 0, 0, 0, -1, 0,
    -1, 0, 0, 0, 1, 0, 1, 1,
]

# ============================================================================
# p-power Endomorphism Constants for G2 Pairing
# ============================================================================

# These are the correct p-power endomorphism coefficients from arkworks
# PSI_X = (u+9)^((p-1)/3) for x-coordinate multiplication after Frobenius
const P_POWER_ENDOMORPHISM_COEFF_0 = Fp2Element(
    bn254_fq(parse(BigInt, "21575463638280843010398324269430826099269044274347216827212613867836435027261")),
    bn254_fq(parse(BigInt, "10307601595873709700152284273816112264069230130616436755625194854815875713954"))
)

# PSI_Y = (u+9)^((p-1)/2) for y-coordinate multiplication after Frobenius
const P_POWER_ENDOMORPHISM_COEFF_1 = Fp2Element(
    bn254_fq(parse(BigInt, "2821565182194536844548159561693502659359617185244120367078079554186484126554")),
    bn254_fq(parse(BigInt, "3505843767911556378687030309984248845540243509899259641013678093033130930403"))
)

# Note: For π²(Q), we compose the endomorphism twice rather than using squared coefficients

"""
    LineCoeffs

Raw coefficients for a line function in the Miller loop, before evaluation at P.
For D-twist, the line has form: a*y + b*x + c where x,y are indeterminates.
Note the order difference from M-twist!
"""
struct LineCoeffs
    a::Fp2Element  # y coefficient (for D-twist)
    b::Fp2Element  # x coefficient (for D-twist)
    c::Fp2Element  # constant term
end

"""
    doubling_step(T::G2Point) -> (G2Point, LineCoeffs)

Compute the doubling step in the Miller loop.
Returns 2T and the raw line coefficients for the tangent line at T.

Uses D-twist formulas for BN254 with Jacobian coordinates.
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

    # Step 1: Compute line coefficients for D-twist
    # For point T = (X:Y:Z) where x = X/Z², y = Y/Z³
    # D-twist coefficients (note the different order from M-twist):
    a = -Fp2Element(2) * Y * Z3                          # y coefficient (D-twist)
    b = Fp2Element(3) * X2 * Z2                          # x coefficient (D-twist)
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

    # Step 1: Compute line coefficients for D-twist BEFORE mutating T
    # For Jacobian coordinates where T = (X:Y:Z) with x = X/Z², y = Y/Z³
    # and Q = (xQ, yQ) in affine
    # D-twist coefficients (note the different order from M-twist):
    a = -Z * H                       # coefficient for y (D-twist)
    b = R                            # coefficient for x (D-twist)
    c = Z * (H * yQ) - R * xQ        # constant term

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
        (zero(BN254Fq), zero(BN254Fq))
    else
        to_affine(P)
    end

    # Convert to Fp2 for computation
    xP = Fp2Element(P_x, zero(BN254Fq))
    yP = Fp2Element(P_y, zero(BN254Fq))

    # D-twist embedding - different from M-twist!
    # For D-twist, the line evaluation gives a sparse Fp12 element
    # with non-zero entries at positions (0, 3, 4) using mul_by_034 format
    # Line ℓ(P) = a*yP + b*xP*w + c*w²
    # This maps to Fp12 coordinates as:
    #   c0 = (a*yP, 0, 0)
    #   c1 = (b*xP, c, 0)

    # Build sparse Fp12 element with correct D-twist mapping
    c0_fp6 = Fp6Element(coeffs.a * yP, zero(Fp2Element), zero(Fp2Element))
    c1_fp6 = Fp6Element(coeffs.b * xP, coeffs.c, zero(Fp2Element))

    return Fp12Element(c0_fp6, c1_fp6)
end

"""
    p_power_endomorphism(Q::G2Point) -> G2Point

Apply the p-power endomorphism to a G2 point.
This is the untwist-Frobenius-twist endomorphism: ψ ∘ π_p ∘ ψ^(-1)
Maps (x,y) → (x^p * PSI_X, y^p * PSI_Y)
"""
function p_power_endomorphism(Q::G2Point)
    if iszero(Q)
        return Q
    end

    # Get affine coordinates
    x, y = to_affine(Q)

    # Apply Frobenius (conjugation in Fp2)
    x_frob = conjugate(x)
    y_frob = conjugate(y)

    # Multiply by the p-power endomorphism coefficients
    x_final = x_frob * P_POWER_ENDOMORPHISM_COEFF_0
    y_final = y_frob * P_POWER_ENDOMORPHISM_COEFF_1

    return G2Point(x_final, y_final)
end

"""
    frobenius_g2(Q::G2Point, power::Int) -> G2Point

Apply powers of the p-power endomorphism to a G2 point.
For power=1: applies p_power_endomorphism once
For power=2: applies it twice (π²)
"""
function frobenius_g2(Q::G2Point, power::Int)
    if iszero(Q)
        return Q
    end

    if power == 1
        # π(Q) = p-power endomorphism
        return p_power_endomorphism(Q)

    elseif power == 2
        # π²(Q) = apply p-power endomorphism twice
        # We need to compose: π²(Q) = π(π(Q))
        intermediate = p_power_endomorphism(Q)
        return p_power_endomorphism(intermediate)

    else
        # For other powers, compose
        result = Q
        for _ in 1:power
            result = p_power_endomorphism(result)
        end
        return result
    end
end

"""
    miller_loop(::BN254Engine, P::G1Point, Q::G2Point) -> Fp12Element

Compute the Miller loop for the optimal ate pairing using the BN254 backend.
This is the main computational component of the pairing.

The loop iterates over u (not 6u+2), with two Frobenius correction lines
at the end to achieve the full optimal ate pairing.
"""
function miller_loop(::BN254Engine, P::G1Point, Q::G2Point)
    # Handle special cases
    if iszero(P) || iszero(Q)
        return one(Fp12Element)
    end

    # Initialize
    T = Q
    f = one(Fp12Element)

    # Process NAF representation from MSB to LSB
    # The loop matches arkworks: iterate from len-1 down to 1, indexing bits at i-1
    for i in (length(ATE_LOOP_COUNT_NAF)):-1:2
        # Square f (except on the first iteration)
        if i != length(ATE_LOOP_COUNT_NAF)
            f = f^2
        end

        # Always perform doubling step
        T_new, line_coeffs = doubling_step(T)
        T = T_new
        f = f * evaluate_line(line_coeffs, P)

        # Addition/subtraction step based on NAF digit at index i-1
        bit = ATE_LOOP_COUNT_NAF[i-1]
        if bit == 1
            T_new, line_coeffs = addition_step(T, Q)
            T = T_new
            f = f * evaluate_line(line_coeffs, P)
        elseif bit == -1
            T_new, line_coeffs = addition_step(T, -Q)
            T = T_new
            f = f * evaluate_line(line_coeffs, P)
        end
        # If bit == 0, we do nothing extra
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

"""
    miller_loop(P::G1Point, Q::G2Point)

Convenience overload that reuses the shared BN254 engine.
"""
miller_loop(P::G1Point, Q::G2Point) = miller_loop(BN254_ENGINE, P, Q)

# Export functions
export LineCoeffs, doubling_step, addition_step, evaluate_line, miller_loop
export frobenius_g2, p_power_endomorphism
export ATE_LOOP_COUNT_NAF
