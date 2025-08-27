using GrothCurves
using Test

# BN254 parameters
const u_positive = BigInt(4965661367192848881)
const u_negative = -BigInt(4965661367192848881)
const p = GrothCurves.BN254_PRIME
const r = parse(BigInt, "21888242871839275222246405745257275088548364400416034343698204186575808495617")

println("=== Testing BN254 with NEGATIVE u parameter ===\n")
println("u (positive) = ", u_positive)
println("u (negative) = ", u_negative)
println("6u+2 (positive) = ", 6 * u_positive + 2)
println("6u+2 (negative) = ", 6 * u_negative + 2)
println()

# Miller loop with negative u
function miller_loop_negative_u(P::G1Point, Q::G2Point)
    if iszero(P) || iszero(Q)
        return one(Fp12Element)
    end

    T = Q
    f = one(Fp12Element)

    # Use POSITIVE u for the loop (we'll handle sign later)
    loop_bits = digits(Bool, u_positive, base=2)

    for i in (length(loop_bits)-1):-1:1
        T_new, line_coeffs = GrothCurves.doubling_step(T)
        T = T_new
        f = f^2 * GrothCurves.evaluate_line(line_coeffs, P)

        if loop_bits[i]
            T_new, line_coeffs = GrothCurves.addition_step(T, Q)
            T = T_new
            f = f * GrothCurves.evaluate_line(line_coeffs, P)
        end
    end

    # Now T = [|u|]Q (positive u)
    # For negative u, we need T = [-u]Q = -[|u|]Q
    T = -T
    # And we need to conjugate f
    f = GrothCurves.conjugate(f)

    # Apply Frobenius corrections for NEGATIVE u
    Q_pi = GrothCurves.frobenius_g2(Q, 1)
    Q_pi2 = GrothCurves.frobenius_g2(Q, 2)

    # For negative u, the corrections might be different
    # Try: line from T to π(Q)
    T1, line1 = GrothCurves.addition_step(T, Q_pi)
    f = f * GrothCurves.evaluate_line(line1, P)

    # And line from T1 to -π²(Q)
    neg_Q_pi2 = -Q_pi2
    T2, line2 = GrothCurves.addition_step(T1, neg_Q_pi2)
    f = f * GrothCurves.evaluate_line(line2, P)

    return f
end

# Alternative: Different correction for negative u
function miller_loop_negative_u_alt(P::G1Point, Q::G2Point)
    if iszero(P) || iszero(Q)
        return one(Fp12Element)
    end

    T = Q
    f = one(Fp12Element)

    loop_bits = digits(Bool, u_positive, base=2)

    for i in (length(loop_bits)-1):-1:1
        T_new, line_coeffs = GrothCurves.doubling_step(T)
        T = T_new
        f = f^2 * GrothCurves.evaluate_line(line_coeffs, P)

        if loop_bits[i]
            T_new, line_coeffs = GrothCurves.addition_step(T, Q)
            T = T_new
            f = f * GrothCurves.evaluate_line(line_coeffs, P)
        end
    end

    # T = [|u|]Q
    # For negative u: T = -[|u|]Q
    T = -T
    f = GrothCurves.conjugate(f)

    # Try different corrections for negative u
    # Based on the formula: f_{-u,Q}(P) * line_{[-u]Q, π²(Q)}(P) * line_{[-u]Q+π²(Q), π(Q)}(P)
    Q_pi = GrothCurves.frobenius_g2(Q, 1)
    Q_pi2 = GrothCurves.frobenius_g2(Q, 2)

    # First correction: add π²(Q)
    T1, line1 = GrothCurves.addition_step(T, Q_pi2)
    f = f * GrothCurves.evaluate_line(line1, P)

    # Second correction: add π(Q)
    T2, line2 = GrothCurves.addition_step(T1, Q_pi)
    f = f * GrothCurves.evaluate_line(line2, P)

    return f
end

# Test if we're using the wrong sign
function verify_u_sign()
    println("=== Verifying correct sign for u ===\n")

    Q = g2_generator()

    # Compute [u]Q and [-u]Q
    Q_u_pos = u_positive * Q
    Q_u_neg = u_negative * Q

    # Get Frobenius images
    Q_pi = GrothCurves.frobenius_g2(Q, 1)
    Q_pi2 = GrothCurves.frobenius_g2(Q, 2)

    # Test various combinations to see which gives [6u+2]Q or [6u-2]Q
    println("Testing with POSITIVE u:")
    target_pos = (6 * u_positive + 2) * Q
    result_pos = Q_u_pos + Q_pi - Q_pi2
    println("  [u]Q + π(Q) - π²(Q) == [6u+2]Q? ", result_pos == target_pos)

    println("\nTesting with NEGATIVE u:")
    target_neg = (6 * u_negative + 2) * Q  # This is actually 6*(-u) + 2 = -6u + 2
    result_neg = Q_u_neg + Q_pi - Q_pi2
    println("  [-u]Q + π(Q) - π²(Q) == [6*(-u)+2]Q? ", result_neg == target_neg)

    # Also test with different correction formulas
    println("\nAlternative formulas:")

    # For negative u, maybe it's: [-u]Q + π²(Q) + π(Q)
    result_alt1 = Q_u_neg + Q_pi2 + Q_pi
    target_alt = (-6 * u_positive + 2) * Q
    println("  [-u]Q + π²(Q) + π(Q) == [-6u+2]Q? ", result_alt1 == target_alt)

    # Or maybe we need to check the absolute value
    target_abs = (6 * abs(u_negative) + 2) * Q
    println("  [-u]Q + π²(Q) + π(Q) == [6|u|+2]Q? ", result_alt1 == target_abs)

    return Q_u_pos, Q_u_neg, Q_pi, Q_pi2
end

# Pairing with negative u
function pairing_negative_u(P::G1Point, Q::G2Point)
    f = miller_loop_negative_u(P, Q)
    return GrothCurves.final_exponentiation(f)
end

function pairing_negative_u_alt(P::G1Point, Q::G2Point)
    f = miller_loop_negative_u_alt(P, Q)
    return GrothCurves.final_exponentiation(f)
end

# Test bilinearity
function test_bilinearity_negative()
    println("\n=== Testing Bilinearity with Negative u ===\n")

    P = g1_generator()
    Q = g2_generator()

    # Test original (positive u)
    println("Original implementation (positive u):")
    e_orig = GrothCurves.pairing(P, Q)
    e_2P_orig = GrothCurves.pairing(2 * P, Q)
    e_P2Q_orig = GrothCurves.pairing(P, 2 * Q)
    println("  e(2P, Q) == e(P, Q)^2: ", e_2P_orig == e_orig^2)
    println("  e(P, 2Q) == e(P, Q)^2: ", e_P2Q_orig == e_orig^2)
    println("  e(2P, Q) == e(P, 2Q): ", e_2P_orig == e_P2Q_orig)

    # Test with negative u
    println("\nNegative u implementation:")
    e_neg = pairing_negative_u(P, Q)
    e_2P_neg = pairing_negative_u(2 * P, Q)
    e_P2Q_neg = pairing_negative_u(P, 2 * Q)
    println("  e(2P, Q) == e(P, Q)^2: ", e_2P_neg == e_neg^2)
    println("  e(P, 2Q) == e(P, Q)^2: ", e_P2Q_neg == e_neg^2)
    println("  e(2P, Q) == e(P, 2Q): ", e_2P_neg == e_P2Q_neg)

    # Test alternative negative u
    println("\nAlternative negative u implementation:")
    e_alt = pairing_negative_u_alt(P, Q)
    e_2P_alt = pairing_negative_u_alt(2 * P, Q)
    e_P2Q_alt = pairing_negative_u_alt(P, 2 * Q)
    println("  e(2P, Q) == e(P, Q)^2: ", e_2P_alt == e_alt^2)
    println("  e(P, 2Q) == e(P, Q)^2: ", e_P2Q_alt == e_alt^2)
    println("  e(2P, Q) == e(P, 2Q): ", e_2P_alt == e_P2Q_alt)

    # Test negation properties
    println("\nNegation properties with negative u:")
    e_negP = pairing_negative_u(-P, Q)
    e_negQ = pairing_negative_u(P, -Q)
    println("  e(-P, Q) == inv(e(P, Q)): ", e_negP == inv(e_neg))
    println("  e(P, -Q) == inv(e(P, Q)): ", e_negQ == inv(e_neg))

    return e_orig, e_neg, e_alt
end

# Check known test vectors if we have any
function check_against_known_values()
    println("\n=== Checking Against Known Values ===\n")

    P = g1_generator()
    Q = g2_generator()

    # Compute with both signs
    e_pos = GrothCurves.pairing(P, Q)
    e_neg = pairing_negative_u(P, Q)
    e_neg_alt = pairing_negative_u_alt(P, Q)

    println("Results are different:")
    println("  Positive u result == Negative u result? ", e_pos == e_neg)
    println("  Positive u result == Alt negative result? ", e_pos == e_neg_alt)
    println("  Negative u result == Alt negative result? ", e_neg == e_neg_alt)

    # The actual pairing result should be deterministic
    # If we had a known correct value, we could check against it
    println("\nNote: Without known test vectors, we can't determine which is correct")
    println("But bilinearity should hold for the correct implementation")
end

# Main execution
println("Starting tests...\n")

# First verify which sign might be correct
Q_u_pos, Q_u_neg, Q_pi, Q_pi2 = verify_u_sign()

# Test bilinearity with both signs
e_orig, e_neg, e_alt = test_bilinearity_negative()

# Check against known values
check_against_known_values()

println("\n=== Conclusion ===")
println("If bilinearity works with negative u but not positive u,")
println("then we've been using the wrong sign for the BN254 parameter!")
println("Many BN254 implementations use u = -4965661367192848881")
