using GrothCurves
using Test

# Modified Miller loop with Frobenius-adjusted line evaluations
function miller_loop_frobenius_adjusted(P::G1Point, Q::G2Point)
    if iszero(P) || iszero(Q)
        return one(Fp12Element)
    end

    T = Q
    f = one(Fp12Element)

    # Main loop over u
    loop_bits = digits(Bool, GrothCurves.ATE_LOOP_COUNT, base=2)

    for i in (length(loop_bits)-1):-1:1
        # Double step
        T_new, line_coeffs = GrothCurves.doubling_step(T)
        T = T_new
        f = f^2 * GrothCurves.evaluate_line(line_coeffs, P)

        # Addition step if bit is 1
        if loop_bits[i]
            T_new, line_coeffs = GrothCurves.addition_step(T, Q)
            T = T_new
            f = f * GrothCurves.evaluate_line(line_coeffs, P)
        end
    end

    # Modified correction steps with Frobenius on line evaluations
    Q_pi = GrothCurves.frobenius_g2(Q, 1)
    Q_pi2 = GrothCurves.frobenius_g2(Q, 2)

    # Correction step 1: Apply Frobenius to the line evaluation
    T1, line1 = GrothCurves.addition_step(T, Q_pi)
    line1_eval = GrothCurves.evaluate_line(line1, P)
    # Apply p-power Frobenius to the line evaluation
    f = f * GrothCurves.frobenius_map(line1_eval, 1)

    # Correction step 2: Apply Frobenius^2 to the line evaluation
    neg_Q_pi2 = -Q_pi2
    T2, line2 = GrothCurves.addition_step(T1, neg_Q_pi2)
    line2_eval = GrothCurves.evaluate_line(line2, P)
    # Apply p^2-power Frobenius to the line evaluation
    f = f * GrothCurves.frobenius_map(line2_eval, 2)

    return f
end

# Alternative: Try without Frobenius on correction lines
function miller_loop_no_frobenius_on_corrections(P::G1Point, Q::G2Point)
    if iszero(P) || iszero(Q)
        return one(Fp12Element)
    end

    T = Q
    f = one(Fp12Element)

    loop_bits = digits(Bool, GrothCurves.ATE_LOOP_COUNT, base=2)

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

    # Correction steps WITHOUT Frobenius on line evaluations (control)
    Q_pi = GrothCurves.frobenius_g2(Q, 1)
    Q_pi2 = GrothCurves.frobenius_g2(Q, 2)

    T1, line1 = GrothCurves.addition_step(T, Q_pi)
    f = f * GrothCurves.evaluate_line(line1, P)

    neg_Q_pi2 = -Q_pi2
    T2, line2 = GrothCurves.addition_step(T1, neg_Q_pi2)
    f = f * GrothCurves.evaluate_line(line2, P)

    return f
end

# Alternative: Try conjugating the line evaluations instead
function miller_loop_conjugate_corrections(P::G1Point, Q::G2Point)
    if iszero(P) || iszero(Q)
        return one(Fp12Element)
    end

    T = Q
    f = one(Fp12Element)

    loop_bits = digits(Bool, GrothCurves.ATE_LOOP_COUNT, base=2)

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

    # Correction steps with conjugated line evaluations
    Q_pi = GrothCurves.frobenius_g2(Q, 1)
    Q_pi2 = GrothCurves.frobenius_g2(Q, 2)

    T1, line1 = GrothCurves.addition_step(T, Q_pi)
    line1_eval = GrothCurves.evaluate_line(line1, P)
    # Try conjugating instead of Frobenius
    f = f * GrothCurves.conjugate(line1_eval)

    neg_Q_pi2 = -Q_pi2
    T2, line2 = GrothCurves.addition_step(T1, neg_Q_pi2)
    line2_eval = GrothCurves.evaluate_line(line2, P)
    # Conjugate this one too
    f = f * GrothCurves.conjugate(line2_eval)

    return f
end

# Test pairing functions
function pairing_frobenius_adjusted(P::G1Point, Q::G2Point)
    f = miller_loop_frobenius_adjusted(P, Q)
    return GrothCurves.final_exponentiation(f)
end

function pairing_no_frobenius(P::G1Point, Q::G2Point)
    f = miller_loop_no_frobenius_on_corrections(P, Q)
    return GrothCurves.final_exponentiation(f)
end

function pairing_conjugate(P::G1Point, Q::G2Point)
    f = miller_loop_conjugate_corrections(P, Q)
    return GrothCurves.final_exponentiation(f)
end

# Verify that T2 = [6u+2]Q after corrections
function verify_final_point(Q::G2Point)
    T = Q

    # Main loop: T = [u]Q
    loop_bits = digits(Bool, GrothCurves.ATE_LOOP_COUNT, base=2)

    for i in (length(loop_bits)-1):-1:1
        T_new, _ = GrothCurves.doubling_step(T)
        T = T_new

        if loop_bits[i]
            T_new, _ = GrothCurves.addition_step(T, Q)
            T = T_new
        end
    end

    println("After main loop: T = [u]Q")

    # Apply corrections
    Q_pi = GrothCurves.frobenius_g2(Q, 1)
    T1, _ = GrothCurves.addition_step(T, Q_pi)

    println("After correction 1: T1 = [u]Q + π(Q)")

    Q_pi2 = GrothCurves.frobenius_g2(Q, 2)
    neg_Q_pi2 = -Q_pi2
    T2, _ = GrothCurves.addition_step(T1, neg_Q_pi2)

    println("After correction 2: T2 = [u]Q + π(Q) - π²(Q)")

    # Check if T2 = [6u+2]Q
    u = BigInt(4965661367192848881)
    expected = (6 * u + 2) * Q

    println("T2 == [6u+2]Q? ", T2 == expected)

    return T2, expected
end

# Main test
println("=== Testing Modified Miller Loop with Frobenius Adjustments ===\n")

P = g1_generator()
Q = g2_generator()

# First verify the final point
println("Verifying final accumulator point:")
T2, expected = verify_final_point(Q)
println()

# Test original implementation (baseline)
println("Testing ORIGINAL implementation:")
e_orig = pairing(P, Q)
e_2P_orig = pairing(2 * P, Q)
e_P2Q_orig = pairing(P, 2 * Q)
e_squared_orig = e_orig^2

println("  e(2P, Q) == e(P, Q)^2: ", e_2P_orig == e_squared_orig)
println("  e(P, 2Q) == e(P, Q)^2: ", e_P2Q_orig == e_squared_orig)
println("  e(2P, Q) == e(P, 2Q): ", e_2P_orig == e_P2Q_orig)
println()

# Test with Frobenius-adjusted corrections
println("Testing with FROBENIUS-ADJUSTED line evaluations:")
e_frob = pairing_frobenius_adjusted(P, Q)
e_2P_frob = pairing_frobenius_adjusted(2 * P, Q)
e_P2Q_frob = pairing_frobenius_adjusted(P, 2 * Q)
e_squared_frob = e_frob^2

println("  e(2P, Q) == e(P, Q)^2: ", e_2P_frob == e_squared_frob)
println("  e(P, 2Q) == e(P, Q)^2: ", e_P2Q_frob == e_squared_frob)
println("  e(2P, Q) == e(P, 2Q): ", e_2P_frob == e_P2Q_frob)
println()

# Test without Frobenius (control)
println("Testing WITHOUT Frobenius on corrections (control):")
e_no_frob = pairing_no_frobenius(P, Q)
e_2P_no_frob = pairing_no_frobenius(2 * P, Q)
e_P2Q_no_frob = pairing_no_frobenius(P, 2 * Q)
e_squared_no_frob = e_no_frob^2

println("  e(2P, Q) == e(P, Q)^2: ", e_2P_no_frob == e_squared_no_frob)
println("  e(P, 2Q) == e(P, Q)^2: ", e_P2Q_no_frob == e_squared_no_frob)
println("  e(2P, Q) == e(P, 2Q): ", e_2P_no_frob == e_P2Q_no_frob)
println()

# Test with conjugated corrections
println("Testing with CONJUGATED line evaluations:")
e_conj = pairing_conjugate(P, Q)
e_2P_conj = pairing_conjugate(2 * P, Q)
e_P2Q_conj = pairing_conjugate(P, 2 * Q)
e_squared_conj = e_conj^2

println("  e(2P, Q) == e(P, Q)^2: ", e_2P_conj == e_squared_conj)
println("  e(P, 2Q) == e(P, Q)^2: ", e_P2Q_conj == e_squared_conj)
println("  e(2P, Q) == e(P, 2Q): ", e_2P_conj == e_P2Q_conj)
println()

# Compare results
println("=== Comparing different approaches ===")
println("Original == No Frobenius? ", e_orig == e_no_frob, " (should be true)")
println("Original == Frobenius-adjusted? ", e_orig == e_frob)
println("Original == Conjugated? ", e_orig == e_conj)

# Test negation properties with each approach
println("\n=== Testing negation properties ===")
for (name, pairing_fn) in [
    ("Original", pairing),
    ("Frobenius-adjusted", pairing_frobenius_adjusted),
    ("Conjugated", pairing_conjugate)
]
    e = pairing_fn(P, Q)
    e_negP = pairing_fn(-P, Q)
    e_negQ = pairing_fn(P, -Q)

    println("$name:")
    println("  e(-P, Q) == inv(e(P, Q)): ", e_negP == inv(e))
    println("  e(P, -Q) == inv(e(P, Q)): ", e_negQ == inv(e))
end
