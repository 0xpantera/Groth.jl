using GrothCurves
using Test

# BN254 parameters
const u = BigInt(4965661367192848881)
const p = GrothCurves.BN254_PRIME
const r = parse(BigInt, "21888242871839275222246405745257275088548364400416034343698204186575808495617")

println("=== Understanding BN254 Optimal Ate Pairing Formula ===\n")
println("BN254 Parameters:")
println("  u = ", u)
println("  p = ", p)
println("  r (group order) = ", r)
println("  t (trace) = 6u + 2 = ", 6 * u + 2)
println()

# Test the Frobenius endomorphism eigenvalue equation
function test_frobenius_eigenvalue(Q::G2Point)
    println("Testing Frobenius eigenvalue equation on G2:")

    # For BN254, the characteristic polynomial of Frobenius on E'(Fp2) is:
    # X^2 - tX + p = 0, where t is the trace
    # For points in the r-torsion, we have a simpler relation

    Q_pi = GrothCurves.frobenius_g2(Q, 1)   # π(Q)
    Q_pi2 = GrothCurves.frobenius_g2(Q, 2)  # π²(Q)

    # Test: π²(Q) - t*π(Q) + p*Q = O ?
    t = 6 * u + 2
    lhs = Q_pi2 - t * Q_pi + p * Q
    println("  π²(Q) - t*π(Q) + p*Q == O? ", iszero(lhs))

    # Actually, for the r-torsion, we should have:
    # π²(Q) + Q = t*π(Q) (mod r)
    # Let's test this
    lhs2 = Q_pi2 + Q
    rhs2 = t * Q_pi
    println("  π²(Q) + Q == t*π(Q)? ", lhs2 == rhs2)

    # Another relation: π(Q) = p*Q in E'[r]
    # But p mod r gives us the eigenvalue
    p_mod_r = p % r
    println("  p mod r = ", p_mod_r)
    test_pQ = p_mod_r * Q
    println("  π(Q) == (p mod r)*Q? ", Q_pi == test_pQ)

    return Q_pi, Q_pi2
end

# Test different formulas for the ate pairing
function test_ate_formulas(P::G1Point, Q::G2Point)
    println("\n=== Testing Different Ate Pairing Formulas ===\n")

    # Formula 1: Just f_{u,Q}(P)
    println("Formula 1: f_{u,Q}(P)")
    f1 = compute_miller_u(P, Q)
    e1 = GrothCurves.final_exponentiation(f1)
    test_bilinearity("Formula 1", e1, P, Q)

    # Formula 2: f_{6u+2,Q}(P)
    println("\nFormula 2: f_{6u+2,Q}(P)")
    f2 = compute_miller_6u_plus_2(P, Q)
    e2 = GrothCurves.final_exponentiation(f2)
    test_bilinearity("Formula 2", e2, P, Q)

    # Formula 3: f_{u,Q}(P) * line_{[u]Q,π(Q)}(P) * line_{[u]Q+π(Q),-π²(Q)}(P)
    println("\nFormula 3: Current implementation")
    f3 = GrothCurves.miller_loop(P, Q)
    e3 = GrothCurves.final_exponentiation(f3)
    test_bilinearity("Formula 3", e3, P, Q)

    # Formula 4: Try a different correction
    # Based on the trace formula: [t]Q = [6u+2]Q = π²(Q) + Q
    # So [6u+2]Q - [u]Q = [6u+2-u]Q = [5u+2]Q
    println("\nFormula 4: Different correction approach")
    f4 = compute_miller_with_alt_correction(P, Q)
    e4 = GrothCurves.final_exponentiation(f4)
    test_bilinearity("Formula 4", e4, P, Q)

    return e1, e2, e3, e4
end

function compute_miller_u(P::G1Point, Q::G2Point)
    if iszero(P) || iszero(Q)
        return one(Fp12Element)
    end

    T = Q
    f = one(Fp12Element)

    loop_bits = digits(Bool, u, base=2)

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

    return f
end

function compute_miller_6u_plus_2(P::G1Point, Q::G2Point)
    if iszero(P) || iszero(Q)
        return one(Fp12Element)
    end

    T = Q
    f = one(Fp12Element)

    loop_count = 6 * u + 2
    loop_bits = digits(Bool, loop_count, base=2)

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

    return f
end

function compute_miller_with_alt_correction(P::G1Point, Q::G2Point)
    if iszero(P) || iszero(Q)
        return one(Fp12Element)
    end

    # Start with f_{u,Q}(P)
    T = Q
    f = one(Fp12Element)

    loop_bits = digits(Bool, u, base=2)

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

    # Now T = [u]Q
    # We need to get to [6u+2]Q
    # [6u+2]Q - [u]Q = [5u+2]Q
    # Let's try a different approach: use the fact that π²(Q) + Q = [t]Q = [6u+2]Q

    Q_pi2 = GrothCurves.frobenius_g2(Q, 2)
    target = Q_pi2 + Q  # This should be [6u+2]Q

    # We need to add [5u+2]Q to T
    # But that's expensive. Let's try something else.

    # Alternative: The optimal ate might be f_{u,π(Q)}(P) instead
    # Let's just try adding what the trace formula gives us
    diff = target - T  # This is [6u+2]Q - [u]Q = [5u+2]Q

    if !iszero(diff)
        # Add a line from T to target
        T_final, line_final = GrothCurves.addition_step(T, diff)
        f = f * GrothCurves.evaluate_line(line_final, P)
    end

    return f
end

function test_bilinearity(formula_name::String, e, P::G1Point, Q::G2Point)
    e_2P = GrothCurves.pairing(2 * P, Q)
    e_P2Q = GrothCurves.pairing(P, 2 * Q)
    e_squared = e^2

    # Use the same pairing function for consistency
    # Actually, let's compute them fresh with the same approach

    println("  Testing bilinearity for $formula_name:")
    println("    e(2P, Q) == e(P, Q)^2: ", e_2P == e_squared, " (using standard pairing for 2P)")
    println("    Note: This test isn't fully consistent - need same formula for all")
end

# Investigate the final accumulator point
function analyze_accumulator_points(Q::G2Point)
    println("\n=== Analyzing Accumulator Points ===\n")

    # Compute various scalar multiples of Q
    Q_u = u * Q
    Q_6u = 6 * u * Q
    Q_6u_plus_2 = (6 * u + 2) * Q

    # Compute Frobenius images
    Q_pi = GrothCurves.frobenius_g2(Q, 1)
    Q_pi2 = GrothCurves.frobenius_g2(Q, 2)

    println("Checking relationships:")
    println("  [u]Q computed correctly? (verify independently)")

    # What does [u]Q + π(Q) - π²(Q) equal?
    result = Q_u + Q_pi - Q_pi2
    println("  [u]Q + π(Q) - π²(Q) == [6u+2]Q? ", result == Q_6u_plus_2)
    println("  [u]Q + π(Q) - π²(Q) == O? ", iszero(result))

    # Check if π²(Q) + Q = [6u+2]Q (trace formula)
    trace_result = Q_pi2 + Q
    println("  π²(Q) + Q == [6u+2]Q? ", trace_result == Q_6u_plus_2)

    # What about π(Q) relationship?
    # For BN curves, we have the 6-th twist
    # The eigenvalue of Frobenius should be related to p

    println("\nFrobenius eigenvalues:")
    # The minimal polynomial of π on G2[r] should give us hints

    # Try to find λ such that π(Q) = λ*Q
    # We know π is not a scalar multiplication, but let's see patterns

    return Q_u, Q_6u_plus_2, Q_pi, Q_pi2
end

# Main testing
println("Starting analysis...\n")

P = g1_generator()
Q = g2_generator()

# Test Frobenius eigenvalue equation
Q_pi, Q_pi2 = test_frobenius_eigenvalue(Q)

# Analyze accumulator points
Q_u, Q_6u_plus_2, _, _ = analyze_accumulator_points(Q)

# Test different ate formulas
e1, e2, e3, e4 = test_ate_formulas(P, Q)

println("\n=== Summary ===")
println("The issue is that the correction formula doesn't produce [6u+2]Q")
println("This suggests we need a different approach to the optimal ate pairing")
println("Possible issues:")
println("1. Wrong correction formula")
println("2. Wrong loop parameter (should we use -u instead?)")
println("3. Missing some other correction factor")
