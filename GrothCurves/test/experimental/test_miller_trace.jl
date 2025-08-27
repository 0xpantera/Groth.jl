# Detailed trace through the Miller loop to debug bilinearity issues
using GrothCurves
using GrothAlgebra

println("="^80)
println("Miller Loop Trace Debugging")
println("="^80)

# Test setup
P = g1_generator()
Q = g2_generator()

# We'll test e([2]P, Q) vs e(P, [2]Q) since these should be equal for bilinearity
P2 = P + P
Q2 = Q + Q

println("\nTest points:")
println("P = G1 generator")
println("Q = G2 generator")
println("[2]P computed")
println("[2]Q computed")

# Custom Miller loop with detailed tracing
function miller_loop_traced(P::G1Point, Q::G2Point, label::String)
    println("\n" * "="^60)
    println("Tracing Miller loop for $label")
    println("="^60)

    if iszero(P) || iszero(Q)
        println("Zero point detected, returning 1")
        return one(Fp12Element)
    end

    # Initialize
    T = Q
    f = one(Fp12Element)

    # Get binary representation
    loop_bits = digits(Bool, GrothCurves.ATE_LOOP_COUNT, base=2)
    println("Loop count (u) = ", GrothCurves.ATE_LOOP_COUNT)
    println("Number of bits = ", length(loop_bits))

    # Count actual iterations
    num_doubles = length(loop_bits) - 1
    num_adds = sum(loop_bits[1:end-1])
    println("Expected doublings: $num_doubles")
    println("Expected additions: $num_adds")

    actual_doubles = 0
    actual_adds = 0

    # Main Miller loop
    for i in (length(loop_bits)-1):-1:1
        # Double step
        T_new, line_coeffs = GrothCurves.doubling_step(T)
        T = T_new
        line_eval = GrothCurves.evaluate_line(line_coeffs, P)
        f = f^2 * line_eval
        actual_doubles += 1

        if i == length(loop_bits) - 1
            # First iteration details
            println("\nFirst iteration (i=$i):")
            println("  After doubling:")
            x_T, y_T = GrothCurves.to_affine(T)
            println("    T.x first component = ", x_T[1])
            println("    T.y first component = ", y_T[1])
            println("    Line evaluation = ", line_eval)
            println("    f accumulator first component = ", f[1][1][1])
        end

        # Addition step if bit is 1
        if loop_bits[i]
            T_new, line_coeffs = GrothCurves.addition_step(T, Q)
            T = T_new
            line_eval = GrothCurves.evaluate_line(line_coeffs, P)
            f = f * line_eval
            actual_adds += 1

            if actual_adds == 1
                println("\nFirst addition (i=$i):")
                x_T, y_T = GrothCurves.to_affine(T)
                println("    T.x first component = ", x_T[1])
                println("    T.y first component = ", y_T[1])
                println("    Line evaluation = ", line_eval)
            end
        end
    end

    println("\nAfter main loop:")
    println("  Actual doublings: $actual_doubles")
    println("  Actual additions: $actual_adds")
    x_T, y_T = GrothCurves.to_affine(T)
    println("  T = [u]Q:")
    println("    x first component = ", x_T[1])
    println("    y first component = ", y_T[1])
    println("  f accumulator first component = ", f[1][1][1])

    # Correction steps
    println("\nCorrection steps:")

    # Step 1: Line through T and π(Q)
    Q_pi = GrothCurves.frobenius_g2(Q, 1)
    T1, line1 = GrothCurves.addition_step(T, Q_pi)
    line1_eval = GrothCurves.evaluate_line(line1, P)
    f = f * line1_eval

    x_Qpi, y_Qpi = GrothCurves.to_affine(Q_pi)
    println("  π(Q):")
    println("    x first component = ", x_Qpi[1])
    println("    y first component = ", y_Qpi[1])
    println("  After adding π(Q) to T:")
    x_T1, y_T1 = GrothCurves.to_affine(T1)
    println("    T1.x first component = ", x_T1[1])
    println("    T1.y first component = ", y_T1[1])
    println("    Line evaluation = ", line1_eval)

    # Step 2: Line through T1 and -π²(Q)
    Q_pi2 = GrothCurves.frobenius_g2(Q, 2)
    neg_Q_pi2 = -Q_pi2
    T2, line2 = GrothCurves.addition_step(T1, neg_Q_pi2)
    line2_eval = GrothCurves.evaluate_line(line2, P)
    f = f * line2_eval

    x_Qpi2, y_Qpi2 = GrothCurves.to_affine(Q_pi2)
    x_neg_Qpi2, y_neg_Qpi2 = GrothCurves.to_affine(neg_Q_pi2)
    println("  π²(Q):")
    println("    x first component = ", x_Qpi2[1])
    println("    y first component = ", y_Qpi2[1])
    println("  -π²(Q):")
    println("    x first component = ", x_neg_Qpi2[1])
    println("    y first component = ", y_neg_Qpi2[1])
    println("  After adding -π²(Q) to T1:")
    x_T2, y_T2 = GrothCurves.to_affine(T2)
    println("    T2.x first component = ", x_T2[1])
    println("    T2.y first component = ", y_T2[1])
    println("    Line evaluation = ", line2_eval)

    println("\nFinal f (before final exp):")
    println("  First component = ", f[1][1][1])

    return f
end

# Run traced Miller loops
println("\n" * "="^80)
println("COMPARING MILLER LOOPS")
println("="^80)

f1 = miller_loop_traced(P2, Q, "e([2]P, Q)")
f2 = miller_loop_traced(P, Q2, "e(P, [2]Q)")

println("\n" * "="^60)
println("Comparison of results (before final exponentiation):")
println("="^60)
println("\nf1 (from e([2]P, Q)):")
println("  First component = ", f1[1][1][1])
println("\nf2 (from e(P, [2]Q)):")
println("  First component = ", f2[1][1][1])
println("\nAre they equal? ", f1 == f2)

if f1 != f2
    println("\nDifference detected in Miller loop output!")
    println("This suggests the issue is in the Miller loop itself, not final exponentiation.")

    # Check if they're related by some simple transformation
    println("\nChecking relationships:")
    println("  f1 / f2 first component = ", (f1/f2)[1][1][1])
    println("  f2 / f1 first component = ", (f2/f1)[1][1][1])

    # Check if conjugate relationship
    f1_conj = GrothCurves.conjugate(f1)
    println("  f1 == conjugate(f2)? ", f1 == GrothCurves.conjugate(f2))
    println("  f2 == conjugate(f1)? ", f2 == f1_conj)
end

# Now apply final exponentiation
println("\n" * "="^60)
println("After final exponentiation:")
println("="^60)

e1 = GrothCurves.final_exponentiation(f1)
e2 = GrothCurves.final_exponentiation(f2)

println("\ne1 = final_exp(f1):")
println("  First component = ", e1[1][1][1])
println("\ne2 = final_exp(f2):")
println("  First component = ", e2[1][1][1])
println("\nAre they equal? ", e1 == e2)

if e1 != e2
    println("\nBilinearity FAILED!")
    println("The pairing is not bilinear, indicating a bug in the implementation.")
else
    println("\nBilinearity SUCCESS!")
    println("e([2]P, Q) = e(P, [2]Q)")
end

# Additional check: compare with direct pairing computation
println("\n" * "="^60)
println("Cross-check with direct pairing function:")
println("="^60)

e1_direct = pairing(P2, Q)
e2_direct = pairing(P, Q2)

println("pairing([2]P, Q) == e1? ", e1_direct == e1)
println("pairing(P, [2]Q) == e2? ", e2_direct == e2)

if e1_direct != e1 || e2_direct != e2
    println("\nWARNING: Direct pairing computation differs from traced version!")
    println("This suggests an issue with the tracing code itself.")
end
