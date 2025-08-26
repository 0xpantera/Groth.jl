using GrothCurves

# Debug Miller loop to see what's happening
function miller_loop_debug(P::G1Point, Q::G2Point, label::String)
    if iszero(P) || iszero(Q)
        return one(Fp12Element)
    end
    
    println("\n=== Miller loop for $label ===")
    
    T = Q
    f = one(Fp12Element)
    
    # Use current ATE_LOOP_COUNT (which is u)
    loop_bits = digits(Bool, GrothCurves.ATE_LOOP_COUNT, base=2)
    
    bit_count = 0
    for i in (length(loop_bits)-1):-1:1
        bit_count += 1
        
        # Double step
        T_new, line_coeffs = GrothCurves.doubling_step(T)
        T = T_new
        line_val = GrothCurves.evaluate_line(line_coeffs, P)
        f = f^2 * line_val
        
        # Addition step if bit is 1
        if loop_bits[i]
            T_new, line_coeffs = GrothCurves.addition_step(T, Q)
            T = T_new
            line_val = GrothCurves.evaluate_line(line_coeffs, P)
            f = f * line_val
        end
    end
    
    println("After main loop (", bit_count, " iterations):")
    println("  T = [u]Q")
    println("  f has been computed")
    
    # Frobenius corrections
    Q_pi = GrothCurves.frobenius_g2(Q, 1)
    T1, line1 = GrothCurves.addition_step(T, Q_pi)
    f = f * GrothCurves.evaluate_line(line1, P)
    println("After correction 1: Added π(Q)")
    
    Q_pi2 = GrothCurves.frobenius_g2(Q, 2)
    neg_Q_pi2 = -Q_pi2
    T2, line2 = GrothCurves.addition_step(T1, neg_Q_pi2)
    f = f * GrothCurves.evaluate_line(line2, P)
    println("After correction 2: Added -π²(Q)")
    
    return f
end

# Test with P and Q
P = g1_generator()
Q = g2_generator()

f_P_Q = miller_loop_debug(P, Q, "e(P, Q)")
f_2P_Q = miller_loop_debug(2*P, Q, "e(2P, Q)")
f_P_2Q = miller_loop_debug(P, 2*Q, "e(P, 2Q)")

# Apply final exponentiation
e_P_Q = GrothCurves.final_exponentiation(f_P_Q)
e_2P_Q = GrothCurves.final_exponentiation(f_2P_Q)
e_P_2Q = GrothCurves.final_exponentiation(f_P_2Q)

println("\n=== Final Results ===")
println("e(2P, Q) == e(P, Q)^2: ", e_2P_Q == e_P_Q^2)
println("e(P, 2Q) == e(P, Q)^2: ", e_P_2Q == e_P_Q^2)
println("e(2P, Q) == e(P, 2Q): ", e_2P_Q == e_P_2Q)

# Check if the Miller values have the right relationship before final exp
println("\n=== Miller function ratios (before final exp) ===")
ratio1 = f_2P_Q / f_P_Q^2
ratio2 = f_P_2Q / f_P_Q^2

# Just print first Fp2 component to see if there's a pattern
println("f(2P, Q) / f(P, Q)^2: first component = ", ratio1[1][1])
println("f(P, 2Q) / f(P, Q)^2: first component = ", ratio2[1][1])