using GrothCurves

# Plan B: Loop over 6u+2 without Frobenius corrections
function miller_loop_plan_b(P::G1Point, Q::G2Point)
    if iszero(P) || iszero(Q)
        return one(Fp12Element)
    end
    
    T = Q
    f = one(Fp12Element)
    
    # Use 6u + 2 instead of just u
    loop_count = 6 * GrothCurves.BN254_U + 2
    loop_bits = digits(Bool, loop_count, base=2)
    
    println("Plan B: Looping over ", loop_count, " (6u+2)")
    println("Number of iterations: ", length(loop_bits)-1)
    
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
    
    # NO Frobenius corrections for Plan B!
    
    return f
end

function pairing_plan_b(P::G1Point, Q::G2Point)
    f = miller_loop_plan_b(P, Q)
    return GrothCurves.final_exponentiation(f)
end

println("=== Testing Plan B: 6u+2 iterations, NO Frobenius corrections ===\n")

P = g1_generator()
Q = g2_generator()

e = pairing_plan_b(P, Q)

# Test negation
println("\nNegation tests:")
e_negP = pairing_plan_b(-P, Q)
e_negQ = pairing_plan_b(P, -Q)
println("  e(-P, Q) == inv(e(P, Q)): ", e_negP == inv(e))
println("  e(P, -Q) == inv(e(P, Q)): ", e_negQ == inv(e))

# Test bilinearity
println("\nBilinearity tests:")
e_2P = pairing_plan_b(2*P, Q)
e_P2Q = pairing_plan_b(P, 2*Q)
e_squared = e^2

println("  e(2P, Q) == e(P, Q)^2: ", e_2P == e_squared)
println("  e(P, 2Q) == e(P, Q)^2: ", e_P2Q == e_squared)
println("  e(2P, Q) == e(P, 2Q): ", e_2P == e_P2Q)

# Test identity
println("\nIdentity tests:")
O1 = zero(G1Point)
O2 = zero(G2Point)
println("  e(O, Q) == 1: ", isone(pairing_plan_b(O1, Q)))
println("  e(P, O) == 1: ", isone(pairing_plan_b(P, O2)))