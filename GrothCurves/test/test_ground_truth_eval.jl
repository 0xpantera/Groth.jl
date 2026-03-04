using GrothCurves

# Ground-truth line evaluation without sparse tricks
function evaluate_line_ground_truth(coeffs::GrothCurves.LineCoeffs, P::G1Point)
    # Affine P over Fp, lift to Fp2 (pure real)
    P_x, P_y = iszero(P) ? (zero(BN254Fq), zero(BN254Fq)) : GrothCurves.to_affine(P)
    xP = Fp2Element(P_x, zero(BN254Fq))
    yP = Fp2Element(P_y, zero(BN254Fq))

    # v = (0, 1, 0) in Fp6
    v_fp6 = Fp6Element(zero(Fp2Element), one(Fp2Element), zero(Fp2Element))

    # Lift a,b,c ∈ Fp2 to Fp6 in c0-slot
    a6 = Fp6Element(coeffs.a, zero(Fp2Element), zero(Fp2Element))
    b6 = Fp6Element(coeffs.b, zero(Fp2Element), zero(Fp2Element))
    c6 = Fp6Element(coeffs.c, zero(Fp2Element), zero(Fp2Element))

    # ψ(P): x̃ = xP * v (in c0),  ỹ = yP * v (will sit with w)
    axv = a6 * Fp6Element(xP, zero(Fp2Element), zero(Fp2Element)) * v_fp6
    byv = b6 * Fp6Element(yP, zero(Fp2Element), zero(Fp2Element)) * v_fp6

    c0 = axv + c6         # (a*xP·v) + c
    c1 = byv              # (b*yP·v) multiplies the w limb

    return Fp12Element(c0, c1)
end

# Miller loop using ground-truth evaluation
function miller_loop_ground_truth(P::G1Point, Q::G2Point)
    if iszero(P) || iszero(Q)
        return one(Fp12Element)
    end
    
    T = Q
    f = one(Fp12Element)
    
    loop_bits = digits(Bool, GrothCurves.ATE_LOOP_COUNT, base=2)
    
    for i in (length(loop_bits)-1):-1:1
        # Double step
        T_new, line_coeffs = GrothCurves.doubling_step(T)
        T = T_new
        f = f^2 * evaluate_line_ground_truth(line_coeffs, P)
        
        # Addition step if bit is 1
        if loop_bits[i]
            T_new, line_coeffs = GrothCurves.addition_step(T, Q)
            T = T_new
            f = f * evaluate_line_ground_truth(line_coeffs, P)
        end
    end
    
    # Frobenius corrections (same as original)
    Q_pi = GrothCurves.frobenius_g2(Q, 1)
    T1, line1 = GrothCurves.addition_step(T, Q_pi)
    f = f * evaluate_line_ground_truth(line1, P)
    
    Q_pi2 = GrothCurves.frobenius_g2(Q, 2)
    neg_Q_pi2 = -Q_pi2
    T2, line2 = GrothCurves.addition_step(T1, neg_Q_pi2)
    f = f * evaluate_line_ground_truth(line2, P)
    
    return f
end

# Test pairing with ground-truth evaluation
function pairing_ground_truth(P::G1Point, Q::G2Point)
    f = miller_loop_ground_truth(P, Q)
    return GrothCurves.final_exponentiation(f)
end

println("=== Testing with GROUND-TRUTH line evaluation (no sparse tricks) ===\n")

P = g1_generator()
Q = g2_generator()

e = pairing_ground_truth(P, Q)

# Test negation first (should still work)
println("Negation tests:")
e_negP = pairing_ground_truth(-P, Q)
e_negQ = pairing_ground_truth(P, -Q)
println("  e(-P, Q) == inv(e(P, Q)): ", e_negP == inv(e))
println("  e(P, -Q) == inv(e(P, Q)): ", e_negQ == inv(e))

# Test bilinearity
println("\nBilinearity tests:")
e_2P = pairing_ground_truth(2*P, Q)
e_P2Q = pairing_ground_truth(P, 2*Q)
e_squared = e^2

println("  e(2P, Q) == e(P, Q)^2: ", e_2P == e_squared)
println("  e(P, 2Q) == e(P, Q)^2: ", e_P2Q == e_squared)
println("  e(2P, Q) == e(P, 2Q): ", e_2P == e_P2Q)

# Test identity
println("\nIdentity tests:")
O1 = zero(G1Point)
O2 = zero(G2Point)
println("  e(O, Q) == 1: ", isone(pairing_ground_truth(O1, Q)))
println("  e(P, O) == 1: ", isone(pairing_ground_truth(P, O2)))

# General bilinearity
println("\nGeneral bilinearity:")
e_3P_5Q = pairing_ground_truth(3*P, 5*Q)
e_15 = e^15
println("  e(3P, 5Q) == e(P, Q)^15: ", e_3P_5Q == e_15)