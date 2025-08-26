using GrothCurves

# Miller loop using FULL Fp12 multiplication - no sparse tricks at all
function miller_loop_full_fp12(P::G1Point, Q::G2Point)
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
        
        # Convert line to FULL Fp12 (no sparse tricks)
        # Get P coordinates
        P_x, P_y = iszero(P) ? (zero(BN254Field), zero(BN254Field)) : GrothCurves.to_affine(P)
        xP = Fp2Element(P_x, zero(BN254Field))
        yP = Fp2Element(P_y, zero(BN254Field))
        
        # Build FULL Fp12 element for line evaluation
        # Using the corrected M-twist mapping: ψ: x→xv, y→yw
        # Line: a*xP*v + b*yP*w + c
        # v = (0,1,0) in Fp6, w is the generator of Fp12/Fp6
        
        # c0 = c + a*xP*v = (c, a*xP, 0) in Fp6
        c0 = Fp6Element(line_coeffs.c, line_coeffs.a * xP, zero(Fp2Element))
        
        # c1 = b*yP in first position = (b*yP, 0, 0) in Fp6
        c1 = Fp6Element(line_coeffs.b * yP, zero(Fp2Element), zero(Fp2Element))
        
        line_fp12 = Fp12Element(c0, c1)
        
        # Use FULL Fp12 multiplication (no sparse optimization)
        f = f * f * line_fp12
        
        # Addition step if bit is 1
        if loop_bits[i]
            T_new, line_coeffs = GrothCurves.addition_step(T, Q)
            T = T_new
            
            # Same full Fp12 conversion
            c0 = Fp6Element(line_coeffs.c, line_coeffs.a * xP, zero(Fp2Element))
            c1 = Fp6Element(line_coeffs.b * yP, zero(Fp2Element), zero(Fp2Element))
            line_fp12 = Fp12Element(c0, c1)
            
            f = f * line_fp12
        end
    end
    
    # Frobenius corrections
    Q_pi = GrothCurves.frobenius_g2(Q, 1)
    T1, line1 = GrothCurves.addition_step(T, Q_pi)
    
    P_x, P_y = GrothCurves.to_affine(P)
    xP = Fp2Element(P_x, zero(BN254Field))
    yP = Fp2Element(P_y, zero(BN254Field))
    c0 = Fp6Element(line1.c, line1.a * xP, zero(Fp2Element))
    c1 = Fp6Element(line1.b * yP, zero(Fp2Element), zero(Fp2Element))
    f = f * Fp12Element(c0, c1)
    
    Q_pi2 = GrothCurves.frobenius_g2(Q, 2)
    neg_Q_pi2 = -Q_pi2
    T2, line2 = GrothCurves.addition_step(T1, neg_Q_pi2)
    
    c0 = Fp6Element(line2.c, line2.a * xP, zero(Fp2Element))
    c1 = Fp6Element(line2.b * yP, zero(Fp2Element), zero(Fp2Element))
    f = f * Fp12Element(c0, c1)
    
    return f
end

function pairing_full_fp12(P::G1Point, Q::G2Point)
    f = miller_loop_full_fp12(P, Q)
    return GrothCurves.final_exponentiation(f)
end

println("=== Testing with FULL Fp12 multiplication (no sparse tricks) ===\n")

P = g1_generator()
Q = g2_generator()

e = pairing_full_fp12(P, Q)

# Test negation first
println("Negation tests:")
e_negP = pairing_full_fp12(-P, Q)
e_negQ = pairing_full_fp12(P, -Q)
println("  e(-P, Q) == inv(e(P, Q)): ", e_negP == inv(e))
println("  e(P, -Q) == inv(e(P, Q)): ", e_negQ == inv(e))

# Test bilinearity
println("\nBilinearity tests:")
e_2P = pairing_full_fp12(2*P, Q)
e_P2Q = pairing_full_fp12(P, 2*Q)
e_squared = e^2

println("  e(2P, Q) == e(P, Q)^2: ", e_2P == e_squared)
println("  e(P, 2Q) == e(P, Q)^2: ", e_P2Q == e_squared)
println("  e(2P, Q) == e(P, 2Q): ", e_2P == e_P2Q)

# Compare with standard pairing
println("\nComparing full vs standard:")
e_standard = pairing(P, Q)
println("  full == standard: ", e == e_standard)