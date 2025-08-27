using GrothCurves

# Vertical line evaluation (Z-free to avoid inversions)
function vertical_line_at(P::G1Point, Xp::Fp2Element, Zp::Fp2Element)
    # Vertical line through T' = (X':Y':Z') evaluated at P
    # v(x,y) = x - X'/Z'^2, but we use Z-free form: Z'^2 * x_P - X'
    
    Px, _ = iszero(P) ? (zero(BN254Field), zero(BN254Field)) : GrothCurves.to_affine(P)
    xP = Fp2Element(Px, zero(BN254Field))  # lift to Fp2
    Z2 = Zp^2
    
    # Compute v = Z'^2 * x_P - X'
    v = Z2 * xP - Xp
    
    # Embed in Fp12 at the c0.c0 position (same slot as b*yP term)
    # This is Fp12(Fp6(v, 0, 0), 0)
    v_fp6 = Fp6Element(v, zero(Fp2Element), zero(Fp2Element))
    v_fp12 = Fp12Element(v_fp6, zero(Fp6Element))
    
    return v_fp12
end

# Modified Miller loop WITH vertical lines
function miller_loop_with_verticals(P::G1Point, Q::G2Point)
    if iszero(P) || iszero(Q)
        return one(Fp12Element)
    end
    
    T = Q
    f = one(Fp12Element)
    
    loop_bits = digits(Bool, GrothCurves.ATE_LOOP_COUNT, base=2)
    
    for i in (length(loop_bits)-1):-1:1
        # Double step
        T_new, line_coeffs = GrothCurves.doubling_step(T)
        f = f^2 * GrothCurves.evaluate_line(line_coeffs, P)
        
        # Divide by vertical at 2T
        v = vertical_line_at(P, GrothCurves.x_coord(T_new), GrothCurves.z_coord(T_new))
        f = f / v
        
        T = T_new
        
        # Addition step if bit is 1
        if loop_bits[i]
            T_new, line_coeffs = GrothCurves.addition_step(T, Q)
            f = f * GrothCurves.evaluate_line(line_coeffs, P)
            
            # Divide by vertical at T+Q
            v = vertical_line_at(P, GrothCurves.x_coord(T_new), GrothCurves.z_coord(T_new))
            f = f / v
            
            T = T_new
        end
    end
    
    # Frobenius corrections (same as before)
    Q_pi = GrothCurves.frobenius_g2(Q, 1)
    T1, line1 = GrothCurves.addition_step(T, Q_pi)
    f = f * GrothCurves.evaluate_line(line1, P)
    v1 = vertical_line_at(P, GrothCurves.x_coord(T1), GrothCurves.z_coord(T1))
    f = f / v1
    
    Q_pi2 = GrothCurves.frobenius_g2(Q, 2)
    neg_Q_pi2 = -Q_pi2
    T2, line2 = GrothCurves.addition_step(T1, neg_Q_pi2)
    f = f * GrothCurves.evaluate_line(line2, P)
    v2 = vertical_line_at(P, GrothCurves.x_coord(T2), GrothCurves.z_coord(T2))
    f = f / v2
    
    return f
end

# Test with verticals
P = g1_generator()
Q = g2_generator()

println("Testing WITH vertical lines (denominator NOT eliminated):")
f_with_v = miller_loop_with_verticals(P, Q)
e_with_v = GrothCurves.final_exponentiation(f_with_v)

f_2P_with_v = miller_loop_with_verticals(2*P, Q)
e_2P_with_v = GrothCurves.final_exponentiation(f_2P_with_v)

f_P2Q_with_v = miller_loop_with_verticals(P, 2*Q)
e_P2Q_with_v = GrothCurves.final_exponentiation(f_P2Q_with_v)

println("e(2P, Q) == e(P, Q)^2: ", e_2P_with_v == e_with_v^2)
println("e(P, 2Q) == e(P, Q)^2: ", e_P2Q_with_v == e_with_v^2)
println("e(2P, Q) == e(P, 2Q): ", e_2P_with_v == e_P2Q_with_v)

println("\nNegation tests:")
f_negP_with_v = miller_loop_with_verticals(-P, Q)
e_negP_with_v = GrothCurves.final_exponentiation(f_negP_with_v)
println("e(-P, Q) == inv(e(P, Q)): ", e_negP_with_v == inv(e_with_v))

f_negQ_with_v = miller_loop_with_verticals(P, -Q)
e_negQ_with_v = GrothCurves.final_exponentiation(f_negQ_with_v)
println("e(P, -Q) == inv(e(P, Q)): ", e_negQ_with_v == inv(e_with_v))