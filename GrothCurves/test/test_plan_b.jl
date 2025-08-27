using GrothCurves
using Test

# Test Plan B: Loop over 6u+2 WITHOUT Frobenius corrections

# Temporarily override the Miller loop
function miller_loop_plan_b(P::G1Point, Q::G2Point)
    if iszero(P) || iszero(Q)
        return one(Fp12Element)
    end
    
    T = Q
    f = one(Fp12Element)
    
    # Use 6u+2 for Plan B
    loop_count = 6 * BigInt(4965661367192848881) + 2
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
    
    # NO FROBENIUS CORRECTIONS for Plan B!
    return f
end

# Test with Plan B
P = g1_generator()
Q = g2_generator()

# Compute pairing with Plan B
f_miller_b = miller_loop_plan_b(P, Q)
e_plan_b = GrothCurves.final_exponentiation(f_miller_b)

# Test bilinearity
f_2P_b = miller_loop_plan_b(2*P, Q)
e_2P_b = GrothCurves.final_exponentiation(f_2P_b)

f_2Q_b = miller_loop_plan_b(P, 2*Q)
e_2Q_b = GrothCurves.final_exponentiation(f_2Q_b)

println("Plan B (6u+2, no corrections):")
println("e(2P, Q) == e(P, Q)^2: ", e_2P_b == e_plan_b^2)
println("e(P, 2Q) == e(P, Q)^2: ", e_2Q_b == e_plan_b^2)
println("e(2P, Q) == e(P, 2Q): ", e_2P_b == e_2Q_b)

# Test negation
f_neg_P = miller_loop_plan_b(-P, Q)
e_neg_P = GrothCurves.final_exponentiation(f_neg_P)

f_neg_Q = miller_loop_plan_b(P, -Q)
e_neg_Q = GrothCurves.final_exponentiation(f_neg_Q)

println("\nNegation:")
println("e(-P, Q) == e(P, Q)^(-1): ", e_neg_P == inv(e_plan_b))
println("e(P, -Q) == e(P, Q)^(-1): ", e_neg_Q == inv(e_plan_b))

if e_2P_b == e_plan_b^2 && e_2Q_b == e_plan_b^2
    println("\n✅ Plan B works! The issue was mixing formulas!")
else
    println("\n❌ Plan B also fails. Something deeper is wrong.")
end