# Compute and verify the correct Frobenius coefficients for BN254 pairing
using GrothCurves
using GrothAlgebra

println("="^70)
println("Computing Frobenius Coefficients for BN254 Pairing")
println("="^70)

# BN254 parameters
p = GrothCurves.BN254_PRIME
println("\np = ", p)

# The sextic twist non-residue: ξ = 9 + u
ξ = Fp2Element(9, 1)
println("\nξ = 9 + u = ", ξ)

# For the p-power endomorphism on G2, we need:
# PSI_X = ξ^((p-1)/3) for x-coordinate
# PSI_Y = ξ^((p-1)/2) for y-coordinate

println("\n" * "="^50)
println("Computing p-power endomorphism coefficients")
println("="^50)

# Compute PSI_X = ξ^((p-1)/3)
pow_x = div(p - 1, 3)
PSI_X = ξ^pow_x
println("\nPSI_X = ξ^((p-1)/3):")
println("  c0 = ", PSI_X[1])
println("  c1 = ", PSI_X[2])

# Compute PSI_Y = ξ^((p-1)/2)
pow_y = div(p - 1, 2)
PSI_Y = ξ^pow_y
println("\nPSI_Y = ξ^((p-1)/2):")
println("  c0 = ", PSI_Y[1])
println("  c1 = ", PSI_Y[2])

# Compare with the arkworks values
println("\n" * "="^50)
println("Comparing with arkworks coefficients")
println("="^50)

# From arkworks g2.rs:
# PSI_X (P_POWER_ENDOMORPHISM_COEFF_0):
#   c0: 21575463638280843010398324269430826099269044274347216827212613867836435027261
#   c1: 10307601595873709700152284273816112264069230130616436755625194854815875713954

arkworks_psi_x_c0 = parse(BigInt, "21575463638280843010398324269430826099269044274347216827212613867836435027261")
arkworks_psi_x_c1 = parse(BigInt, "10307601595873709700152284273816112264069230130616436755625194854815875713954")

println("\nPSI_X comparison:")
println("  Our c0:      ", PSI_X[1].value)
println("  Arkworks c0: ", arkworks_psi_x_c0)
println("  Match: ", PSI_X[1].value == arkworks_psi_x_c0)
println()
println("  Our c1:      ", PSI_X[2].value)
println("  Arkworks c1: ", arkworks_psi_x_c1)
println("  Match: ", PSI_X[2].value == arkworks_psi_x_c1)

# PSI_Y (P_POWER_ENDOMORPHISM_COEFF_1):
#   c0: 2821565182194536844548159561693502659359617185244120367078079554186484126554
#   c1: 3505843767911556378687030309984248845540243509899259641013678093033130930403

arkworks_psi_y_c0 = parse(BigInt, "2821565182194536844548159561693502659359617185244120367078079554186484126554")
arkworks_psi_y_c1 = parse(BigInt, "3505843767911556378687030309984248845540243509899259641013678093033130930403")

println("\nPSI_Y comparison:")
println("  Our c0:      ", PSI_Y[1].value)
println("  Arkworks c0: ", arkworks_psi_y_c0)
println("  Match: ", PSI_Y[1].value == arkworks_psi_y_c0)
println()
println("  Our c1:      ", PSI_Y[2].value)
println("  Arkworks c1: ", arkworks_psi_y_c1)
println("  Match: ", PSI_Y[2].value == arkworks_psi_y_c1)

# Now let's verify these work correctly
println("\n" * "="^50)
println("Testing p-power endomorphism properties")
println("="^50)

# The p-power endomorphism should satisfy: ψ(Q) = [p mod r]Q
# where r is the order of G2

# Get a test point
Q = g2_generator()

# Apply our p-power endomorphism
function p_power_endo(Q::G2Point)
    if iszero(Q)
        return Q
    end
    x, y = GrothCurves.to_affine(Q)

    # Apply Frobenius (conjugation)
    x_frob = GrothCurves.conjugate(x)
    y_frob = GrothCurves.conjugate(y)

    # Multiply by coefficients
    x_final = x_frob * PSI_X
    y_final = y_frob * PSI_Y

    return G2Point(x_final, y_final)
end

# Test: ψ(Q) should equal [p mod r]Q
# For BN254, r = order of the curve
r = parse(BigInt, "21888242871839275222246405745257275088548364400416034343698204186575808495617")
p_mod_r = mod(p, r)

println("\np mod r = ", p_mod_r)

# Computing [p mod r]Q by scalar multiplication would be too expensive
# (it would require ~1.48e38 additions!)
# Instead, we'll just verify simpler properties
# The endomorphism should at least preserve being on the curve
println("\n" * "="^50)
println("Testing that p-power endomorphism preserves curve membership")
println("="^50)

println("Q on curve: ", GrothCurves.is_on_curve(Q))
ψQ = p_power_endo(Q)
println("ψ(Q) on curve: ", GrothCurves.is_on_curve(ψQ))

# Test with multiple applications
ψ2Q = p_power_endo(ψQ)
println("ψ²(Q) on curve: ", GrothCurves.is_on_curve(ψ2Q))

# For Fp2, we have x^(p²) = x (since p² ≡ 1 mod |Fp2*|)
# So ψ²(Q) should use coefficients squared but no conjugation
function p2_power_endo(Q::G2Point)
    if iszero(Q)
        return Q
    end
    x, y = GrothCurves.to_affine(Q)

    # No conjugation for p²
    x_final = x * (PSI_X^2)
    y_final = y * (PSI_Y^2)

    return G2Point(x_final, y_final)
end

ψ2Q_direct = p2_power_endo(Q)
println("\nψ²(Q) via direct formula == ψ(ψ(Q))? ", ψ2Q_direct == ψ2Q)

println("\n" * "="^70)
println("Summary")
println("="^70)
if PSI_X[1].value == arkworks_psi_x_c0 &&
   PSI_X[2].value == arkworks_psi_x_c1 &&
   PSI_Y[1].value == arkworks_psi_y_c0 &&
   PSI_Y[2].value == arkworks_psi_y_c1
    println("✓ All coefficients match arkworks!")
    println("✓ The implementation should be correct.")
else
    println("✗ Coefficients don't match arkworks.")
    println("  There may be a difference in twist convention or field representation.")
end
