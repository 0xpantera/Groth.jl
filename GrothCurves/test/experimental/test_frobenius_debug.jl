# Debug test for Frobenius endomorphism composition on G2
using GrothCurves
using GrothAlgebra

println("="^70)
println("Debugging Frobenius Endomorphism on G2")
println("="^70)

# Get the generator
Q = g2_generator()
println("\nTesting with G2 generator:")
x, y = GrothCurves.to_affine(Q)
println("  x = ", x)
println("  y = ", y)

# Test 1: Verify that Frobenius on Fp2 works as expected
println("\n" * "="^50)
println("Test 1: Frobenius on Fp2 elements")
println("="^50)

test_elem = Fp2Element(bn254_fq(123), bn254_fq(456))
println("\nOriginal: ", test_elem)
conj1 = GrothCurves.conjugate(test_elem)
println("After 1 Frobenius: ", conj1)
conj2 = GrothCurves.conjugate(conj1)
println("After 2 Frobenius: ", conj2)
println("Back to original? ", conj2 == test_elem)

# Test 2: Check the p-power endomorphism
println("\n" * "="^50)
println("Test 2: p-power endomorphism")
println("="^50)

π_Q = GrothCurves.p_power_endomorphism(Q)
println("\nπ(Q) computed")
x_pi, y_pi = GrothCurves.to_affine(π_Q)
println("  x = ", x_pi)
println("  y = ", y_pi)
println("  Still on curve? ", GrothCurves.is_on_curve(π_Q))

# Test 3: Check composition π(π(Q))
println("\n" * "="^50)
println("Test 3: Composition π(π(Q))")
println("="^50)

π2_Q_composed = GrothCurves.p_power_endomorphism(π_Q)
println("\nπ(π(Q)) computed via composition")
x_pi2_comp, y_pi2_comp = GrothCurves.to_affine(π2_Q_composed)
println("  x = ", x_pi2_comp)
println("  y = ", y_pi2_comp)
println("  Still on curve? ", GrothCurves.is_on_curve(π2_Q_composed))

# Test 4: Check direct π²(Q) formula
println("\n" * "="^50)
println("Test 4: Direct π²(Q) formula")
println("="^50)

π2_Q_direct = GrothCurves.frobenius_g2(Q, 2)
println("\nπ²(Q) computed via direct formula")
x_pi2_direct, y_pi2_direct = GrothCurves.to_affine(π2_Q_direct)
println("  x = ", x_pi2_direct)
println("  y = ", y_pi2_direct)
println("  Still on curve? ", GrothCurves.is_on_curve(π2_Q_direct))

# Test 5: Compare the two methods
println("\n" * "="^50)
println("Test 5: Comparing π(π(Q)) vs π²(Q) direct")
println("="^50)

println("\nπ(π(Q)) == π²(Q) direct? ", π2_Q_composed == π2_Q_direct)
if π2_Q_composed != π2_Q_direct
    println("MISMATCH DETECTED!")
    println("\nDifference in x:")
    println("  π(π(Q)).x = ", x_pi2_comp)
    println("  π²(Q).x   = ", x_pi2_direct)
    println("\nDifference in y:")
    println("  π(π(Q)).y = ", y_pi2_comp)
    println("  π²(Q).y   = ", y_pi2_direct)
end

# Test 6: Check coefficient squaring
println("\n" * "="^50)
println("Test 6: Checking coefficient squaring")
println("="^50)

# Get the coefficients from our constants
PSI_X = GrothCurves.P_POWER_ENDOMORPHISM_COEFF_0
PSI_Y = GrothCurves.P_POWER_ENDOMORPHISM_COEFF_1
PSI_X_SQ = GrothCurves.P_POWER_ENDOMORPHISM_COEFF_0_SQ
PSI_Y_SQ = GrothCurves.P_POWER_ENDOMORPHISM_COEFF_1_SQ

println("\nPSI_X^2 == PSI_X_SQ? ", PSI_X^2 == PSI_X_SQ)
println("PSI_Y^2 == PSI_Y_SQ? ", PSI_Y^2 == PSI_Y_SQ)

if PSI_X^2 != PSI_X_SQ
    println("ERROR: PSI_X squaring mismatch!")
    println("  PSI_X^2 = ", PSI_X^2)
    println("  PSI_X_SQ = ", PSI_X_SQ)
end

if PSI_Y^2 != PSI_Y_SQ
    println("ERROR: PSI_Y squaring mismatch!")
    println("  PSI_Y^2 = ", PSI_Y^2)
    println("  PSI_Y_SQ = ", PSI_Y_SQ)
end

# Test 7: Manual computation of π²(Q)
println("\n" * "="^50)
println("Test 7: Manual computation of π²(Q)")
println("="^50)

# For π²(Q), we should:
# 1. NOT apply conjugation (since Frob^2 = identity on Fp2)
# 2. Multiply by squared coefficients
x_manual_pi2 = x * (PSI_X^2)
y_manual_pi2 = y * (PSI_Y^2)
Q_manual_pi2 = G2Point(x_manual_pi2, y_manual_pi2)

println("\nManual π²(Q) computed")
println("  x = ", x_manual_pi2)
println("  y = ", y_manual_pi2)
println("  Still on curve? ", GrothCurves.is_on_curve(Q_manual_pi2))
println("  Matches π(π(Q))? ", Q_manual_pi2 == π2_Q_composed)
println("  Matches direct π²(Q)? ", Q_manual_pi2 == π2_Q_direct)

# Test 8: Check if coefficients satisfy expected relations
println("\n" * "="^50)
println("Test 8: Coefficient relations")
println("="^50)

# The coefficients are ξ^((p-1)/3) and ξ^((p-1)/2)
# Let's verify some properties
ξ = Fp2Element(bn254_fq(9), bn254_fq(1))
p = GrothCurves.BN254_PRIME

# PSI_X should be ξ^((p-1)/3)
expected_PSI_X = ξ^(div(p - 1, 3))
println("\nPSI_X matches ξ^((p-1)/3)? ", PSI_X == expected_PSI_X)

# PSI_Y should be ξ^((p-1)/2)
expected_PSI_Y = ξ^(div(p - 1, 2))
println("PSI_Y matches ξ^((p-1)/2)? ", PSI_Y == expected_PSI_Y)

# For PSI_Y, we have a special property: PSI_Y^2 should equal ξ^(p-1)
# And ξ^(p-1) should be in Fp (conjugate to itself)
psi_y_squared = PSI_Y^2
xi_p_minus_1 = ξ^(p - 1)
println("\nPSI_Y^2 = ", psi_y_squared)
println("ξ^(p-1) = ", xi_p_minus_1)
println("PSI_Y^2 == ξ^(p-1)? ", psi_y_squared == xi_p_minus_1)

# Check if ξ^(p-1) is in Fp (imaginary part should be 0)
println("ξ^(p-1) is in Fp? ", xi_p_minus_1[2] == bn254_fq(0))

# Summary
println("\n" * "="^70)
println("Summary")
println("="^70)

all_good = true
if π2_Q_composed != π2_Q_direct
    println("✗ π(π(Q)) != π²(Q) direct - COMPOSITION ERROR")
    all_good = false
end
if PSI_X^2 != PSI_X_SQ || PSI_Y^2 != PSI_Y_SQ
    println("✗ Coefficient squaring mismatch - CONSTANT ERROR")
    all_good = false
end
if Q_manual_pi2 != π2_Q_composed
    println("✗ Manual computation doesn't match composition - LOGIC ERROR")
    all_good = false
end

if all_good
    println("✓ All Frobenius tests passed!")
else
    println("\nErrors detected. The Frobenius implementation needs fixing.")
end
