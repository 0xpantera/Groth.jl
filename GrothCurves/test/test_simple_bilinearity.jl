# Simple bilinearity tests to isolate the issue
using GrothCurves
using GrothAlgebra

println("="^70)
println("Simple Bilinearity Tests")
println("="^70)

# Get generators
P = g1_generator()
Q = g2_generator()

# Test 1: Identity elements
println("\n" * "="^50)
println("Test 1: Identity Elements")
println("="^50)

O_G1 = zero(G1Point)
O_G2 = zero(G2Point)

e_O_Q = pairing(O_G1, Q)
e_P_O = pairing(P, O_G2)
e_O_O = pairing(O_G1, O_G2)

println("e(O, Q) = ", e_O_Q)
println("e(P, O) = ", e_P_O)
println("e(O, O) = ", e_O_O)
println()
println("e(O, Q) == 1? ", e_O_Q == one(Fp12Element))
println("e(P, O) == 1? ", e_P_O == one(Fp12Element))
println("e(O, O) == 1? ", e_O_O == one(Fp12Element))

# Test 2: Basic pairing value
println("\n" * "="^50)
println("Test 2: Basic Pairing")
println("="^50)

e_P_Q = pairing(P, Q)
println("e(P, Q) first component: ", e_P_Q[1][1][1])
println("e(P, Q) != 1? ", e_P_Q != one(Fp12Element))

# Test 3: Negation
println("\n" * "="^50)
println("Test 3: Negation Properties")
println("="^50)

neg_P = -P
neg_Q = -Q

e_negP_Q = pairing(neg_P, Q)
e_P_negQ = pairing(P, neg_Q)
e_negP_negQ = pairing(neg_P, neg_Q)

println("e(-P, Q) first component: ", e_negP_Q[1][1][1])
println("e(P, -Q) first component: ", e_P_negQ[1][1][1])
println("e(-P, -Q) first component: ", e_negP_negQ[1][1][1])
println()
println("e(-P, Q) == inv(e(P, Q))? ", e_negP_Q == inv(e_P_Q))
println("e(P, -Q) == inv(e(P, Q))? ", e_P_negQ == inv(e_P_Q))
println("e(-P, -Q) == e(P, Q)? ", e_negP_negQ == e_P_Q)

# Test 4: Simple doubling
println("\n" * "="^50)
println("Test 4: Doubling")
println("="^50)

P2 = P + P
Q2 = Q + Q

e_2P_Q = pairing(P2, Q)
e_P_2Q = pairing(P, Q2)
e_P_Q_squared = e_P_Q^2

println("e([2]P, Q) first component: ", e_2P_Q[1][1][1])
println("e(P, [2]Q) first component: ", e_P_2Q[1][1][1])
println("e(P, Q)^2 first component: ", e_P_Q_squared[1][1][1])
println()
println("e([2]P, Q) == e(P, [2]Q)? ", e_2P_Q == e_P_2Q)
println("e([2]P, Q) == e(P, Q)^2? ", e_2P_Q == e_P_Q_squared)
println("e(P, [2]Q) == e(P, Q)^2? ", e_P_2Q == e_P_Q_squared)

if e_2P_Q != e_P_2Q
    println("\nBILINEARITY FAILURE: e([2]P, Q) != e(P, [2]Q)")

    # Try to understand the relationship
    ratio1 = e_2P_Q / e_P_2Q
    ratio2 = e_P_2Q / e_2P_Q

    println("\nInvestigating the relationship:")
    println("e([2]P, Q) / e(P, [2]Q) first component: ", ratio1[1][1][1])
    println("e(P, [2]Q) / e([2]P, Q) first component: ", ratio2[1][1][1])
end

# Test 5: Tripling
println("\n" * "="^50)
println("Test 5: Tripling")
println("="^50)

P3 = P + P + P
Q3 = Q + Q + Q

e_3P_Q = pairing(P3, Q)
e_P_3Q = pairing(P, Q3)
e_P_Q_cubed = e_P_Q^3

println("e([3]P, Q) first component: ", e_3P_Q[1][1][1])
println("e(P, [3]Q) first component: ", e_P_3Q[1][1][1])
println("e(P, Q)^3 first component: ", e_P_Q_cubed[1][1][1])
println()
println("e([3]P, Q) == e(P, [3]Q)? ", e_3P_Q == e_P_3Q)
println("e([3]P, Q) == e(P, Q)^3? ", e_3P_Q == e_P_Q_cubed)
println("e(P, [3]Q) == e(P, Q)^3? ", e_P_3Q == e_P_Q_cubed)

# Test 6: Addition formula
println("\n" * "="^50)
println("Test 6: Addition Formula")
println("="^50)

# Test e(P1 + P2, Q) = e(P1, Q) * e(P2, Q)
e_P_plus_P_Q = pairing(P2, Q)  # P + P = [2]P
e_P_Q_times_e_P_Q = e_P_Q * e_P_Q

println("e(P + P, Q) first component: ", e_P_plus_P_Q[1][1][1])
println("e(P, Q) * e(P, Q) first component: ", e_P_Q_times_e_P_Q[1][1][1])
println()
println("e(P + P, Q) == e(P, Q) * e(P, Q)? ", e_P_plus_P_Q == e_P_Q_times_e_P_Q)

# Test 7: Cross products
println("\n" * "="^50)
println("Test 7: Cross Products")
println("="^50)

# Test e([2]P, [3]Q) == e(P, Q)^6
P2 = P + P
Q3 = Q + Q + Q
e_2P_3Q = pairing(P2, Q3)
e_P_Q_to_6 = e_P_Q^6

println("e([2]P, [3]Q) first component: ", e_2P_3Q[1][1][1])
println("e(P, Q)^6 first component: ", e_P_Q_to_6[1][1][1])
println()
println("e([2]P, [3]Q) == e(P, Q)^6? ", e_2P_3Q == e_P_Q_to_6)

# Summary
println("\n" * "="^70)
println("Summary of Bilinearity Tests")
println("="^70)

identity_ok = e_O_Q == one(Fp12Element) && e_P_O == one(Fp12Element)
non_degenerate = e_P_Q != one(Fp12Element)
negation_ok = e_negP_Q == inv(e_P_Q) && e_P_negQ == inv(e_P_Q) && e_negP_negQ == e_P_Q
doubling_bilinear = e_2P_Q == e_P_2Q
doubling_power = e_2P_Q == e_P_Q_squared && e_P_2Q == e_P_Q_squared
tripling_bilinear = e_3P_Q == e_P_3Q
tripling_power = e_3P_Q == e_P_Q_cubed && e_P_3Q == e_P_Q_cubed
addition_formula_ok = e_P_plus_P_Q == e_P_Q_times_e_P_Q
cross_product_ok = e_2P_3Q == e_P_Q_to_6

println("✓ Identity elements work: ", identity_ok)
println("✓ Non-degenerate: ", non_degenerate)
println("✓ Negation properties: ", negation_ok)
println("✓ Doubling bilinearity: ", doubling_bilinear)
println("✓ Doubling power rule: ", doubling_power)
println("✓ Tripling bilinearity: ", tripling_bilinear)
println("✓ Tripling power rule: ", tripling_power)
println("✓ Addition formula: ", addition_formula_ok)
println("✓ Cross product: ", cross_product_ok)

if doubling_bilinear && doubling_power && tripling_bilinear && tripling_power
    println("\n✓✓✓ ALL BILINEARITY TESTS PASSED! ✓✓✓")
else
    println("\n✗✗✗ BILINEARITY TESTS FAILED! ✗✗✗")
    println("The pairing implementation has bugs that need to be fixed.")
end
