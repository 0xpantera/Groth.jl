using GrothCurves

# BN254 group order r
const BN254_GROUP_ORDER = parse(BigInt, "21888242871839275222246405745257275088548364400416034343698204186575808495617")

# Test final exponentiation properties
P = g1_generator()
Q = g2_generator()

# Get a pairing result
f = GrothCurves.miller_loop(P, Q)
result = GrothCurves.final_exponentiation(f)

println("=== Testing Final Exponentiation Properties ===\n")

# Test 1: Subgroup check - result^r should equal 1
println("Test 1: Subgroup check (result^r == 1)")
result_to_r = result^BN254_GROUP_ORDER
is_in_subgroup = isone(result_to_r)
println("  result^r == 1? ", is_in_subgroup)
if !is_in_subgroup
    println("  ERROR: Result not in correct subgroup!")
end

# Test 2: Cyclotomic check - for GT elements, a^(p^6) == a^(-1)
println("\nTest 2: Cyclotomic/unitary check (a^(p^6) == a^(-1))")
result_p6 = GrothCurves.frobenius_map(result, 6)
result_inv = inv(result)
is_cyclotomic = result_p6 == result_inv
println("  frobenius_map(result, 6) == inv(result)? ", is_cyclotomic)
if !is_cyclotomic
    println("  ERROR: Result not in cyclotomic subgroup!")
end

# Test 3: Compare with naive final exponentiation
println("\nTest 3: Naive final exponentiation comparison")

# Compute the full exponent (p^12 - 1) / r
p = GrothCurves.BN254_PRIME
p12_minus_1 = p^12 - 1
full_exponent = div(p12_minus_1, BN254_GROUP_ORDER)

println("  Computing naive f^((p^12-1)/r) using binary exponentiation...")
println("  This will take a moment...")

# Naive exponentiation
function naive_final_exp(f::Fp12Element)
    return f^full_exponent
end

result_naive = naive_final_exp(f)

println("  Optimized result == Naive result? ", result == result_naive)

# Now test bilinearity with naive final exp
println("\n=== Testing Bilinearity with NAIVE Final Exponentiation ===\n")

function pairing_naive(P::G1Point, Q::G2Point)
    f = GrothCurves.miller_loop(P, Q)
    return naive_final_exp(f)
end

e_naive = pairing_naive(P, Q)
e_2P_naive = pairing_naive(2*P, Q)
e_P2Q_naive = pairing_naive(P, 2*Q)
e_squared_naive = e_naive^2

println("With naive final exp:")
println("  e(2P, Q) == e(P, Q)^2: ", e_2P_naive == e_squared_naive)
println("  e(P, 2Q) == e(P, Q)^2: ", e_P2Q_naive == e_squared_naive)
println("  e(2P, Q) == e(P, 2Q): ", e_2P_naive == e_P2Q_naive)

# Also test negation with naive
e_negP_naive = pairing_naive(-P, Q)
e_negQ_naive = pairing_naive(P, -Q)
println("\nNegation with naive final exp:")
println("  e(-P, Q) == inv(e(P, Q)): ", e_negP_naive == inv(e_naive))
println("  e(P, -Q) == inv(e(P, Q)): ", e_negQ_naive == inv(e_naive))