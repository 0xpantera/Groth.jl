using GrothCurves

println("=== Testing G2 Frobenius Endomorphism ===\n")

Q = g2_generator()
xQ, yQ = GrothCurves.to_affine(Q)

# Test Frobenius properties
Q_pi = GrothCurves.frobenius_g2(Q, 1)
Q_pi2 = GrothCurves.frobenius_g2(Q, 2)
Q_pi3 = GrothCurves.frobenius_g2(Q, 3)

xQ_pi, yQ_pi = GrothCurves.to_affine(Q_pi)
xQ_pi2, yQ_pi2 = GrothCurves.to_affine(Q_pi2)

println("Original Q coordinates:")
println("  x = ", xQ[1])
println("  y = ", yQ[1])

println("\nπ(Q) coordinates:")
println("  x = ", xQ_pi[1])
println("  y = ", yQ_pi[1])

println("\nπ²(Q) coordinates:")
println("  x = ", xQ_pi2[1])
println("  y = ", yQ_pi2[1])

# Check gamma constants
println("\n=== Gamma Constants ===")
println("γ₁[2] = ", GrothCurves.G2_GAMMA1[2])
println("γ₁[3] = ", GrothCurves.G2_GAMMA1[3])
println("γ₂[2] = ", GrothCurves.G2_GAMMA2[2])
println("γ₂[3] = ", GrothCurves.G2_GAMMA3[3])

# Verify Frobenius is an endomorphism
println("\n=== Frobenius Endomorphism Properties ===")

# π^p should be identity on G2
p = GrothCurves.BN254_PRIME
println("π^p(Q) == Q? (should be true for p-power Frobenius)")

# Check if points are on curve
println("π(Q) on curve? ", GrothCurves.is_on_curve(Q_pi))
println("π²(Q) on curve? ", GrothCurves.is_on_curve(Q_pi2))

# Check order
r = parse(BigInt, "21888242871839275222246405745257275088548364400416034343698204186575808495617")
println("r*π(Q) == O? ", iszero(r*Q_pi))
println("r*π²(Q) == O? ", iszero(r*Q_pi2))

# Check the correction formula: π²(Q) - p*π(Q) + Q should equal O
# This is the key relation for BN254 optimal ate pairing
println("\n=== Key BN254 Relation ===")
# Note: p*π(Q) is a huge scalar mult, let's check a different way

# The relation for optimal ate is actually about the trace map
# For BN254: π²(Q) + Q = t*π(Q) where t is the trace
# Let's verify the endomorphism structure

# First, let's check if our Frobenius implementation matches the expected formula
println("\nVerifying Frobenius formulas:")

# For M-twist, π¹: (x,y) → (x̄·γ₁[2], ȳ·γ₁[3])
x_pi_expected = GrothCurves.conjugate(xQ) * GrothCurves.G2_GAMMA1[2]
y_pi_expected = GrothCurves.conjugate(yQ) * GrothCurves.G2_GAMMA1[3]

println("π(x) matches expected? ", xQ_pi == x_pi_expected)
println("π(y) matches expected? ", yQ_pi == y_pi_expected)

# For π²: (x,y) → (x·γ₂[2], y·γ₂[3]) - no conjugation for even powers
x_pi2_expected = xQ * GrothCurves.G2_GAMMA2[2]
y_pi2_expected = yQ * GrothCurves.G2_GAMMA2[3]

println("π²(x) matches expected? ", xQ_pi2 == x_pi2_expected)
println("π²(y) matches expected? ", yQ_pi2 == y_pi2_expected)