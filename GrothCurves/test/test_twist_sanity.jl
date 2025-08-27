using GrothCurves

const r = parse(BigInt,
 "21888242871839275222246405745257275088548364400416034343698204186575808495617")

P = g1_generator()
Q = g2_generator()

println("=== Critical Sanity Checks ===\n")

println("G1 on-curve: ", GrothCurves.is_on_curve(P))
println("G2 on-curve (current b_twist): ", GrothCurves.is_on_curve(Q))

println("r*P == O: ", iszero(r*P))
println("r*Q == O: ", iszero(r*Q))

println("\n=== Testing M-twist vs D-twist ===\n")

# Try D-twist temporarily to test the hypothesis:
# y^2 = x^3 + 3*ξ  instead of 3/ξ
ξ = Fp2Element(9, 1)
b_twist_D = Fp2Element(3) * ξ
b_twist_M = Fp2Element(3) / ξ

# Quick probe with the current Q:
X, Y, Z = GrothCurves.x_coord(Q), GrothCurves.y_coord(Q), GrothCurves.z_coord(Q)
Z2, Z4, Z6 = Z^2, Z^4, Z^6

lhs = Y^2
rhs_M = X^3 + b_twist_M * Z6
rhs_D = X^3 + b_twist_D * Z6

println("Q satisfies M-twist (y² = x³ + 3/ξ)? ", lhs == rhs_M)
println("Q satisfies D-twist (y² = x³ + 3*ξ)? ", lhs == rhs_D)

# Let's also check what our code thinks b_twist is
println("\nOur G2_B_TWIST constant: ", GrothCurves.G2_B_TWIST)
println("Expected M-twist (3/ξ): ", b_twist_M)
println("Expected D-twist (3*ξ): ", b_twist_D)
println("Do they match? M-twist: ", GrothCurves.G2_B_TWIST == b_twist_M)

# Check if the generator is actually on ANY curve
println("\n=== Debugging G2 generator ===")
xQ, yQ = GrothCurves.to_affine(Q)
println("G2 generator affine coords:")
println("  x = ", xQ[1])  # Real part
println("  y = ", yQ[1])  # Real part

# Check the curve equation directly in affine
lhs_affine = yQ^2
rhs_M_affine = xQ^3 + b_twist_M
rhs_D_affine = xQ^3 + b_twist_D

println("\nIn affine coordinates:")
println("  y² = ", lhs_affine[1])
println("  x³ + 3/ξ = ", rhs_M_affine[1])
println("  x³ + 3*ξ = ", rhs_D_affine[1])
println("  M-twist match? ", lhs_affine == rhs_M_affine)
println("  D-twist match? ", lhs_affine == rhs_D_affine)