# Test D-twist line evaluation for BN254 pairing
using GrothCurves
using GrothAlgebra

println("="^70)
println("Testing D-twist Line Evaluation")
println("="^70)

# Get test points
P = g1_generator()
Q = g2_generator()

# Get affine coordinates
P_x, P_y = GrothCurves.to_affine(P)
Q_x, Q_y = GrothCurves.to_affine(Q)

println("\nTest points:")
println("P (G1):")
println("  x = ", P_x)
println("  y = ", P_y)
println("\nQ (G2):")
println("  x = ", Q_x)
println("  y = ", Q_y)

# Test 1: Line evaluation structure
println("\n" * "="^50)
println("Test 1: Line Evaluation Structure")
println("="^50)

# Perform a doubling step to get line coefficients
T, line_coeffs = GrothCurves.doubling_step(Q)

println("\nLine coefficients from doubling:")
println("  a (y-coeff) = ", line_coeffs.a)
println("  b (x-coeff) = ", line_coeffs.b)
println("  c (constant) = ", line_coeffs.c)

# Evaluate the line at P
line_eval = GrothCurves.evaluate_line(line_coeffs, P)

println("\nLine evaluation at P:")
println("  Result is Fp12 element")

# Check the structure - for D-twist, we should have sparse multiplication format (0,3,4)
# This means c0 has non-zero at position 0, c1 has non-zero at positions 0,1
c0 = line_eval[1]  # First Fp6 component
c1 = line_eval[2]  # Second Fp6 component

println("\nFp12 structure (c0, c1):")
println("  c0[0] = ", c0[1])  # Should be a*yP
println("  c0[1] = ", c0[2])  # Should be 0
println("  c0[2] = ", c0[3])  # Should be 0
println("  c1[0] = ", c1[1])  # Should be b*xP
println("  c1[1] = ", c1[2])  # Should be c
println("  c1[2] = ", c1[3])  # Should be 0

# Verify the sparse structure
is_034_sparse = (c0[2] == zero(Fp2Element)) &&
                (c0[3] == zero(Fp2Element)) &&
                (c1[3] == zero(Fp2Element))

println("\nIs (0,3,4) sparse? ", is_034_sparse)

# Test 2: Verify D-twist formula
println("\n" * "="^50)
println("Test 2: D-twist Formula Verification")
println("="^50)

# Convert P coordinates to Fp2 for computation
xP_fp2 = Fp2Element(P_x, zero(BN254Fq))
yP_fp2 = Fp2Element(P_y, zero(BN254Fq))

# According to D-twist, the line evaluation should be:
# a*yP at position (0,0,0) in Fp12
# b*xP at position (1,0,0) in Fp12
# c at position (1,1,0) in Fp12

expected_c0_0 = line_coeffs.a * yP_fp2
expected_c1_0 = line_coeffs.b * xP_fp2
expected_c1_1 = line_coeffs.c

println("\nExpected values:")
println("  c0[0] should be a*yP = ", expected_c0_0)
println("  c1[0] should be b*xP = ", expected_c1_0)
println("  c1[1] should be c = ", expected_c1_1)

println("\nActual values:")
println("  c0[0] = ", c0[1])
println("  c1[0] = ", c1[1])
println("  c1[1] = ", c1[2])

println("\nMatches:")
println("  c0[0] matches? ", c0[1] == expected_c0_0)
println("  c1[0] matches? ", c1[1] == expected_c1_0)
println("  c1[1] matches? ", c1[2] == expected_c1_1)

# Test 3: Compare with addition step
println("\n" * "="^50)
println("Test 3: Addition Step Line Evaluation")
println("="^50)

# Perform an addition step
T2, line_coeffs2 = GrothCurves.addition_step(T, Q)

println("\nLine coefficients from addition:")
println("  a (y-coeff) = ", line_coeffs2.a)
println("  b (x-coeff) = ", line_coeffs2.b)
println("  c (constant) = ", line_coeffs2.c)

# Evaluate the line at P
line_eval2 = GrothCurves.evaluate_line(line_coeffs2, P)

# Check the structure
c0_2 = line_eval2[1]
c1_2 = line_eval2[2]

is_034_sparse_2 = (c0_2[2] == zero(Fp2Element)) &&
                  (c0_2[3] == zero(Fp2Element)) &&
                  (c1_2[3] == zero(Fp2Element))

println("\nIs addition line (0,3,4) sparse? ", is_034_sparse_2)

# Test 4: Special cases
println("\n" * "="^50)
println("Test 4: Special Cases")
println("="^50)

# Test with point at infinity
O_G1 = zero(G1Point)
line_eval_O = GrothCurves.evaluate_line(line_coeffs, O_G1)

println("\nLine evaluation at O (point at infinity):")
c0_O = line_eval_O[1]
c1_O = line_eval_O[2]

# For point at infinity, P_x = P_y = 0, so we should get:
# c0[0] = a*0 = 0
# c1[0] = b*0 = 0
# c1[1] = c (unchanged)

println("  c0[0] = ", c0_O[1], " (should be 0)")
println("  c1[0] = ", c1_O[1], " (should be 0)")
println("  c1[1] = ", c1_O[2], " (should be c = ", line_coeffs.c, ")")

is_correct_at_O = (c0_O[1] == zero(Fp2Element)) &&
                  (c1_O[1] == zero(Fp2Element)) &&
                  (c1_O[2] == line_coeffs.c)

println("\nCorrect at O? ", is_correct_at_O)

# Test 5: Multiplication properties
println("\n" * "="^50)
println("Test 5: Sparse Multiplication Test")
println("="^50)

# Create a random Fp12 element
one_fp12 = one(Fp12Element)
result = one_fp12 * line_eval

println("\nMultiplying 1 by line evaluation:")
println("  Result == line_eval? ", result == line_eval)

# Square the line evaluation
line_eval_squared = line_eval^2

println("\nSquaring line evaluation:")
println("  Result computed successfully")

# The sparse structure should be preserved in some sense
# though squaring will generally make it dense

# Summary
println("\n" * "="^70)
println("Summary")
println("="^70)

all_correct = is_034_sparse &&
              (c0[1] == expected_c0_0) &&
              (c1[1] == expected_c1_0) &&
              (c1[2] == expected_c1_1) &&
              is_034_sparse_2 &&
              is_correct_at_O

if all_correct
    println("✓ All D-twist line evaluation tests passed!")
    println("The D-twist formulas appear to be correctly implemented.")
else
    println("✗ Some D-twist tests failed!")
    println("Issues detected:")
    if !is_034_sparse
        println("  - Line evaluation does not have (0,3,4) sparse structure")
    end
    if c0[1] != expected_c0_0
        println("  - c0[0] component mismatch")
    end
    if c1[1] != expected_c1_0
        println("  - c1[0] component mismatch")
    end
    if c1[2] != expected_c1_1
        println("  - c1[1] component mismatch")
    end
    if !is_034_sparse_2
        println("  - Addition line not (0,3,4) sparse")
    end
    if !is_correct_at_O
        println("  - Incorrect evaluation at point at infinity")
    end
end
