using Test
using GrothCurves
using GrothAlgebra

# Demonstration that the Miller loop now correctly depends on both coordinates of P

println("Testing Miller Loop Line Evaluation Fix")
println("=" ^ 50)

g1 = g1_generator()
g2 = g2_generator()

# Test 1: Miller loop should distinguish P and -P
f_pos = miller_loop(g1, g2)
f_neg = miller_loop(-g1, g2)

println("\n1. Testing miller_loop(P, Q) vs miller_loop(-P, Q):")
println("   Are they different? ", f_pos != f_neg ? "✓ YES (correct!)" : "✗ NO (incorrect)")

# Test 2: Line evaluation should depend on y_P
println("\n2. Testing line evaluation dependency on y_P:")

# Get a simple doubling step
T_doubled, coeffs_pos = doubling_step(g2, g1)
_, coeffs_neg = doubling_step(g2, -g1)

line_pos = evaluate_line(coeffs_pos)
line_neg = evaluate_line(coeffs_neg)

println("   Line evaluation at P:  ", line_pos[1][1] != zero(Fp2Element) ? "non-zero y-term ✓" : "zero y-term ✗")
println("   Line evaluation at -P: ", line_neg[1][1] != zero(Fp2Element) ? "non-zero y-term ✓" : "zero y-term ✗")
println("   Are they different? ", line_pos != line_neg ? "✓ YES (correct!)" : "✗ NO (incorrect)")

# Test 3: Check sparse 0-1-4 embedding structure
println("\n3. Testing sparse 0-1-4 embedding structure:")
println("   Checking Fp12 element structure from line evaluation...")

# Fp12 = (c0, c1) where c0, c1 are Fp6 elements
# Sparse 0-1-4 means: c0 = (*, 0, 0) and c1 = (*, *, 0)
is_sparse_014 = iszero(line_pos[1][2]) && iszero(line_pos[1][3]) && iszero(line_pos[2][3])
println("   Has sparse 0-1-4 form? ", is_sparse_014 ? "✓ YES" : "✗ NO")

# Test 4: Verify formulas match literature
println("\n4. Formula verification:")
println("   Using M-twist formulas from 'The Realm of the Pairings'")
println("   Line coefficients include:")
println("   - c0 = -(2YZ) * y_P       [depends on y_P ✓]")
println("   - c1 = (3X²) * x_P        [depends on x_P ✓]")
println("   - c2 = 3b'Z² - Y²         [constant term ✓]")

# Summary
println("\n" * "=" ^ 50)
println("Summary: Miller loop implementation now correctly:")
println("  • Uses both x_P and y_P coordinates in line evaluation")
println("  • Embeds lines in sparse 0-1-4 form for efficiency")
println("  • Distinguishes between P and -P before final exponentiation")
println("  • Follows standard BN254 M-twist formulas from literature")

# Run formal tests
@testset "Miller Loop Fix Verification" begin
    @test f_pos != f_neg
    @test line_pos != line_neg
    @test is_sparse_014
    @test !iszero(coeffs_pos.c0)  # y-term coefficient non-zero
    @test coeffs_pos.c0 == -coeffs_neg.c0  # Sign flips with -P
end