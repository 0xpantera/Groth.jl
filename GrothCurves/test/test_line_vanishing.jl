using GrothCurves

println("Testing that lines vanish at the correct points...")

Q = g2_generator()
T = 2 * Q  # Some other point

# Test addition line
T_plus_Q, line_add = GrothCurves.addition_step(T, Q)

# Get affine coordinates
xQ, yQ = GrothCurves.to_affine(Q)
xT, yT = GrothCurves.to_affine(T)

# Check that the line vanishes at both Q and T
at_Q = line_add.a * xQ + line_add.b * yQ + line_add.c
at_T = line_add.a * xT + line_add.b * yT + line_add.c

println("\nAddition line (T + Q):")
println("At Q: ", at_Q, " (should be zero)")
println("At T: ", at_T, " (should be zero)")
println("Both zero? ", iszero(at_Q) && iszero(at_T))

# Test doubling line
T2, line_double = GrothCurves.doubling_step(T)

# Check that the tangent line vanishes at T
at_T_double = line_double.a * xT + line_double.b * yT + line_double.c

println("\nDoubling line (tangent at T):")
println("At T: ", at_T_double, " (should be zero)")
println("Is zero? ", iszero(at_T_double))

# Also test with the generator
Q2, line_double_Q = GrothCurves.doubling_step(Q)
at_Q_double = line_double_Q.a * xQ + line_double_Q.b * yQ + line_double_Q.c

println("\nDoubling line (tangent at Q):")
println("At Q: ", at_Q_double, " (should be zero)")
println("Is zero? ", iszero(at_Q_double))

if iszero(at_Q) && iszero(at_T) && iszero(at_T_double) && iszero(at_Q_double)
    println("\n✅ All lines vanish at the correct points!")
else
    println("\n❌ Some lines don't vanish where they should")
end