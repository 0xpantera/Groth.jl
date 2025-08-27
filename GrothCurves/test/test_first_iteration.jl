using GrothCurves

P = g1_generator()
Q = g2_generator()

println("=== First Miller Loop Iteration Trace ===\n")

# Initial values
T = Q
f = one(Fp12Element)

# Get P coordinates
Px, Py = GrothCurves.to_affine(P)
println("P coordinates:")
println("  x = ", Px.value)
println("  y = ", Py.value)

# Get Q coordinates  
xQ, yQ = GrothCurves.to_affine(Q)
println("\nQ coordinates:")
println("  x[0] = ", xQ[1].value)
println("  x[1] = ", xQ[2].value)
println("  y[0] = ", yQ[1].value)
println("  y[1] = ", yQ[2].value)

# First iteration is always a doubling (MSB is 1, we skip it)
println("\n=== First Doubling Step ===")

T_new, line_coeffs = GrothCurves.doubling_step(T)

println("\nLine coefficients:")
println("  a[0] = ", line_coeffs.a[1].value)
println("  a[1] = ", line_coeffs.a[2].value)
println("  b[0] = ", line_coeffs.b[1].value)
println("  b[1] = ", line_coeffs.b[2].value)
println("  c[0] = ", line_coeffs.c[1].value)
println("  c[1] = ", line_coeffs.c[2].value)

# Check line vanishing
xT, yT = GrothCurves.to_affine(T)
val_at_T = line_coeffs.a * xT + line_coeffs.b * yT + line_coeffs.c
println("\nLine vanishes at T? ", iszero(val_at_T))

# Evaluate line at P
line_eval = GrothCurves.evaluate_line(line_coeffs, P)
println("\nLine evaluation (sparse Fp12):")
println("  c0.c0 (slot 0, b*yP) = ", line_eval[1][1])
println("  c1.c0 (slot 1, a*xP) = ", line_eval[2][1])
println("  c1.c1 (slot 4, c) = ", line_eval[2][2])

# Update f
f_new = f^2 * line_eval
println("\nAfter first doubling:")
println("  f^2 * line_eval first component = ", f_new[1][1][1].value)

# Get new T coordinates
xT_new, yT_new = GrothCurves.to_affine(T_new)
println("\n2T coordinates:")
println("  x[0] = ", xT_new[1].value)
println("  x[1] = ", xT_new[2].value)
println("  y[0] = ", yT_new[1].value)
println("  y[1] = ", yT_new[2].value)

# Check if 2T is correct
T_check = 2 * Q
xT_check, yT_check = GrothCurves.to_affine(T_check)
println("\n2T matches 2*Q? ", xT_new == xT_check && yT_new == yT_check)

# Let's also check with a known test vector if possible
# For BN254 with standard generators, after first doubling:
# We should have specific values - but we need a reference to compare against

println("\n=== Compare current 0-1-4 mapping vs ground truth ===")

# Our current mapping
println("\nCurrent 0-1-4 mapping:")
println("  slot 0 (c0.c0) ← b*yP")
println("  slot 1 (c1.c0) ← a*xP")  
println("  slot 4 (c1.c1) ← c")

# What it should be for M-twist
println("\nExpected for M-twist with ψ: x→xv, y→yw:")
println("  Line ℓ(ψ(P)) = (a*(xP*v) + c) + (b*yP)*w")
println("  = Fp12((a*xP*v + c), b*yP)")
println("  where v = (0,1,0) in Fp6")

# So the correct sparse form should be:
# c0 contains: a*xP*v + c
# c1 contains: b*yP

println("\nThis means:")
println("  c0 should have a*xP in slot v (position 1 of Fp6) and c in slot 0")
println("  c1 should have b*yP in slot 0")
println("  Our current mapping may be wrong!")