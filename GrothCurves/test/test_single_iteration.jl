using GrothCurves

# Trace a single iteration to debug
P = g1_generator()
Q = g2_generator()

println("=== Single Miller Loop Iteration Debug ===")

# Initialize
T = Q
f = one(Fp12Element)

# Get loop bits
loop_bits = digits(Bool, GrothCurves.ATE_LOOP_COUNT, base=2)

println("\nLoop uses ", length(loop_bits)-1, " iterations")
println("First few bits (MSB first): ", reverse(loop_bits[end-5:end]))

# Do just the first iteration (MSB)
println("\n--- First iteration (doubling only, MSB is always 1 so we skip it) ---")

# Actually, let's find the first iteration with both doubling and addition
global first_add_bit = 0
for i in (length(loop_bits)-1):-1:1
    if loop_bits[i]
        global first_add_bit = i
        break
    end
end

println("First iteration with addition at bit position: ", first_add_bit)

# Process up to and including that iteration
global iter_count = 0
for i in (length(loop_bits)-1):-1:first_add_bit
    global iter_count += 1
    global f, T
    
    # Doubling step
    f_before_double = f
    T_new, line_coeffs = GrothCurves.doubling_step(T)
    line_eval = GrothCurves.evaluate_line(line_coeffs, P)
    f = f^2 * line_eval
    
    if i == first_add_bit
        println("\n=== Iteration ", iter_count, " (bit position ", i, ") ===")
        println("DOUBLING:")
        println("  Line coeffs: a=", line_coeffs.a[1], " (first Fp2 component)")
        println("              b=", line_coeffs.b[1])
        println("              c=", line_coeffs.c[1])
        println("  Line eval sparse Fp12:")
        println("    c0.c0 (slot 0) = ", line_eval[1][1])  # b*yP
        println("    c1.c0 (slot 1) = ", line_eval[2][1])  # a*xP
        println("    c1.c1 (slot 4) = ", line_eval[2][2])  # c
        println("  f changed: ", f != f_before_double^2 * line_eval)
    end
    
    T = T_new
    
    # Addition step if bit is 1
    if loop_bits[i]
        f_before_add = f
        T_new, line_coeffs = GrothCurves.addition_step(T, Q)
        line_eval = GrothCurves.evaluate_line(line_coeffs, P)
        f = f * line_eval
        
        if i == first_add_bit
            println("\nADDITION:")
            println("  Line coeffs: a=", line_coeffs.a[1], " (first Fp2 component)")
            println("              b=", line_coeffs.b[1])
            println("              c=", line_coeffs.c[1])
            println("  Line eval sparse Fp12:")
            println("    c0.c0 (slot 0) = ", line_eval[1][1])  # b*yP
            println("    c1.c0 (slot 1) = ", line_eval[2][1])  # a*xP
            println("    c1.c1 (slot 4) = ", line_eval[2][2])  # c
            println("  f changed: ", f != f_before_add * line_eval)
        end
        
        T = T_new
    end
end

println("\n=== Checking line vanishing ===")
# Check that the lines actually vanish where they should
xT, yT = GrothCurves.to_affine(T)
xQ, yQ = GrothCurves.to_affine(Q)

# Get a fresh doubling line
T_test = 2*Q
T_double, line_double = GrothCurves.doubling_step(T_test)
xT_test, yT_test = GrothCurves.to_affine(T_test)
val_at_T = line_double.a * xT_test + line_double.b * yT_test + line_double.c
println("Doubling line vanishes at T: ", iszero(val_at_T))

# Get a fresh addition line
T_add, line_add = GrothCurves.addition_step(T_test, Q)
val_at_Q_add = line_add.a * xQ + line_add.b * yQ + line_add.c
val_at_T_add = line_add.a * xT_test + line_add.b * yT_test + line_add.c
println("Addition line vanishes at Q: ", iszero(val_at_Q_add))
println("Addition line vanishes at T: ", iszero(val_at_T_add))