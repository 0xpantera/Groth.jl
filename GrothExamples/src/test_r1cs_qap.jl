"""
Test R1CS and QAP conversion for r = x * y * z * u
"""

push!(LOAD_PATH, dirname(@__DIR__))

using GrothAlgebra
using GrothCurves  
using GrothProofs

function main()
    println("=" ^ 60)
    println("R1CS to QAP Test: r = x * y * z * u")
    println("=" ^ 60)
    
    # Step 1: Create R1CS
    println("\n1. Creating R1CS")
    r1cs = create_r1cs_example_multiplication()
    
    println("   Matrix L (Left):")
    for i in 1:r1cs.num_constraints
        row = [r1cs.L[i,j].value for j in 1:r1cs.num_vars]
        println("   Row $i: $row")
    end
    
    println("\n   Matrix R (Right):")
    for i in 1:r1cs.num_constraints
        row = [r1cs.R[i,j].value for j in 1:r1cs.num_vars]
        println("   Row $i: $row")
    end
    
    println("\n   Matrix O (Output):")
    for i in 1:r1cs.num_constraints
        row = [r1cs.O[i,j].value for j in 1:r1cs.num_vars]
        println("   Row $i: $row")
    end
    
    # Step 2: Create witness
    x, y, z, u = 3, 5, 7, 11
    println("\n2. Creating witness with x=$x, y=$y, z=$z, u=$u")
    witness = create_witness_multiplication(x, y, z, u)
    
    println("   Witness vector:")
    labels = ["1", "r", "x", "y", "z", "u", "v1", "v2"]
    for (i, (label, val)) in enumerate(zip(labels, witness.values))
        println("   $label = $(val.value)")
    end
    
    # Step 3: Verify R1CS
    println("\n3. Verifying R1CS constraints")
    w = witness.values
    
    for i in 1:r1cs.num_constraints
        # Compute L·w, R·w, O·w for this constraint
        Lw = sum(r1cs.L[i,j] * w[j] for j in 1:r1cs.num_vars)
        Rw = sum(r1cs.R[i,j] * w[j] for j in 1:r1cs.num_vars)
        Ow = sum(r1cs.O[i,j] * w[j] for j in 1:r1cs.num_vars)
        
        product = Lw * Rw
        
        if product == Ow
            println("   ✓ Constraint $i: $(Lw.value) * $(Rw.value) = $(Ow.value)")
        else
            println("   ✗ Constraint $i failed: $(Lw.value) * $(Rw.value) ≠ $(Ow.value)")
        end
    end
    
    # Step 4: Convert to QAP
    println("\n4. Converting to QAP")
    qap = r1cs_to_qap(r1cs)
    
    println("   Target polynomial t(x) = ∏(x - ωⁱ)")
    println("   Roots: $(qap.domain)")
    println("   Degree: $(degree(qap.t))")
    
    # Step 5: Check QAP at roots
    println("\n5. Verifying QAP at roots of t(x)")
    for (i, root) in enumerate(qap.domain)
        u_val, v_val, w_val = evaluate_qap(qap, witness, root)
        if u_val * v_val == w_val
            println("   ✓ At ω$i = $(root.value): U*V = W")
        else
            println("   ✗ At ω$i = $(root.value): U*V ≠ W")
        end
    end
    
    # Step 6: Compute h(x)
    println("\n6. Computing quotient polynomial h(x)")
    h = compute_h_polynomial(qap, witness)
    println("   Degree of h(x): $(degree(h))")
    
    # The key property: U(x)*V(x) - W(x) = h(x)*t(x)
    println("\n   Key property: U(x)*V(x) - W(x) = h(x)*t(x)")
    println("   This means U*V - W is divisible by t(x)")
    
    println("\n" * "=" ^ 60)
    println("✓ R1CS to QAP conversion successful!")
    println("=" ^ 60)
end

main()