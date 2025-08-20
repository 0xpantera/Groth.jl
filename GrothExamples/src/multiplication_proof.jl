"""
End-to-end example: Groth16 proof for r = x * y * z * u

This example demonstrates the complete workflow:
1. Create R1CS for the computation
2. Convert R1CS to QAP
3. Generate trusted setup
4. Create a proof
5. Verify the proof
"""

# Add the parent directory to the load path
push!(LOAD_PATH, dirname(@__DIR__))

using GrothAlgebra
using GrothCurves
using GrothProofs

function main()
    println("=" ^ 60)
    println("Groth16 Proof System Example: r = x * y * z * u")
    println("=" ^ 60)
    
    # Step 1: Create the R1CS
    println("\n1. Creating R1CS for r = x * y * z * u")
    println("   Variables: [1, r, x, y, z, u, v1, v2]")
    println("   Constraints:")
    println("   - v1 = x * y")
    println("   - v2 = z * u")
    println("   - r = v1 * v2")
    
    r1cs = create_r1cs_example_multiplication()
    println("   ✓ R1CS created with $(r1cs.num_constraints) constraints and $(r1cs.num_vars) variables")
    
    # Step 2: Create a witness with specific values
    x, y, z, u = 3, 5, 7, 11
    println("\n2. Creating witness with values:")
    println("   x = $x, y = $y, z = $z, u = $u")
    
    witness = create_witness_multiplication(x, y, z, u)
    r = x * y * z * u
    println("   Expected result: r = $r")
    
    # Verify the witness satisfies the R1CS
    if is_satisfied(r1cs, witness)
        println("   ✓ Witness satisfies R1CS constraints")
    else
        println("   ✗ Witness does NOT satisfy R1CS constraints")
        return
    end
    
    # Step 3: Convert R1CS to QAP
    println("\n3. Converting R1CS to QAP")
    qap = r1cs_to_qap(r1cs)
    println("   ✓ QAP created with degree-$(degree(qap.t)) target polynomial")
    println("   Evaluation domain: $(qap.domain)")
    
    # Step 4: Generate trusted setup
    println("\n4. Generating trusted setup (CRS)")
    crs = setup(qap)
    println("   ✓ Trusted setup complete")
    println("   - Generated $(length(crs.tau_powers_g1)) powers of τ in G1")
    println("   - Generated $(length(crs.tau_powers_g2)) powers of τ in G2")
    
    # Step 5: Generate proof
    println("\n5. Generating Groth16 proof")
    proof = prove(crs, qap, witness)
    println("   ✓ Proof generated:")
    println("   - [A]₁ ∈ G1")
    println("   - [B]₂ ∈ G2")
    println("   - [C]₁ ∈ G1")
    
    # Show proof elements in affine coordinates
    A_affine = to_affine(proof.A)
    C_affine = to_affine(proof.C)
    println("\n   Proof elements (affine coordinates):")
    println("   A.x = $(A_affine[1].value)")
    println("   A.y = $(A_affine[2].value)")
    println("   C.x = $(C_affine[1].value)")
    println("   C.y = $(C_affine[2].value)")
    
    # Step 6: Verify proof
    println("\n6. Verifying proof")
    
    # Extract public inputs (first 6 elements of witness: 1, r, x, y, z, u)
    public_inputs = witness.values[1:r1cs.num_public]
    
    # Verify using pairing check
    println("   Checking pairing equation: e(A, B) = e(C, [1]₂)")
    
    # For demonstration, we'll show the pairing computation
    # (Note: our pairing is simplified for demonstration)
    is_valid = verify(crs, proof, public_inputs)
    
    if is_valid
        println("   ✓ Proof is VALID!")
    else
        println("   ✗ Proof is INVALID!")
    end
    
    # Additional verification details
    println("\n7. Verification Summary")
    println("   Public inputs verified:")
    println("   - r = $(witness.values[2].value)")
    println("   - x = $(witness.values[3].value)")
    println("   - y = $(witness.values[4].value)")
    println("   - z = $(witness.values[5].value)")
    println("   - u = $(witness.values[6].value)")
    println("   Private witnesses (not revealed to verifier):")
    println("   - v1 = $(witness.values[7].value) (= x*y)")
    println("   - v2 = $(witness.values[8].value) (= z*u)")
    
    println("\n" * "=" ^ 60)
    println("Example completed successfully!")
    println("=" ^ 60)
end

# Run the example
main()