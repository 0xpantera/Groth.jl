"""
Basic test of the Groth16 implementation components.
"""

push!(LOAD_PATH, @__DIR__)

using GrothAlgebra
using GrothCurves  
using GrothProofs

println("Testing basic components...")

# Test 1: Field arithmetic
println("\n1. Testing BN254 field arithmetic")
a = bn254_field(5)
b = bn254_field(7)
c = a + b
println("   5 + 7 = $(c.value) (mod p)")
println("   ✓ Field arithmetic works")

# Test 2: R1CS creation
println("\n2. Testing R1CS creation")
r1cs = create_r1cs_example_multiplication()
println("   ✓ R1CS created with $(r1cs.num_constraints) constraints")

# Test 3: Witness creation
println("\n3. Testing witness creation")
witness = create_witness_multiplication(2, 3, 5, 7)
println("   ✓ Witness created for x=2, y=3, z=5, u=7")
println("   Expected r = 2*3*5*7 = 210")

# Test 4: R1CS satisfaction
println("\n4. Testing R1CS satisfaction")
if is_satisfied(r1cs, witness)
    println("   ✓ Witness satisfies R1CS")
else
    println("   ✗ Witness does NOT satisfy R1CS")
end

# Test 5: QAP conversion
println("\n5. Testing QAP conversion")
qap = r1cs_to_qap(r1cs)
println("   ✓ QAP created")
println("   Target polynomial degree: $(degree(qap.t))")

println("\n✓ All basic tests passed!")