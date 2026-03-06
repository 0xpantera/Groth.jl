module GrothProofs

using GrothAlgebra
using GrothCurves

# Include proof system components
include("R1CS.jl")
include("QAP.jl")
include("Groth16.jl")

# Re-export types and functions
export R1CS, Witness, is_satisfied
export create_r1cs_example_multiplication, create_witness_multiplication
export QAP, r1cs_to_qap, evaluate_qap, compute_h_polynomial
export Groth16Proof
export ProvingKey, VerificationKey, Keypair
export setup_full, prove_full, verify_full
export PreparedVerificationKey, prepare_verifying_key, prepare_inputs, verify_with_prepared
export validate_witness_shape, public_inputs_from_witness
export setup, prove, process_vk, verify, verify_prepared

end
