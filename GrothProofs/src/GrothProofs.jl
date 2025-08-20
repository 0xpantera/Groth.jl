module GrothProofs

using GrothAlgebra
using GrothCurves

# Include proof system components
include("R1CS.jl")
include("QAP.jl")
# include("Groth16.jl")  # Temporarily disabled until curves are ready

# Re-export types and functions
export R1CS, Witness, is_satisfied
export create_r1cs_example_multiplication, create_witness_multiplication
export QAP, r1cs_to_qap, evaluate_qap, compute_h_polynomial
# export TrustedSetup, Groth16Proof
# export setup, prove, verify

end
