"""
GrothExamples - Educational examples and demonstrations for the Groth zkSNARK ecosystem

This package contains various examples demonstrating the use of:
- GrothAlgebra: Finite field arithmetic and polynomial operations
- GrothCurves: Elliptic curve operations
- GrothProofs: Zero-knowledge proof systems
"""
module GrothExamples

using GrothAlgebra
using GrothCurves
using GrothProofs

# Include example files
include("multiplication_proof.jl")
include("test_r1cs_qap.jl")

export demonstrate_r1cs_qap, multiplication_proof_example

end # module GrothExamples