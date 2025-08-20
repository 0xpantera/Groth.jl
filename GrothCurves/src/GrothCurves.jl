module GrothCurves

using GrothAlgebra
using StaticArrays

# Include BN254 curve implementation
include("BN254Fields.jl")
# include("BN254Curve.jl")  # Temporarily comment out
# include("BN254Pairing.jl")  # Temporarily comment out

# Re-export from submodules
export BN254Field, bn254_field, BN254_PRIME
export Fp2Element, conjugate, norm, frobenius
# export BN254Curve, G1Point, G2Point
# export to_affine, is_on_curve, double
# export g1_generator, g2_generator
# export optimal_ate_pairing, pairing

end
