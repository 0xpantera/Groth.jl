module GrothCurves

using GrothAlgebra: BN254Field, bn254_field, AbstractCurve, GroupElem, prime
using StaticArrays: SVector

# Include BN254 curve implementation
include("BN254Fp2.jl")  # Quadratic extension Fp2/Fp
include("BN254Fp6_3over2.jl")  # Cubic extension Fp6/Fp2
include("BN254Fp12_2over6.jl")  # Quadratic extension Fp12/Fp6
include("BN254Curve.jl")  # Curve operations G1, G2
include("BN254MillerLoop.jl")  # Miller loop implementation
include("BN254FinalExp.jl")  # Final exponentiation
include("BN254Pairing.jl")  # Complete pairing

# Re-export from submodules
export BN254Field, bn254_field, BN254_PRIME
export Fp2Element, conjugate, norm, frobenius, real, imag
export Fp6Element, Fp12Element, GTElement, square
export BN254Curve, G1Point, G2Point
export to_affine, is_on_curve, double
export g1_generator, g2_generator
export x_coord, y_coord, z_coord
export LineCoeffs, doubling_step, addition_step, evaluate_line, miller_loop, frobenius_g2
export frobenius_map, frobenius_p1, frobenius_p2, frobenius_p3
export final_exponentiation_easy, final_exponentiation_hard
export final_exponentiation, exp_by_u
export optimal_ate_pairing, pairing, pairing_batch

end
