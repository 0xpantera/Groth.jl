module GrothCurves

using GrothAlgebra: BN254Fq, bn254_fq, AbstractCurve, GroupElem, prime
using StaticArrays: SVector

"""
    AbstractPairingEngine{Curve}

Interface tag for pairing backends tied to a particular curve family.
Concrete engines implement the generic pairing operations declared below.
"""
abstract type AbstractPairingEngine{Curve} end

# Pairing interface function declarations; concrete engines provide methods.
function miller_loop end
function final_exponentiation end
function pairing end
function pairing_batch end

# Include BN254 curve implementation
include("BN254Fp2.jl")  # Quadratic extension Fp2/Fp
include("BN254Fp6_3over2.jl")  # Cubic extension Fp6/Fp2
include("BN254Fp12_2over6.jl")  # Quadratic extension Fp12/Fp6
include("ProjectivePoint.jl")  # Generic projective point helpers
include("BN254Curve.jl")  # Curve operations G1, G2
include("BN254Engine.jl")  # Pairing engine wrapper
include("BN254MillerLoop.jl")  # Miller loop implementation
include("BN254FinalExp.jl")  # Final exponentiation
include("BN254Pairing.jl")  # Complete pairing

# Re-export from submodules
export AbstractPairingEngine
export BN254Fq, bn254_fq, BN254_PRIME
export Fp2Element, conjugate, norm, frobenius, real, imag
export Fp6Element, Fp12Element, GTElement, square
export BN254Curve, G1Point, G2Point, BN254Engine, BN254_ENGINE
export to_affine, is_on_curve, double, batch_to_affine!
export g1_generator, g2_generator
export x_coord, y_coord, z_coord
export LineCoeffs, doubling_step, addition_step, evaluate_line, miller_loop, frobenius_g2
export frobenius_map, frobenius_p1, frobenius_p2, frobenius_p3
export final_exponentiation_easy, final_exponentiation_hard
export final_exponentiation, exp_by_u
export optimal_ate_pairing, pairing, pairing_batch

end
