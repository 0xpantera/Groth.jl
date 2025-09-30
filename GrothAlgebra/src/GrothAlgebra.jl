module GrothAlgebra

# Import dependencies
using BitIntegers

# Include submodules
include("FiniteFields.jl")  # Active finite field implementation using BigInt

include("Group.jl")
include("Polynomial.jl")

# Export field element types and operations
export FiniteFieldElement, FieldElem
export BN254Field, bn254_field
export Secp256k1Field, secp256k1_field
export prime, is_zero, is_one, is_unity

# Export group element types and operations
export AbstractCurve, GroupElem
export scalar_mul, group_identity, inv, iszero, isone
export multi_scalar_mul, wnaf_encode, scalar_mul_wnaf
export FixedBaseTable, build_fixed_table, mul_fixed, batch_mul
export double, triple, order, is_on_curve

# Export polynomial operations
export Polynomial, degree, evaluate, interpolate
export leading_coefficient, is_constant, is_monic
export constant_polynomial, monomial, derivative
export roots_of_unity, primitive_root_of_unity, EvaluationDomain
export get_coset, coset_offset, coset_offset_inv, coset_offset_pow_size
export ntt!, fft, ifft, interpolate_fft, fft_polynomial_multiply

# Export utility functions
export isequal

end
