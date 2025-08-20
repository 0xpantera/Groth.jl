module GrothAlgebra

# Import dependencies
using BitIntegers

# Include submodules
include("FiniteFields.jl")  # Active finite field implementation using BigInt

# Note: FieldElement.jl is kept for reference but not used
# It contains the original parametric type implementation with UInt256
# which caused LLVM compilation issues

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
export double, triple, order, is_on_curve

# Export polynomial operations
export Polynomial, degree, evaluate, interpolate
export leading_coefficient, is_constant, is_monic
export constant_polynomial, monomial, derivative

# Export utility functions
export isequal

end
