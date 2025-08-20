"""
R1CS (Rank-1 Constraint System) implementation.

An R1CS consists of three matrices L, R, O and a witness vector w such that:
L·w ∘ R·w = O·w (where ∘ is element-wise multiplication)
"""

using GrothCurves

"""
    R1CS{F}

Rank-1 Constraint System with field elements of type F.
"""
struct R1CS{F}
    # Number of variables (including the constant 1)
    num_vars::Int
    # Number of constraints
    num_constraints::Int
    # Number of public inputs (including the constant 1)
    num_public::Int
    # Constraint matrices
    L::Matrix{F}  # Left matrix
    R::Matrix{F}  # Right matrix
    O::Matrix{F}  # Output matrix
end

"""
    Witness{F}

Witness vector for an R1CS.
The first element is always 1, followed by public inputs, then private witnesses.
"""
struct Witness{F}
    values::Vector{F}
end

"""
    is_satisfied(r1cs::R1CS, witness::Witness)

Check if the witness satisfies all constraints in the R1CS.
"""
function is_satisfied(r1cs::R1CS{F}, witness::Witness{F}) where F
    w = witness.values
    
    # Check dimensions
    if length(w) != r1cs.num_vars
        error("Witness length $(length(w)) doesn't match number of variables $(r1cs.num_vars)")
    end
    
    # Check that first element is 1
    if !isone(w[1])
        error("First element of witness must be 1")
    end
    
    # Compute L·w, R·w, O·w
    Lw = r1cs.L * w
    Rw = r1cs.R * w
    Ow = r1cs.O * w
    
    # Check L·w ∘ R·w = O·w for each constraint
    for i in 1:r1cs.num_constraints
        if Lw[i] * Rw[i] != Ow[i]
            return false
        end
    end
    
    return true
end

"""
    create_r1cs_example_multiplication()

Create the R1CS for r = x * y * z * u as described in the example.
Variables: [1, r, x, y, z, u, v1, v2] where v1 = x*y, v2 = z*u, r = v1*v2
"""
function create_r1cs_example_multiplication()
    F = BN254Field
    
    # 8 variables: [1, r, x, y, z, u, v1, v2]
    num_vars = 8
    # 3 constraints
    num_constraints = 3
    # 6 public inputs: [1, r, x, y, z, u]
    num_public = 6
    
    # Initialize matrices with zeros
    L = zeros(F, num_constraints, num_vars)
    R = zeros(F, num_constraints, num_vars)
    O = zeros(F, num_constraints, num_vars)
    
    # Constraint 1: v1 = x * y
    # L: x (column 3), R: y (column 4), O: v1 (column 7)
    L[1, 3] = one(F)  # x
    R[1, 4] = one(F)  # y
    O[1, 7] = one(F)  # v1
    
    # Constraint 2: v2 = z * u
    # L: z (column 5), R: u (column 6), O: v2 (column 8)
    L[2, 5] = one(F)  # z
    R[2, 6] = one(F)  # u
    O[2, 8] = one(F)  # v2
    
    # Constraint 3: r = v1 * v2
    # L: v1 (column 7), R: v2 (column 8), O: r (column 2)
    L[3, 7] = one(F)  # v1
    R[3, 8] = one(F)  # v2
    O[3, 2] = one(F)  # r
    
    return R1CS{F}(num_vars, num_constraints, num_public, L, R, O)
end

"""
    create_witness_multiplication(x, y, z, u)

Create a witness for r = x * y * z * u.
"""
function create_witness_multiplication(x::Integer, y::Integer, z::Integer, u::Integer)
    F = BN254Field
    
    # Compute intermediate values
    v1 = x * y
    v2 = z * u
    r = v1 * v2
    
    # Create witness vector: [1, r, x, y, z, u, v1, v2]
    witness_values = [
        one(F),
        bn254_field(r),
        bn254_field(x),
        bn254_field(y),
        bn254_field(z),
        bn254_field(u),
        bn254_field(v1),
        bn254_field(v2)
    ]
    
    return Witness{F}(witness_values)
end

# Export types and functions
export R1CS, Witness, is_satisfied
export create_r1cs_example_multiplication, create_witness_multiplication