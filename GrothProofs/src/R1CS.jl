"""
R1CS (Rank-1 Constraint System) implementation.

An R1CS consists of three matrices L, R, O and a witness vector w such that:
``L \\cdot w \\circ R \\cdot w = O \\cdot w`` (where ``\\circ`` is element-wise multiplication).
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
    F = BN254Fr

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
    F = BN254Fr

    # Compute intermediate values
    v1 = x * y
    v2 = z * u
    r = v1 * v2

    # Create witness vector: [1, r, x, y, z, u, v1, v2]
    witness_values = [
        one(F),
        bn254_fr(r),
        bn254_fr(x),
        bn254_fr(y),
        bn254_fr(z),
        bn254_fr(u),
        bn254_fr(v1),
        bn254_fr(v2)
    ]

    return Witness{F}(witness_values)
end

# Export types and functions
export R1CS, Witness, is_satisfied
export create_r1cs_example_multiplication, create_witness_multiplication
export create_r1cs_example_sum_of_products, create_witness_sum_of_products
export create_r1cs_example_affine_product, create_witness_affine_product
export create_r1cs_example_square_offset, create_witness_square_offset

"""
    create_r1cs_example_sum_of_products()

Create the R1CS for r = x*y + z*u.
Variables: [1, r, x, y, z, u, v1, v2] where v1 = x*y, v2 = z*u, r = v1 + v2
"""
function create_r1cs_example_sum_of_products()
    F = BN254Fr

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
    L[1, 3] = one(F)  # x
    R[1, 4] = one(F)  # y
    O[1, 7] = one(F)  # v1

    # Constraint 2: v2 = z * u
    L[2, 5] = one(F)  # z
    R[2, 6] = one(F)  # u
    O[2, 8] = one(F)  # v2

    # Constraint 3: r = v1 + v2 -> (r - v1 - v2) * 1 = 0
    L[3, 2] = one(F)   # r
    L[3, 7] = -one(F)  # -v1
    L[3, 8] = -one(F)  # -v2
    R[3, 1] = one(F)   # multiply by 1
    # O row remains zero (equals 0)

    return R1CS{F}(num_vars, num_constraints, num_public, L, R, O)
end

"""
    create_witness_sum_of_products(x, y, z, u)

Create a witness for r = x*y + z*u.
"""
function create_witness_sum_of_products(x::Integer, y::Integer, z::Integer, u::Integer)
    F = BN254Fr

    v1 = x * y
    v2 = z * u
    r = v1 + v2

    witness_values = [
        one(F),
        bn254_fr(r),
        bn254_fr(x),
        bn254_fr(y),
        bn254_fr(z),
        bn254_fr(u),
        bn254_fr(v1),
        bn254_fr(v2),
    ]

    return Witness{F}(witness_values)
end

"""
    create_r1cs_example_affine_product()

Create the R1CS for r = (x + y) * (z + u).
Variables: [1, r, x, y, z, u, s1, s2] with helpers s1 = x + y, s2 = z + u.
"""
function create_r1cs_example_affine_product()
    F = BN254Fr

    num_vars = 8
    num_constraints = 3
    num_public = 6  # [1, r, x, y, z, u]

    L = zeros(F, num_constraints, num_vars)
    R = zeros(F, num_constraints, num_vars)
    O = zeros(F, num_constraints, num_vars)

    # Constraint 1: s1 = x + y -> (x + y) * 1 = s1
    L[1, 3] = one(F)  # x
    L[1, 4] = one(F)  # y
    R[1, 1] = one(F)  # multiply by constant 1
    O[1, 7] = one(F)  # s1

    # Constraint 2: s2 = z + u -> (z + u) * 1 = s2
    L[2, 5] = one(F)  # z
    L[2, 6] = one(F)  # u
    R[2, 1] = one(F)
    O[2, 8] = one(F)  # s2

    # Constraint 3: r = s1 * s2
    L[3, 7] = one(F)  # s1
    R[3, 8] = one(F)  # s2
    O[3, 2] = one(F)  # r

    return R1CS{F}(num_vars, num_constraints, num_public, L, R, O)
end

"""
    create_witness_affine_product(x, y, z, u)

Witness for r = (x + y) * (z + u).
"""
function create_witness_affine_product(x::Integer, y::Integer, z::Integer, u::Integer)
    F = BN254Fr

    s1 = x + y
    s2 = z + u
    r = s1 * s2

    witness_values = [
        one(F),
        bn254_fr(r),
        bn254_fr(x),
        bn254_fr(y),
        bn254_fr(z),
        bn254_fr(u),
        bn254_fr(s1),
        bn254_fr(s2),
    ]

    return Witness{F}(witness_values)
end

"""
    create_r1cs_example_square_offset()

Create the R1CS for r = x^2 + y + c with c public.
Variables: [1, r, x, y, c, x_sq].
"""
function create_r1cs_example_square_offset()
    F = BN254Fr

    num_vars = 6
    num_constraints = 2
    num_public = 5  # [1, r, x, y, c]

    L = zeros(F, num_constraints, num_vars)
    R = zeros(F, num_constraints, num_vars)
    O = zeros(F, num_constraints, num_vars)

    # Constraint 1: x_sq = x * x
    L[1, 3] = one(F)  # x
    R[1, 3] = one(F)  # x
    O[1, 6] = one(F)  # x_sq

    # Constraint 2: r = x_sq + y + c -> (x_sq + y + c) * 1 = r
    L[2, 6] = one(F)  # x_sq
    L[2, 4] = one(F)  # y
    L[2, 5] = one(F)  # c
    R[2, 1] = one(F)
    O[2, 2] = one(F)  # r

    return R1CS{F}(num_vars, num_constraints, num_public, L, R, O)
end

"""
    create_witness_square_offset(x, y, c)

Witness for r = x^2 + y + c.
"""
function create_witness_square_offset(x::Integer, y::Integer, c::Integer)
    F = BN254Fr

    x_sq = x * x
    r = x_sq + y + c

    witness_values = [
        one(F),
        bn254_fr(r),
        bn254_fr(x),
        bn254_fr(y),
        bn254_fr(c),
        bn254_fr(x_sq),
    ]

    return Witness{F}(witness_values)
end
