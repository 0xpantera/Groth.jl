# BN254 final exponentiation implementation.
#
# Raises the Miller loop output to `(p^12 - 1) / r` to map into the cyclotomic
# subgroup. The routine is split into easy and hard parts for efficiency.
#
# References:
# - Beuchat et al., "High-Speed Software Implementation of the Optimal Ate
#   Pairing over BN curves", 2010. https://eprint.iacr.org/2010/354.pdf
# - Devegili et al., "Implementing Cryptographic Pairings over Barreto-Naehrig
#   Curves", 2007. https://eprint.iacr.org/2007/390.pdf
# - LambdaClass, "How we implemented the BN254 Ate pairing in lambdaworks",
#   2024. https://blog.lambdaclass.com/how-we-implemented-the-bn254-ate-pairing-in-lambdaworks/
# - HackMD, "Computing the Optimal Ate Pairing over the BN254 Curve".

# using GrothAlgebra

# Import what we need
import GrothCurves: conjugate, BN254_PRIME

# ============================================================================
# Frobenius Constants for BN254 Sextic Twist
# ============================================================================

# ξ = 9 + u ∈ Fp2 (non-residue for the sextic twist)
# Rename to avoid redefinition warning with BN254MillerLoop.jl
const XI_FP2_FINALEXP = Fp2Element(9, 1)

# Precompute the gamma constants for Frobenius maps
# These arise from the sextic twist and tower construction
function compute_frobenius_constants()
    # Power for computing roots of unity
    pow = div(BN254_PRIME - 1, 6)

    # γ1[j] = ξ^((p-1)j/6) for j = 1..5
    gamma1 = Vector{Fp2Element}(undef, 5)
    for j in 1:5
        gamma1[j] = XI_FP2_FINALEXP^(pow * BigInt(j))
    end

    # γ2[j] = γ1[j] * conj(γ1[j])
    gamma2 = Vector{Fp2Element}(undef, 5)
    for j in 1:5
        gamma2[j] = gamma1[j] * conjugate(gamma1[j])
    end

    # γ3[j] = γ1[j] * γ2[j]
    gamma3 = Vector{Fp2Element}(undef, 5)
    for j in 1:5
        gamma3[j] = gamma1[j] * gamma2[j]
    end

    return gamma1, gamma2, gamma3
end

# Compute once at module load
const (GAMMA1, GAMMA2, GAMMA3) = compute_frobenius_constants()

# ============================================================================
# Frobenius Maps on Fp12
# ============================================================================

"""
    frobenius_p1(a::Fp12Element)

Apply the p-power Frobenius to an Fp12 element.
Uses precomputed gamma constants for the sextic twist.
"""
function frobenius_p1(a::Fp12Element)
    a0, a1 = a[1], a[2]  # Fp6 components

    # Extract Fp2 coefficients from each Fp6
    a00, a01, a02 = a0[1], a0[2], a0[3]
    a10, a11, a12 = a1[1], a1[2], a1[3]

    # Apply Frobenius to each Fp2 coefficient (conjugation)
    # and multiply by appropriate gamma factors
    c0 = Fp6Element(
        conjugate(a00),                    # v^0 term
        conjugate(a01) * GAMMA1[2],         # v^1 term
        conjugate(a02) * GAMMA1[4]          # v^2 term
    )

    c1 = Fp6Element(
        conjugate(a10) * GAMMA1[1],         # w * v^0 term
        conjugate(a11) * GAMMA1[3],         # w * v^1 term
        conjugate(a12) * GAMMA1[5]          # w * v^2 term
    )

    return Fp12Element(c0, c1)
end

"""
    frobenius_p2(a::Fp12Element)

Apply the p^2-power Frobenius to an Fp12 element.
"""
function frobenius_p2(a::Fp12Element)
    a0, a1 = a[1], a[2]

    a00, a01, a02 = a0[1], a0[2], a0[3]
    a10, a11, a12 = a1[1], a1[2], a1[3]

    # p^2 Frobenius: Fp2 elements unchanged (conjugate twice = identity)
    # Just apply gamma2 factors
    c0 = Fp6Element(
        a00,                               # v^0 term
        a01 * GAMMA2[2],                   # v^1 term
        a02 * GAMMA2[4]                    # v^2 term
    )

    c1 = Fp6Element(
        a10 * GAMMA2[1],                   # w * v^0 term
        a11 * GAMMA2[3],                   # w * v^1 term
        a12 * GAMMA2[5]                    # w * v^2 term
    )

    return Fp12Element(c0, c1)
end

"""
    frobenius_p3(a::Fp12Element)

Apply the p^3-power Frobenius to an Fp12 element.
"""
function frobenius_p3(a::Fp12Element)
    a0, a1 = a[1], a[2]

    a00, a01, a02 = a0[1], a0[2], a0[3]
    a10, a11, a12 = a1[1], a1[2], a1[3]

    # p^3 Frobenius: conjugate Fp2 elements and apply gamma3 factors
    c0 = Fp6Element(
        conjugate(a00),                    # v^0 term
        conjugate(a01) * GAMMA3[2],         # v^1 term
        conjugate(a02) * GAMMA3[4]          # v^2 term
    )

    c1 = Fp6Element(
        conjugate(a10) * GAMMA3[1],         # w * v^0 term
        conjugate(a11) * GAMMA3[3],         # w * v^1 term
        conjugate(a12) * GAMMA3[5]          # w * v^2 term
    )

    return Fp12Element(c0, c1)
end

"""
    frobenius_map(a::Fp12Element, power::Int)

Apply the p^power Frobenius endomorphism to an Fp12 element.
"""
function frobenius_map(a::Fp12Element, power::Int)
    p = mod(power, 12)  # Frobenius^12 = identity in Fp12

    if p == 0
        return a
    elseif p == 1
        return frobenius_p1(a)
    elseif p == 2
        return frobenius_p2(a)
    elseif p == 3
        return frobenius_p3(a)
    elseif p == 6
        # p^6 is conjugation in Fp12
        return conjugate(a)
    else
        # For other powers, compose
        result = a
        for _ in 1:p
            result = frobenius_p1(result)
        end
        return result
    end
end

# ============================================================================
# Final Exponentiation - Easy Part
# ============================================================================

"""
    final_exponentiation_easy(f::Fp12Element)

Compute the easy part of the final exponentiation: f^((p^6 - 1)(p^2 + 1))

This eliminates elements not in the cyclotomic subgroup and can be
computed efficiently using Frobenius maps and one inversion.
"""
function final_exponentiation_easy(f::Fp12Element)
    # Step 1: f^(p^6 - 1) = conj(f) / f
    f_conj = conjugate(f)
    f_inv = inv(f)
    f1 = f_conj * f_inv

    # Step 2: f1^(p^2 + 1) = f1^(p^2) * f1
    f1_p2 = frobenius_p2(f1)
    result = f1_p2 * f1

    return result
end

# ============================================================================
# Final Exponentiation - Hard Part
# ============================================================================

"""
    exp_by_u(a::Fp12Element)

Compute `a^u` where `u = 4965661367192848881` is the BN parameter.
"""
function exp_by_u(a::Fp12Element)
    result = one(Fp12Element)
    base = a
    u = UInt64(4965661367192848881)

    while u != 0
        if (u & UInt64(1)) == UInt64(1)
            result = result * base
        end
        base = base * base
        u >>= 1
    end

    return result
end

"""
    final_exponentiation_hard(m::Fp12Element)

Compute the hard part of the final exponentiation using the standard BN254 formula.

After the easy part, m is in the cyclotomic subgroup where inversion = conjugation.
The hard part computes m^((p^4 - p^2 + 1)/r) using a specific decomposition
involving u-powers and Frobenius maps.
"""
function final_exponentiation_hard(m::Fp12Element)
    # Compute u-powers using cyclotomic exponentiation.
    mx = exp_by_u(m)
    mx2 = exp_by_u(mx)
    mx3 = exp_by_u(mx2)

    # Compute Frobenius images.
    mp = frobenius_p1(m)
    mp2 = frobenius_p2(m)
    mp3 = frobenius_p3(m)

    mxp = frobenius_p1(mx)
    mx2p = frobenius_p1(mx2)
    mx3p = frobenius_p1(mx3)
    mx2p2 = frobenius_p2(mx2)

    y0 = mp * mp2 * mp3
    y1 = cyclotomic_inverse(m)
    y2 = mx2p2
    y3 = cyclotomic_inverse(mxp)
    y4 = cyclotomic_inverse(mx * mx2p)
    y5 = cyclotomic_inverse(mx2)
    y6 = cyclotomic_inverse(mx3 * mx3p)

    y1_2 = y1 * y1
    y2_2 = y2 * y2
    y2_6 = (y2_2 * y2_2) * y2_2
    y3_3 = (y3 * y3) * y3
    y3_6 = y3_3 * y3_3
    y3_12 = y3_6 * y3_6
    y4_2 = y4 * y4
    y4_4 = y4_2 * y4_2
    y4_8 = y4_4 * y4_4
    y4_16 = y4_8 * y4_8
    y4_18 = y4_16 * y4_2
    y5_2 = y5 * y5
    y5_4 = y5_2 * y5_2
    y5_8 = y5_4 * y5_4
    y5_16 = y5_8 * y5_8
    y5_30 = y5_16 * y5_8 * y5_4 * y5_2
    y6_2 = y6 * y6
    y6_4 = y6_2 * y6_2
    y6_8 = y6_4 * y6_4
    y6_16 = y6_8 * y6_8
    y6_32 = y6_16 * y6_16
    y6_36 = y6_32 * y6_4

    result = y0
    result = result * y1_2
    result = result * y2_6
    result = result * y3_12
    result = result * y4_18
    result = result * y5_30
    result = result * y6_36

    return result
end

# ============================================================================
# Complete Final Exponentiation
# ============================================================================

"""
    final_exponentiation(f::Fp12Element)

Complete final exponentiation: raise f to the power (p^12 - 1)/r.

This maps elements from the image of the Miller loop into the cyclotomic
subgroup GT, ensuring the pairing has the correct mathematical properties.
"""
function final_exponentiation(f::Fp12Element)
    # Check for special cases
    if iszero(f) || isone(f)
        return f
    end

    # Easy part: f^((p^6 - 1)(p^2 + 1))
    m = final_exponentiation_easy(f)

    # Hard part: m^((p^4 - p^2 + 1)/r)
    result = final_exponentiation_hard(m)

    return result
end

# Export functions
export frobenius_map, frobenius_p1, frobenius_p2, frobenius_p3
export final_exponentiation_easy, final_exponentiation_hard
export final_exponentiation, exp_by_u

"""
    final_exponentiation(::BN254Engine, f::Fp12Element)

Engine-aware overload that forwards to the BN254 final exponentiation.
"""
final_exponentiation(::BN254Engine, f::Fp12Element) = final_exponentiation(f)
