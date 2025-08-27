"""
BN254 pairing implementation.

This module implements the complete BN254 optimal ate pairing by combining
the Miller loop with final exponentiation.

The pairing e: G1 × G2 → GT satisfies bilinearity:
- e(aP, bQ) = e(P, Q)^(ab)
- e(P + R, Q) = e(P, Q) · e(R, Q)
- e(P, Q + S) = e(P, Q) · e(P, S)
"""

using GrothAlgebra

# Import necessary functions
import GrothCurves: miller_loop, final_exponentiation

# GT is the target group - cyclotomic subgroup of Fp12*
const GTElement = Fp12Element

"""
    optimal_ate_pairing(P::G1Point, Q::G2Point)

Compute the optimal ate pairing e(P, Q).

This performs:
1. Miller loop to compute f_{6u+2,Q}(P)
2. Final exponentiation to map to the cyclotomic subgroup

The result is an element in GT ⊂ Fp12* of order r.
"""
function optimal_ate_pairing(P::G1Point, Q::G2Point)
    # Handle special cases
    if iszero(P) || iszero(Q)
        return one(GTElement)
    end
    
    # Step 1: Miller loop
    f = miller_loop(P, Q)
    
    # Step 2: Final exponentiation
    result = final_exponentiation(f)
    
    return result
end

"""
    pairing(P::G1Point, Q::G2Point)

Alias for optimal_ate_pairing.
"""
pairing(P::G1Point, Q::G2Point) = optimal_ate_pairing(P, Q)

"""
    pairing_batch(P_vec::Vector{G1Point}, Q_vec::Vector{G2Point})

Compute the product of pairings: ∏ e(Pᵢ, Qᵢ)

This is more efficient than computing individual pairings and multiplying,
as we can share the final exponentiation.
"""
function pairing_batch(P_vec::Vector{G1Point}, Q_vec::Vector{G2Point})
    if length(P_vec) != length(Q_vec)
        throw(ArgumentError("Input vectors must have the same length"))
    end
    
    if isempty(P_vec)
        return one(GTElement)
    end
    
    # Accumulate Miller loop results
    f = one(Fp12Element)
    for (P, Q) in zip(P_vec, Q_vec)
        if !iszero(P) && !iszero(Q)
            f = f * miller_loop(P, Q)
        end
    end
    
    # Single final exponentiation
    return final_exponentiation(f)
end

# Export pairing functions
export GTElement, optimal_ate_pairing, pairing, pairing_batch