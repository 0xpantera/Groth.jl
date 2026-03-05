# BN254 pairing implementation.
#
# Combines the Miller loop with final exponentiation to realise the optimal ate
# pairing e : G1 × G2 → GT.

# using GrothAlgebra

# Import necessary functions
import GrothCurves: miller_loop, final_exponentiation, pairing, pairing_batch

# GT is the target group - cyclotomic subgroup of Fp12*
const GTElement = Fp12Element

"""
    optimal_ate_pairing(P::G1Point, Q::G2Point)

Compute the optimal ate pairing ``e(P, Q)``.

This performs:
1. Miller loop to compute ``f_{6u+2,Q}(P)``
2. Final exponentiation to map to the cyclotomic subgroup

The result is an element in ``GT \\subset \\mathbb{F}_{p^{12}}^*`` of order ``r``.
"""
function optimal_ate_pairing(engine::BN254Engine, P::G1Point, Q::G2Point)
    # Handle special cases
    if iszero(P) || iszero(Q)
        return one(GTElement)
    end

    # Step 1: Miller loop
    f = miller_loop(engine, P, Q)

    # Step 2: Final exponentiation
    result = final_exponentiation(engine, f)

    return result
end

optimal_ate_pairing(P::G1Point, Q::G2Point) = optimal_ate_pairing(BN254_ENGINE, P, Q)

"""
    pairing(P::G1Point, Q::G2Point)

Call `optimal_ate_pairing`.
"""
pairing(engine::BN254Engine, P::G1Point, Q::G2Point) = optimal_ate_pairing(engine, P, Q)
pairing(P::G1Point, Q::G2Point) = pairing(BN254_ENGINE, P, Q)

"""
    pairing_batch(P_vec::Vector{G1Point}, Q_vec::Vector{G2Point})

Compute the product of pairings:

```math
\\prod_i e(P_i, Q_i)
```

This is more efficient than computing individual pairings and multiplying,
as we can share the final exponentiation.
"""
function pairing_batch(engine::BN254Engine, P_vec::Vector{G1Point}, Q_vec::Vector{G2Point})
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
            f = f * miller_loop(engine, P, Q)
        end
    end

    # Single final exponentiation
    return final_exponentiation(engine, f)
end

pairing_batch(P_vec::Vector{G1Point}, Q_vec::Vector{G2Point}) =
    pairing_batch(BN254_ENGINE, P_vec, Q_vec)

# Export pairing functions and engine
export GTElement, BN254Engine, BN254_ENGINE
export optimal_ate_pairing, pairing, pairing_batch
